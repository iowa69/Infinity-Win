//! Which genomes are selected, and what a click does to that.
//!
//! The usual desktop rules, the same in the genome list, the tree and the map:
//! a click selects one genome, Ctrl-click (Cmd on a Mac) adds or removes one,
//! Shift-click selects the run of genomes between the last one clicked and this
//! one, in the order the clicked view shows them, and Ctrl+A selects every
//! genome the list shows. There is no limit on how many.
//!
//! Pure data, no egui: the rules are tested here, and the panels only translate
//! pointer modifiers into [`ClickMods`] and pass their own row order along.

use std::collections::{HashMap, HashSet};

use infinity_core::engine::RunResults;

use crate::fmt;
use crate::state::UiModel;

/// The modifier keys held during a click.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClickMods {
    /// Ctrl (Cmd on a Mac): add or remove one genome.
    pub toggle: bool,
    /// Shift: select the run from the anchor to here.
    pub range: bool,
}

impl ClickMods {
    pub const PLAIN: ClickMods = ClickMods {
        toggle: false,
        range: false,
    };
    pub const TOGGLE: ClickMods = ClickMods {
        toggle: true,
        range: false,
    };
    #[cfg(test)]
    pub const RANGE: ClickMods = ClickMods {
        toggle: false,
        range: true,
    };

    pub fn from_modifiers(m: &egui::Modifiers) -> Self {
        Self {
            toggle: m.command,
            range: m.shift,
        }
    }
}

/// The selected genomes, in the order they were chosen.
#[derive(Debug, Clone, Default)]
pub struct Selection {
    ids: Vec<u32>,
    set: HashSet<u32>,
    /// Where the next Shift-click range starts: the last genome clicked
    /// without Shift.
    anchor: Option<u32>,
}

impl PartialEq<Vec<u32>> for Selection {
    fn eq(&self, other: &Vec<u32>) -> bool {
        &self.ids == other
    }
}

impl Selection {
    pub fn ids(&self) -> &[u32] {
        &self.ids
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    pub fn contains(&self, id: u32) -> bool {
        self.set.contains(&id)
    }

    #[cfg(test)]
    pub fn anchor(&self) -> Option<u32> {
        self.anchor
    }

    pub fn clear(&mut self) {
        self.ids.clear();
        self.set.clear();
        self.anchor = None;
    }

    fn add(&mut self, id: u32) {
        if self.set.insert(id) {
            self.ids.push(id);
        }
    }

    fn remove(&mut self, id: u32) {
        if self.set.remove(&id) {
            self.ids.retain(|&g| g != id);
        }
    }

    /// Exactly this genome.
    pub fn set_one(&mut self, id: u32) {
        self.clear();
        self.add(id);
        self.anchor = Some(id);
    }

    /// Add or remove one genome; it becomes the anchor either way.
    pub fn toggle(&mut self, id: u32) {
        if self.contains(id) {
            self.remove(id);
        } else {
            self.add(id);
        }
        self.anchor = Some(id);
    }

    /// Every genome of `order`, replacing the selection. The anchor moves to
    /// the first of them, so a following Shift-click still makes sense.
    pub fn select_all(&mut self, order: &[u32]) {
        self.clear();
        for &id in order {
            self.add(id);
        }
        self.anchor = order.first().copied();
    }

    /// Select many at once, keeping what is selected already.
    pub fn extend(&mut self, ids: impl IntoIterator<Item = u32>) {
        for id in ids {
            self.add(id);
        }
    }

    /// Apply a click on `id`. `order` is the row order of the view that was
    /// clicked; it is only needed for a Shift-click, and a range whose anchor
    /// is not in it falls back to a plain click.
    pub fn click(&mut self, id: u32, mods: ClickMods, order: &[u32]) {
        if mods.range {
            let from = self.anchor.and_then(|a| order.iter().position(|&g| g == a));
            let to = order.iter().position(|&g| g == id);
            if let (Some(from), Some(to)) = (from, to) {
                let (lo, hi) = if from <= to { (from, to) } else { (to, from) };
                let anchor = self.anchor;
                if !mods.toggle {
                    self.clear();
                }
                for &g in &order[lo..=hi] {
                    self.add(g);
                }
                // The anchor stays put, so Shift-clicking again re-shapes the
                // same range, as in every file manager.
                self.anchor = anchor;
                return;
            }
        }
        if mods.toggle {
            self.toggle(id);
        } else {
            self.set_one(id);
        }
    }
}

/// The numbers shown when more than two genomes are selected.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Summary {
    pub n: usize,
    /// `(cluster, selected members)`, largest share first, then by id.
    pub clusters: Vec<(u32, usize)>,
    /// Selected genomes not in any cluster yet.
    pub unclustered: usize,
    /// SNP range over the selected pairs that share a cluster, when any could
    /// be measured.
    pub snps: Option<SnpRange>,
    /// Pairs that share a cluster.
    pub pairs: usize,
    /// Of those, the ones with no distance yet (below Focus 1, or a cluster
    /// too large for its matrix to travel in the event stream).
    pub missing: usize,
}

/// The smallest and largest SNP distance among measured pairs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SnpRange {
    pub min: f64,
    pub max: f64,
    /// At least one of the numbers is an estimate (Focus 1).
    pub approximate: bool,
}

impl Summary {
    /// One line for the tree header.
    pub fn text(&self) -> String {
        let mut parts = vec![format!("{} genomes selected", fmt::count(self.n as u64))];
        if !self.clusters.is_empty() {
            let shown: Vec<String> = self
                .clusters
                .iter()
                .take(6)
                .map(|(c, k)| format!("c{c} ×{k}"))
                .collect();
            let more = self.clusters.len().saturating_sub(6);
            parts.push(format!(
                "{} cluster{}: {}{}",
                self.clusters.len(),
                if self.clusters.len() == 1 { "" } else { "s" },
                shown.join(", "),
                if more > 0 {
                    format!(" +{more} more")
                } else {
                    String::new()
                }
            ));
        }
        if self.unclustered > 0 {
            parts.push(format!("{} not clustered yet", self.unclustered));
        }
        match (&self.snps, self.pairs) {
            (_, 0) => parts.push("no two share a cluster: no exact distance".to_string()),
            (Some(r), pairs) => {
                let m = fmt::approx_mark(r.approximate);
                let scope = if self.clusters.len() == 1 {
                    "between them".to_string()
                } else {
                    format!("within clusters ({} pairs)", fmt::count(pairs as u64))
                };
                let range = if (r.max - r.min).abs() < 0.5 {
                    format!("{m}{}", fmt::snps(r.min))
                } else {
                    format!("{m}{}–{}", fmt::snps(r.min), fmt::snps(r.max))
                };
                let missing = if self.missing > 0 {
                    format!(", {} not measured yet", fmt::count(self.missing as u64))
                } else {
                    String::new()
                };
                parts.push(format!("{range} SNPs {scope}{missing}"));
            }
            (None, _) => parts.push("SNP distances not measured yet".to_string()),
        }
        parts.join(" · ")
    }
}

/// Summarise a selection: clusters involved and the SNP range among pairs
/// that share a cluster. Distances come from the finished run's matrices when
/// there are any, else from the matrices the engine streamed per cluster
/// (`ClusterDetail`); the two hold the same numbers (SNPs after the SNP rules,
/// or Focus 1 estimates).
pub fn summarise(model: &UiModel, results: Option<&RunResults>, ids: &[u32]) -> Summary {
    let mut by_cluster: HashMap<u32, Vec<u32>> = HashMap::new();
    let mut unclustered = 0;
    for &id in ids {
        match model.cluster_of(id) {
            Some(c) => by_cluster.entry(c).or_default().push(id),
            None => unclustered += 1,
        }
    }
    let mut clusters: Vec<(u32, usize)> = by_cluster.iter().map(|(&c, v)| (c, v.len())).collect();
    clusters.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

    let mut out = Summary {
        n: ids.len(),
        clusters,
        unclustered,
        ..Default::default()
    };
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let mut approximate = false;
    for (&cid, members) in &by_cluster {
        if members.len() < 2 {
            continue;
        }
        let k = members.len();
        out.pairs += k * (k - 1) / 2;
        let lookup = matrix_lookup(model, results, cid, members);
        let Some((rows, values, approx)) = lookup else {
            out.missing += k * (k - 1) / 2;
            continue;
        };
        for i in 0..k {
            for j in (i + 1)..k {
                match (rows[i], rows[j]) {
                    (Some(a), Some(b)) => {
                        let v = values(a, b);
                        lo = lo.min(v);
                        hi = hi.max(v);
                        approximate |= approx;
                    }
                    _ => out.missing += 1,
                }
            }
        }
    }
    if lo.is_finite() && hi.is_finite() {
        out.snps = Some(SnpRange {
            min: lo,
            max: hi,
            approximate,
        });
    }
    out
}

/// Row lookup into whichever matrix exists for `cid`, as `(row per member,
/// value(row, row), approximate)`.
#[allow(clippy::type_complexity)]
fn matrix_lookup<'a>(
    model: &'a UiModel,
    results: Option<&'a RunResults>,
    cid: u32,
    members: &[u32],
) -> Option<(
    Vec<Option<usize>>,
    Box<dyn Fn(usize, usize) -> f64 + 'a>,
    bool,
)> {
    if let Some(m) = results.and_then(|r| r.matrices.get(&cid)) {
        let rows = members
            .iter()
            .map(|g| m.ids.iter().position(|x| x == g))
            .collect();
        let exact = matches!(m.basis, infinity_core::profile::TreeBasis::Snps);
        return Some((
            rows,
            Box::new(move |a, b| {
                if exact {
                    m.snps[a][b] as f64
                } else {
                    m.estimates
                        .get(a)
                        .and_then(|r| r.get(b))
                        .map(|&v| v as f64)
                        .unwrap_or(m.snps[a][b] as f64)
                }
            }),
            !exact,
        ));
    }
    let d = model.details.get(&cid)?;
    let matrix = d.matrix.as_ref()?;
    let index: HashMap<&str, usize> = d
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let rows = members
        .iter()
        .map(|&g| {
            model
                .genome(g)
                .and_then(|row| index.get(row.name.as_str()).copied())
        })
        .collect();
    Some((
        rows,
        Box::new(move |a, b| matrix[a][b]),
        d.focus < fmt::FOCUS_EXACT,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::{ClusterDetailInfo, Event};

    fn order() -> Vec<u32> {
        vec![5, 3, 9, 1, 7]
    }

    #[test]
    fn a_plain_click_selects_one_and_moves_the_anchor() {
        let mut s = Selection::default();
        s.click(3, ClickMods::PLAIN, &order());
        s.click(9, ClickMods::PLAIN, &order());
        assert_eq!(s, vec![9]);
        assert_eq!(s.anchor(), Some(9));
    }

    #[test]
    fn ctrl_click_adds_and_removes_without_any_limit() {
        let mut s = Selection::default();
        for id in [5, 3, 9, 1, 7] {
            s.click(id, ClickMods::TOGGLE, &order());
        }
        assert_eq!(s, vec![5, 3, 9, 1, 7], "no cap on the count");
        s.click(9, ClickMods::TOGGLE, &order());
        assert_eq!(s, vec![5, 3, 1, 7]);
        assert!(!s.contains(9));
        assert_eq!(s.len(), 4);
    }

    #[test]
    fn shift_click_selects_the_run_in_the_views_order_either_way() {
        let mut s = Selection::default();
        s.click(3, ClickMods::PLAIN, &order());
        s.click(1, ClickMods::RANGE, &order());
        assert_eq!(s, vec![3, 9, 1]);
        // From the same anchor, backwards: the range is re-shaped.
        s.click(5, ClickMods::RANGE, &order());
        assert_eq!(s, vec![5, 3]);
        assert_eq!(s.anchor(), Some(3), "the anchor stays");
    }

    #[test]
    fn ctrl_shift_click_adds_a_range_to_what_is_selected() {
        let mut s = Selection::default();
        s.click(7, ClickMods::PLAIN, &order());
        s.click(5, ClickMods::TOGGLE, &order());
        s.click(
            9,
            ClickMods {
                toggle: true,
                range: true,
            },
            &order(),
        );
        assert_eq!(s, vec![7, 5, 3, 9]);
    }

    #[test]
    fn a_range_without_an_anchor_in_view_is_a_plain_click() {
        let mut s = Selection::default();
        s.click(1, ClickMods::RANGE, &order());
        assert_eq!(s, vec![1]);
        s.click(42, ClickMods::PLAIN, &[]);
        s.click(7, ClickMods::RANGE, &order());
        assert_eq!(s, vec![7], "anchor 42 is not in this view");
    }

    #[test]
    fn select_all_and_clear() {
        let mut s = Selection::default();
        s.select_all(&order());
        assert_eq!(s.len(), 5);
        assert!(order().iter().all(|&g| s.contains(g)));
        s.clear();
        assert!(s.is_empty());
        assert_eq!(s.anchor(), None);
    }

    fn model() -> UiModel {
        let mut m = UiModel::default();
        for i in 0..6u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2, 3]), (1, 4, &[4, 5])]),
        );
        m
    }

    fn detail(cluster: u32, names: &[&str], matrix: Vec<Vec<f64>>, focus: u8) -> Event {
        Event::ClusterDetail(Box::new(ClusterDetailInfo {
            cluster,
            focus,
            ids: (0..names.len() as u32).collect(),
            names: names.iter().map(|s| s.to_string()).collect(),
            matrix: Some(matrix),
            raw_matrix: None,
            core_bp: 0,
            pass_snps: vec![],
            masked_by_rule: vec![],
            snp_thresholds: vec![],
            snp_groups: vec![],
        }))
    }

    #[test]
    fn the_summary_names_clusters_and_the_snp_range_of_pairs_sharing_one() {
        let mut m = model();
        let mat = vec![
            vec![0.0, 3.0, 8.0, 12.0],
            vec![3.0, 0.0, 5.0, 9.0],
            vec![8.0, 5.0, 0.0, 41.0],
            vec![12.0, 9.0, 41.0, 0.0],
        ];
        m.apply(0.0, detail(0, &["g0", "g1", "g2", "g3"], mat, 2));
        let s = summarise(&m, None, &[1, 2, 3, 4]);
        assert_eq!(s.n, 4);
        assert_eq!(s.clusters, vec![(0, 3), (1, 1)]);
        assert_eq!(s.pairs, 3);
        assert_eq!(s.missing, 0);
        let r = s.snps.unwrap();
        assert_eq!((r.min, r.max), (5.0, 41.0));
        assert!(!r.approximate);
        let t = s.text();
        assert!(t.starts_with("4 genomes selected"), "{t}");
        assert!(t.contains("c0 ×3, c1 ×1"), "{t}");
        assert!(t.contains("5–41 SNPs within clusters (3 pairs)"), "{t}");
    }

    #[test]
    fn the_summary_says_when_nothing_can_be_measured() {
        let m = model();
        let s = summarise(&m, None, &[0, 4]);
        assert_eq!(s.pairs, 0);
        assert!(s.text().contains("no two share a cluster"));
        let s = summarise(&m, None, &[0, 1, 2]);
        assert_eq!((s.pairs, s.missing), (3, 3));
        assert!(s.snps.is_none());
        assert!(s.text().contains("not measured yet"));
    }

    #[test]
    fn focus_one_matrices_are_marked_as_estimates() {
        let mut m = model();
        m.apply(
            0.0,
            detail(1, &["g4", "g5"], vec![vec![0.0, 7.0], vec![7.0, 0.0]], 1),
        );
        let s = summarise(&m, None, &[4, 5, 0]);
        let r = s.snps.unwrap();
        assert!(r.approximate);
        assert!(s.text().contains("≈ 7 SNPs"), "{}", s.text());
    }
}
