//! The map panel's geometry (section 8): the pivot as a horizontal axis with its
//! contigs side by side, and where one member's differences sit on it.
//!
//! Everything here is arithmetic on 0-based half-open coordinates. The only
//! place a position becomes 1-based is [`Axis::label`], which is what the user
//! reads when hovering (section 5).

use infinity_core::events::MapUpdateInfo;
use infinity_core::interval::Iv;

use crate::layout;
use crate::state::UiModel;

/// The pivot's contigs laid end to end.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Axis {
    pub contigs: Vec<String>,
    pub lens: Vec<u64>,
    /// Start of each contig in flat coordinates.
    pub offsets: Vec<u64>,
    pub total: u64,
}

impl Axis {
    pub fn new(contigs: Vec<String>, lens: Vec<u64>) -> Self {
        let offsets = crate::pairwise::contig_offsets(&lens);
        let total = lens.iter().copied().sum();
        Self {
            contigs,
            lens,
            offsets,
            total,
        }
    }

    pub fn from_map(m: &MapUpdateInfo) -> Self {
        Self::new(m.pivot_contigs.clone(), m.pivot_contig_lens.clone())
    }

    /// Which contig a flat position falls in, and its 0-based offset inside it.
    pub fn contig_at(&self, flat: u64) -> Option<(usize, u64)> {
        if flat >= self.total {
            return None;
        }
        // Last contig whose offset is not past the position.
        let i = self
            .offsets
            .partition_point(|&o| o <= flat)
            .saturating_sub(1);
        if i >= self.lens.len() {
            return None;
        }
        let local = flat - self.offsets[i];
        if local < self.lens[i] {
            Some((i, local))
        } else {
            None
        }
    }

    /// A position as the user reads it: `contig:1-based`. A position past the
    /// end of the pivot has no label rather than a wrong one.
    pub fn label(&self, flat: u64) -> String {
        match self.contig_at(flat) {
            Some((i, local)) => format!("{}:{}", self.contigs[i], local + 1),
            None => "-".to_string(),
        }
    }

    /// A "nice" gridline step and the positions it produces, for the axis
    /// ruler. `target` is roughly how many gridlines should be visible.
    pub fn grid(&self, target: usize) -> Vec<u64> {
        let Some(step) = layout::nice_below(self.total as f64 / target.max(1) as f64) else {
            return Vec::new();
        };
        let step = step as u64;
        if step == 0 {
            return Vec::new();
        }
        let mut out = Vec::new();
        let mut p = step;
        while p < self.total && out.len() < 4096 {
            out.push(p);
            p += step;
        }
        out
    }
}

/// The complement of a merged interval list inside `[0, total)`: the stretches
/// that *are* covered. The map draws these as the thin grey line of section 8
/// and leaves the rest as a gap.
pub fn complement(uncovered: &[Iv], total: u64) -> Vec<Iv> {
    let mut out = Vec::new();
    let mut cursor = 0u64;
    for &(s, e) in uncovered {
        if s > cursor {
            out.push((cursor, s.min(total)));
        }
        cursor = cursor.max(e);
        if cursor >= total {
            break;
        }
    }
    if cursor < total {
        out.push((cursor, total));
    }
    out
}

/// Why the map is showing the rows it is showing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RowReason {
    /// The genomes the user selected first (restricted to this cluster), then
    /// the rest of the cluster.
    Selection,
    /// Nothing selected, so the cluster's members, up to a cap.
    WholeCluster,
    /// No cluster chosen, or the cluster has no map data yet.
    #[default]
    NoData,
}

/// Which members get a row.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RowPlan {
    pub genomes: Vec<u32>,
    pub reason: RowReason,
    /// How many of the first rows are selected genomes.
    pub selected: usize,
    /// True when rows had to be dropped to stay drawable.
    pub truncated: bool,
    pub available: usize,
}

/// Members of `cluster` that have map data: the selected ones first, in the
/// order they were selected, then the others in id order. Capped so that a
/// cluster of a thousand genomes cannot ask for a thousand rows in a 180 px
/// panel; selected rows are never dropped by the cap.
pub fn plan_rows(model: &UiModel, cluster: Option<u32>, selection: &[u32], cap: usize) -> RowPlan {
    let Some(cid) = cluster else {
        return RowPlan::default();
    };
    let mapped = model.mapped_members(cid);
    if mapped.is_empty() {
        return RowPlan::default();
    }
    let available = mapped.len();
    let mut genomes: Vec<u32> = selection
        .iter()
        .filter(|g| mapped.binary_search(g).is_ok())
        .copied()
        .collect();
    let selected = genomes.len();
    let limit = cap.max(selected);
    for g in mapped {
        if genomes.len() >= limit {
            break;
        }
        if !genomes[..selected].contains(&g) {
            genomes.push(g);
        }
    }
    RowPlan {
        truncated: genomes.len() < available,
        genomes,
        reason: if selected > 0 {
            RowReason::Selection
        } else {
            RowReason::WholeCluster
        },
        selected,
        available,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::PieceView;
    use infinity_core::interval::merge;
    use infinity_core::variant::{Variant, VariantKind};

    fn axis3() -> Axis {
        Axis::new(
            vec!["c0".to_string(), "c1".to_string(), "c2".to_string()],
            vec![100, 50, 25],
        )
    }

    fn member(cluster: u32, genome: u32, uncovered: Vec<Iv>) -> MapUpdateInfo {
        MapUpdateInfo {
            cluster,
            genome,
            genome_name: format!("g{genome}"),
            pivot_name: "p".to_string(),
            pieces: vec![PieceView {
                pivot_contig: 0,
                pivot_contig_name: "c0".to_string(),
                start: 10,
                end: 20,
                resolved: true,
            }],
            variants: vec![Variant {
                kind: VariantKind::Snp,
                pivot_contig: 0,
                pivot_start: 12,
                pivot_end: 13,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                genome_contig: 0,
                genome_start: 12,
                genome_end: 13,
                strand: '+',
                flags: 0,
            }],
            uncovered,
            pivot_contigs: vec!["c0".into(), "c1".into(), "c2".into()],
            pivot_contig_lens: vec![100, 50, 25],
        }
    }

    #[test]
    fn the_axis_lays_contigs_end_to_end() {
        let a = axis3();
        assert_eq!(a.offsets, vec![0, 100, 150]);
        assert_eq!(a.total, 175);
        assert_eq!(Axis::new(vec![], vec![]).total, 0);
    }

    #[test]
    fn a_flat_position_finds_its_contig_and_offset() {
        let a = axis3();
        assert_eq!(a.contig_at(0), Some((0, 0)));
        assert_eq!(a.contig_at(99), Some((0, 99)));
        assert_eq!(a.contig_at(100), Some((1, 0)));
        assert_eq!(a.contig_at(149), Some((1, 49)));
        assert_eq!(a.contig_at(150), Some((2, 0)));
        assert_eq!(a.contig_at(174), Some((2, 24)));
        assert_eq!(a.contig_at(175), None, "past the end of the pivot");
    }

    #[test]
    fn a_zero_length_contig_is_never_reported() {
        let a = Axis::new(vec!["a".into(), "b".into()], vec![0, 10]);
        assert_eq!(a.contig_at(0), Some((1, 0)));
    }

    #[test]
    fn hover_positions_are_one_based_and_named() {
        let a = axis3();
        assert_eq!(a.label(0), "c0:1");
        assert_eq!(a.label(99), "c0:100");
        assert_eq!(a.label(100), "c1:1");
        assert_eq!(a.label(174), "c2:25");
        assert_eq!(a.label(175), "-");
    }

    #[test]
    fn gridlines_are_evenly_spaced_and_inside_the_pivot() {
        let a = axis3();
        let g = a.grid(4);
        assert!(!g.is_empty());
        assert!(g.iter().all(|&p| p > 0 && p < a.total));
        let step = g[1] - g[0];
        assert!(g.windows(2).all(|w| w[1] - w[0] == step));
        assert!(Axis::new(vec![], vec![0]).grid(4).is_empty());
    }

    #[test]
    fn the_complement_of_the_gaps_is_what_gets_a_grey_line() {
        assert_eq!(complement(&[], 100), vec![(0, 100)]);
        assert_eq!(complement(&[(0, 100)], 100), Vec::<Iv>::new());
        assert_eq!(complement(&[(10, 20)], 100), vec![(0, 10), (20, 100)]);
        assert_eq!(
            complement(&[(0, 10), (90, 100)], 100),
            vec![(10, 90)],
            "gaps at both ends leave one covered stretch"
        );
        // Intervals are merged before use, so a gap past the end is clipped.
        assert_eq!(complement(&[(50, 200)], 100), vec![(0, 50)]);
    }

    #[test]
    fn the_complement_and_the_gaps_partition_the_pivot() {
        let unc = merge(vec![(10, 20), (19, 30), (90, 95)]);
        let cov = complement(&unc, 100);
        let total: u64 = unc.iter().map(|(s, e)| e - s).sum::<u64>()
            + cov.iter().map(|(s, e)| e - s).sum::<u64>();
        assert_eq!(total, 100);
    }

    #[test]
    fn selected_members_come_first_then_the_rest_of_the_cluster() {
        let mut m = UiModel::default();
        for g in 0..6 {
            m.apply(0.0, crate::state::scan_event(g, &format!("g{g}"), 100, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2, 3]), (1, 4, &[4, 5])]),
        );
        for g in 0..4 {
            m.apply(
                0.0,
                infinity_core::events::Event::MapUpdate(Box::new(member(0, g, vec![]))),
            );
        }
        // Genome 5 is selected too, but belongs to another cluster.
        let plan = plan_rows(&m, Some(0), &[3, 5, 1], 10);
        assert_eq!(plan.genomes, vec![3, 1, 0, 2]);
        assert_eq!(plan.selected, 2);
        assert_eq!(plan.reason, RowReason::Selection);
        assert_eq!(plan.available, 4);
        assert!(!plan.truncated);
        // The cap drops unselected rows, never selected ones.
        let plan = plan_rows(&m, Some(0), &[3, 1, 2], 2);
        assert_eq!(plan.genomes, vec![3, 1, 2]);
        assert!(plan.truncated);
    }

    #[test]
    fn with_no_selection_the_whole_cluster_is_shown_up_to_a_cap() {
        let mut m = UiModel::default();
        for g in 0..6 {
            m.apply(0.0, crate::state::scan_event(g, &format!("g{g}"), 100, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2, 3, 4, 5])]),
        );
        for g in 0..6 {
            m.apply(
                0.0,
                infinity_core::events::Event::MapUpdate(Box::new(member(0, g, vec![]))),
            );
        }
        let plan = plan_rows(&m, Some(0), &[], 4);
        assert_eq!(plan.genomes, vec![0, 1, 2, 3]);
        assert!(plan.truncated);
        assert_eq!(plan.available, 6);
        assert_eq!(plan.reason, RowReason::WholeCluster);
    }

    #[test]
    fn a_selection_of_genomes_from_another_cluster_shows_nothing() {
        let mut m = UiModel::default();
        for g in 0..2 {
            m.apply(0.0, crate::state::scan_event(g, &format!("g{g}"), 100, 1));
        }
        m.apply(0.0, crate::state::clustered_event(&[(0, 0, &[0, 1])]));
        assert_eq!(plan_rows(&m, Some(0), &[], 10).reason, RowReason::NoData);
        assert_eq!(plan_rows(&m, None, &[0], 10).reason, RowReason::NoData);
    }
}
