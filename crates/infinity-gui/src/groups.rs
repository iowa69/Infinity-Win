//! SNP groups: inside each exact cluster, the genomes linked by chains of at
//! most N SNPs (single linkage), recomputed live as the slider in the tree
//! header moves.
//!
//! The pairwise SNPs come from the engine. A `ClusterDetail` event carries the
//! filtered matrix of a cluster at Focus 2 for clusters of up to
//! [`infinity_core::events::DETAIL_MATRIX_MAX`] genomes; above that, and for a
//! finished run whose recording predates the event, the run folder's
//! `cluster_<id>.snp_matrix.tsv` is read instead (a square table with genome
//! names as row and column headers, `NA` where a pair shares no bases).
//!
//! Single linkage at N is exactly what the engine reports in its
//! `snp_clusters.tsv` for its own thresholds (5, 10, 20): two genomes are in
//! one group when a chain of genomes joins them in which every step is at most
//! N SNPs. A group is therefore not "everything within N of everything", and
//! the side panel shows each group's widest pair (its diameter) next to N.
//!
//! Pure data, no egui: everything here is tested without a window.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};

use infinity_core::events::ClusterDetailInfo;

use crate::layout::{DisplayNode, NodeKind};
use crate::state::UiModel;

/// A square SNP matrix of one cluster, in the model's genome ids.
#[derive(Debug, Clone, PartialEq)]
pub struct SnpMatrix {
    pub cluster: u32,
    /// Genome ids in matrix order.
    pub ids: Vec<u32>,
    /// Row-major, `n * n`; NaN where the pair was not comparable.
    d: Vec<f32>,
}

impl SnpMatrix {
    /// A matrix from rows of values (`None` = not comparable). Fails unless
    /// it is square and matches the ids.
    pub fn new(cluster: u32, ids: Vec<u32>, rows: &[Vec<Option<f64>>]) -> Result<Self, String> {
        let n = ids.len();
        if rows.len() != n || rows.iter().any(|r| r.len() != n) {
            return Err(format!(
                "cluster {cluster}: the SNP matrix is not {n} × {n}"
            ));
        }
        let d = rows
            .iter()
            .flat_map(|r| {
                r.iter().map(|v| match v {
                    Some(x) if x.is_finite() => *x as f32,
                    _ => f32::NAN,
                })
            })
            .collect();
        Ok(Self { cluster, ids, d })
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// SNPs between rows `i` and `j`, `None` when not comparable.
    pub fn get(&self, i: usize, j: usize) -> Option<f64> {
        let v = self.d[i * self.ids.len() + j];
        (!v.is_nan()).then_some(f64::from(v))
    }
}

/// The matrix a `ClusterDetail` event carries, in model ids. `None` unless it
/// is an exact (Focus 2) matrix, it is present (small clusters only), and
/// every genome name resolves.
pub fn from_detail(d: &ClusterDetailInfo, ids_by_name: &HashMap<String, u32>) -> Option<SnpMatrix> {
    if d.focus < crate::fmt::FOCUS_EXACT {
        return None;
    }
    let m = d.matrix.as_ref()?;
    let ids: Option<Vec<u32>> = d
        .names
        .iter()
        .map(|n| ids_by_name.get(n).copied())
        .collect();
    let rows: Vec<Vec<Option<f64>>> = m
        .iter()
        .map(|r| r.iter().map(|&v| Some(v)).collect())
        .collect();
    SnpMatrix::new(d.cluster, ids?, &rows).ok()
}

/// Parse the run folder's `cluster_<id>.snp_matrix.tsv`.
pub fn parse_tsv(
    text: &str,
    cluster: u32,
    ids_by_name: &HashMap<String, u32>,
) -> Result<SnpMatrix, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines
        .next()
        .ok_or_else(|| format!("cluster {cluster}: empty SNP matrix"))?;
    let names: Vec<&str> = header.split('\t').skip(1).collect();
    let resolve = |n: &str| {
        ids_by_name
            .get(n)
            .copied()
            .ok_or_else(|| format!("cluster {cluster}: unknown genome {n} in the SNP matrix"))
    };
    let ids: Vec<u32> = names.iter().map(|n| resolve(n)).collect::<Result<_, _>>()?;
    let mut rows: Vec<Vec<Option<f64>>> = Vec::with_capacity(ids.len());
    for (k, line) in lines.enumerate() {
        let mut cells = line.split('\t');
        let row_name = cells.next().unwrap_or_default();
        if names.get(k) != Some(&row_name) {
            return Err(format!(
                "cluster {cluster}: row {} is {row_name}, the header says {}",
                k + 1,
                names.get(k).unwrap_or(&"nothing")
            ));
        }
        rows.push(
            cells
                .map(|c| match c.trim() {
                    "NA" | "" => None,
                    v => v.parse::<f64>().ok(),
                })
                .collect(),
        );
    }
    SnpMatrix::new(cluster, ids, &rows)
}

/// The run folder's matrix file for a cluster.
pub fn tsv_path(dir: &Path, cluster: u32) -> std::path::PathBuf {
    dir.join(format!("cluster_{cluster}.snp_matrix.tsv"))
}

/// Read and parse a cluster's matrix from the run folder.
pub fn read_tsv(
    dir: &Path,
    cluster: u32,
    ids_by_name: &HashMap<String, u32>,
) -> Result<SnpMatrix, String> {
    let p = tsv_path(dir, cluster);
    let text = std::fs::read_to_string(&p).map_err(|e| format!("{}: {e}", p.display()))?;
    parse_tsv(&text, cluster, ids_by_name)
}

/// Single-linkage components at `threshold`: every row index in exactly one
/// component, each component sorted, components ordered by their first row.
pub fn components(m: &SnpMatrix, threshold: f64) -> Vec<Vec<usize>> {
    let n = m.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    let limit = threshold + 1e-9;
    for i in 0..n {
        for j in (i + 1)..n {
            if m.get(i, j).is_some_and(|v| v <= limit) {
                let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                if a != b {
                    parent[a.max(b)] = a.min(b);
                }
            }
        }
    }
    let mut by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = root(&mut parent, i);
        by_root.entry(r).or_default().push(i);
    }
    // Roots are the smallest index of their component, so BTreeMap order is
    // "by first row".
    by_root.into_values().collect()
}

/// One group of two or more genomes.
#[derive(Debug, Clone, PartialEq)]
pub struct SnpGroup {
    /// 1-based, in tree order: G1 is the topmost group of the full tree.
    pub number: usize,
    pub cluster: u32,
    /// Genome ids, in tree order.
    pub members: Vec<u32>,
    /// The widest pair inside the group, in SNPs (single linkage can chain
    /// past the threshold).
    pub diameter: f64,
}

impl SnpGroup {
    pub fn label(&self) -> String {
        format!("G{}", self.number)
    }
}

/// Every group at one threshold, plus what was covered.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GroupSet {
    pub threshold: u32,
    pub groups: Vec<SnpGroup>,
    /// Exact clusters whose matrix was available, with their genomes.
    pub covered: BTreeMap<u32, Vec<u32>>,
    /// Exact clusters of two or more genomes with no matrix to group by.
    pub missing: Vec<u32>,
}

impl GroupSet {
    /// Genomes that share a group with at least one other genome.
    pub fn grouped_genomes(&self) -> usize {
        self.groups.iter().map(|g| g.members.len()).sum()
    }

    /// Genome id to group index (into `groups`).
    pub fn index(&self) -> HashMap<u32, usize> {
        let mut out = HashMap::new();
        for (i, g) in self.groups.iter().enumerate() {
            for &m in &g.members {
                out.insert(m, i);
            }
        }
        out
    }
}

/// Group every matrix at `threshold`. `rank` orders genomes as the fully
/// opened tree shows them, top to bottom, so the numbering reads down the
/// tree and does not change when clusters are opened or closed.
pub fn compute(
    matrices: &[&SnpMatrix],
    threshold: u32,
    rank: &HashMap<u32, usize>,
    missing: Vec<u32>,
) -> GroupSet {
    let key = |g: u32| rank.get(&g).copied().unwrap_or(usize::MAX);
    let mut groups: Vec<SnpGroup> = Vec::new();
    let mut covered = BTreeMap::new();
    for m in matrices {
        covered.insert(m.cluster, m.ids.clone());
        for comp in components(m, f64::from(threshold)) {
            if comp.len() < 2 {
                continue;
            }
            let mut diameter = 0.0f64;
            for (a, &i) in comp.iter().enumerate() {
                for &j in &comp[a + 1..] {
                    if let Some(v) = m.get(i, j) {
                        diameter = diameter.max(v);
                    }
                }
            }
            let mut members: Vec<u32> = comp.iter().map(|&i| m.ids[i]).collect();
            members.sort_by_key(|&g| (key(g), g));
            groups.push(SnpGroup {
                number: 0,
                cluster: m.cluster,
                members,
                diameter,
            });
        }
    }
    groups.sort_by_key(|g| (key(g.members[0]), g.cluster, g.members[0]));
    for (i, g) in groups.iter_mut().enumerate() {
        g.number = i + 1;
    }
    GroupSet {
        threshold,
        groups,
        covered,
        missing,
    }
}

/// The groups as CSV: one line per genome of every grouped cluster, the group
/// left empty for a genome that is linked to nobody at this threshold.
pub fn csv(set: &GroupSet, name_of: impl Fn(u32) -> String) -> String {
    let mut out =
        String::from("genome,cluster,snp_group,group_size,threshold_snps,group_diameter_snps\n");
    let index = set.index();
    let mut seen: HashSet<u32> = HashSet::new();
    let quote = |s: String| {
        if s.contains([',', '"', '\n']) {
            format!("\"{}\"", s.replace('"', "\"\""))
        } else {
            s
        }
    };
    for g in &set.groups {
        for &m in &g.members {
            seen.insert(m);
            out.push_str(&format!(
                "{},c{},{},{},{},{}\n",
                quote(name_of(m)),
                g.cluster,
                g.label(),
                g.members.len(),
                set.threshold,
                crate::fmt::snps(g.diameter)
            ));
        }
    }
    for (cluster, ids) in &set.covered {
        for &m in ids {
            if !seen.contains(&m) && !index.contains_key(&m) {
                out.push_str(&format!(
                    "{},c{cluster},,1,{},0\n",
                    quote(name_of(m)),
                    set.threshold
                ));
            }
        }
    }
    out
}

/// The genomes of the fully opened tree, top to bottom: the order groups are
/// numbered in.
pub fn tree_rank(model: &UiModel) -> HashMap<u32, usize> {
    let all: HashSet<u32> = model.clusters.iter().map(|c| c.id).collect();
    let mut out = HashMap::new();
    let Some(root) = model.display_tree(&all) else {
        return out;
    };
    let mut stack: Vec<&DisplayNode> = vec![&root];
    while let Some(n) = stack.pop() {
        if let NodeKind::Genome { id, .. } = n.kind {
            let k = out.len();
            out.entry(id).or_insert(k);
        }
        // Reversed, so the first child is visited first.
        stack.extend(n.children.iter().rev());
    }
    out
}

/// The SNP groups the window shows, kept current as the model, the threshold
/// and the run folder change. Matrices are read once per cluster membership.
#[derive(Debug, Default)]
pub struct GroupCache {
    pub set: GroupSet,
    key: Option<(u64, u64, u32, Option<PathBuf>)>,
    /// Cluster id to (sorted member names, matrix).
    matrices: HashMap<u32, (Vec<String>, std::sync::Arc<SnpMatrix>)>,
    /// Clusters whose matrix file was tried and not usable: (membership,
    /// whether the run had finished). Tried again only when either changes.
    tried_tsv: HashMap<u32, (Vec<String>, bool)>,
}

impl GroupCache {
    /// Recompute when anything the groups depend on changed. Returns true
    /// when it did.
    pub fn refresh(&mut self, model: &UiModel, run_dir: Option<&Path>, threshold: u32) -> bool {
        let key = (
            model.version,
            model.detail_version,
            threshold,
            run_dir.map(Path::to_path_buf),
        );
        if self.key.as_ref() == Some(&key) {
            return false;
        }
        let ids = model.tree_name_index();
        let mut mats: Vec<std::sync::Arc<SnpMatrix>> = Vec::new();
        let mut missing = Vec::new();
        for c in &model.clusters {
            if c.size() < 2 || c.focus < crate::fmt::FOCUS_EXACT || c.approximate {
                continue;
            }
            let mut sig = c.member_names.clone();
            sig.sort_unstable();
            let members: HashSet<u32> = c.members.iter().copied().collect();
            let fits = |m: &SnpMatrix| {
                m.ids.len() == members.len() && m.ids.iter().all(|g| members.contains(g))
            };
            if let Some((s, m)) = self.matrices.get(&c.id) {
                if *s == sig {
                    mats.push(std::sync::Arc::clone(m));
                    continue;
                }
            }
            let mut found = model
                .details
                .get(&c.id)
                .and_then(|d| from_detail(d, &ids))
                .filter(|m| fits(m));
            if found.is_none() {
                if let Some(dir) = run_dir {
                    let attempt = (sig.clone(), model.finished);
                    if self.tried_tsv.get(&c.id) != Some(&attempt) {
                        found = read_tsv(dir, c.id, &ids).ok().filter(|m| fits(m));
                        if found.is_none() {
                            self.tried_tsv.insert(c.id, attempt);
                        }
                    }
                }
            }
            match found {
                Some(m) => {
                    let m = std::sync::Arc::new(m);
                    self.matrices.insert(c.id, (sig, std::sync::Arc::clone(&m)));
                    mats.push(m);
                }
                None => missing.push(c.id),
            }
        }
        let rank = tree_rank(model);
        let refs: Vec<&SnpMatrix> = mats.iter().map(|m| m.as_ref()).collect();
        self.set = compute(&refs, threshold, &rank, missing);
        self.key = Some(key);
        true
    }
}

/// The slider's scale: fine steps at the low end, where outbreak thresholds
/// live (0–20 one by one), coarser above. Maps a slider position in `[0, 1]`
/// to SNPs and back.
pub fn slider_to_snps(t: f64) -> u32 {
    let t = if t.is_finite() {
        t.clamp(0.0, 1.0)
    } else {
        0.0
    };
    // Two linear pieces: the first 60 % of the track is 0..=20.
    let v = if t <= 0.6 {
        t / 0.6 * 20.0
    } else {
        20.0 + (t - 0.6) / 0.4 * 80.0
    };
    v.round().clamp(0.0, 100.0) as u32
}

pub fn snps_to_slider(n: u32) -> f64 {
    let n = f64::from(n.min(100));
    if n <= 20.0 {
        n / 20.0 * 0.6
    } else {
        0.6 + (n - 20.0) / 80.0 * 0.4
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Five genomes: a–b 2 SNPs, b–c 4, c–d 30, d–e 1; ids 10..14.
    fn chain() -> SnpMatrix {
        let d = [
            [0.0, 2.0, 6.0, 36.0, 37.0],
            [2.0, 0.0, 4.0, 34.0, 35.0],
            [6.0, 4.0, 0.0, 30.0, 31.0],
            [36.0, 34.0, 30.0, 0.0, 1.0],
            [37.0, 35.0, 31.0, 1.0, 0.0],
        ];
        let rows: Vec<Vec<Option<f64>>> = d
            .iter()
            .map(|r| r.iter().map(|&v| Some(v)).collect())
            .collect();
        SnpMatrix::new(7, vec![10, 11, 12, 13, 14], &rows).unwrap()
    }

    #[test]
    fn single_linkage_chains_through_close_steps() {
        let m = chain();
        assert_eq!(components(&m, 5.0), vec![vec![0, 1, 2], vec![3, 4]]);
        // At 3, b–c (4 SNPs) no longer links: c is alone.
        assert_eq!(components(&m, 3.0), vec![vec![0, 1], vec![2], vec![3, 4]]);
        // At 0 nothing links; at 30 everything does.
        assert_eq!(components(&m, 0.0).len(), 5);
        assert_eq!(components(&m, 30.0), vec![vec![0, 1, 2, 3, 4]]);
    }

    #[test]
    fn groups_are_numbered_down_the_tree_and_singletons_are_left_out() {
        let m = chain();
        // The tree shows the d–e pair above a–b–c.
        let rank: HashMap<u32, usize> = [(13, 0), (14, 1), (10, 2), (11, 3), (12, 4)]
            .into_iter()
            .collect();
        let s = compute(&[&m], 5, &rank, vec![]);
        assert_eq!(s.groups.len(), 2);
        assert_eq!(s.groups[0].label(), "G1");
        assert_eq!(s.groups[0].members, vec![13, 14]);
        assert_eq!(s.groups[1].members, vec![10, 11, 12]);
        // a–c are 6 apart although every link is at most 5.
        assert_eq!(s.groups[1].diameter, 6.0);
        assert_eq!(s.grouped_genomes(), 5);
        let s3 = compute(&[&m], 3, &rank, vec![]);
        assert_eq!(s3.groups.len(), 2, "c drops out, singletons are not groups");
        assert_eq!(s3.grouped_genomes(), 4);
        assert_eq!(s3.index().get(&12), None);
    }

    #[test]
    fn a_pair_that_shares_no_bases_never_links() {
        let rows = vec![
            vec![Some(0.0), None, Some(9.0)],
            vec![None, Some(0.0), Some(9.0)],
            vec![Some(9.0), Some(9.0), Some(0.0)],
        ];
        let m = SnpMatrix::new(1, vec![1, 2, 3], &rows).unwrap();
        // 1 and 2 share no bases, but both are 9 from 3, which joins them.
        assert_eq!(components(&m, 100.0), vec![vec![0, 1, 2]]);
        assert_eq!(components(&m, 5.0).len(), 3);
        assert!(SnpMatrix::new(1, vec![1, 2], &rows).is_err(), "not square");
    }

    #[test]
    fn the_run_folders_tsv_parses_with_na_and_names() {
        let tsv = "genome\ta\tb\tc\na\t0\t3\tNA\nb\t3\t0\t12\nc\tNA\t12\t0\n";
        let ids: HashMap<String, u32> = [("a", 5u32), ("b", 6), ("c", 7)]
            .iter()
            .map(|(n, i)| (n.to_string(), *i))
            .collect();
        let m = parse_tsv(tsv, 3, &ids).unwrap();
        assert_eq!(m.ids, vec![5, 6, 7]);
        assert_eq!(m.get(0, 1), Some(3.0));
        assert_eq!(m.get(0, 2), None);
        assert_eq!(components(&m, 5.0), vec![vec![0, 1], vec![2]]);
        assert_eq!(components(&m, 12.0), vec![vec![0, 1, 2]]);
        // An unknown name or a shuffled row is refused, not guessed at.
        assert!(parse_tsv("genome\tz\nz\t0\n", 3, &ids).is_err());
        assert!(parse_tsv("genome\ta\tb\nb\t0\t1\na\t1\t0\n", 3, &ids).is_err());
        assert!(parse_tsv("", 3, &ids).is_err());
    }

    #[test]
    fn a_detail_event_gives_a_matrix_only_when_exact_and_present() {
        let ids: HashMap<String, u32> = [("a", 1u32), ("b", 2)]
            .iter()
            .map(|(n, i)| (n.to_string(), *i))
            .collect();
        let mut d = ClusterDetailInfo {
            cluster: 4,
            focus: 2,
            ids: vec![0, 1],
            names: vec!["a".into(), "b".into()],
            matrix: Some(vec![vec![0.0, 3.0], vec![3.0, 0.0]]),
            raw_matrix: None,
            core_bp: 0,
            pass_snps: vec![],
            masked_by_rule: vec![],
            snp_thresholds: vec![],
            snp_groups: vec![],
        };
        let m = from_detail(&d, &ids).unwrap();
        assert_eq!(
            (m.cluster, m.ids.clone(), m.get(0, 1)),
            (4, vec![1, 2], Some(3.0))
        );
        d.focus = 1;
        assert!(
            from_detail(&d, &ids).is_none(),
            "Focus 1 values are estimates"
        );
        d.focus = 2;
        d.matrix = None;
        assert!(
            from_detail(&d, &ids).is_none(),
            "large clusters carry no matrix"
        );
    }

    #[test]
    fn the_csv_lists_every_genome_with_its_group() {
        let m = chain();
        let rank: HashMap<u32, usize> = HashMap::new();
        let s = compute(&[&m], 3, &rank, vec![]);
        let text = csv(&s, |g| format!("g{g}"));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "genome,cluster,snp_group,group_size,threshold_snps,group_diameter_snps"
        );
        assert_eq!(lines.len(), 6, "{text}");
        assert!(lines.contains(&"g10,c7,G1,2,3,2"));
        assert!(lines.contains(&"g12,c7,,1,3,0"), "{text}");
        assert!(text.contains("g13,c7,G2,2,3,1"));
    }

    /// A model with one exact cluster of four genomes (a 2-SNP pair and a
    /// 3-SNP pair, 20 apart) and its `ClusterDetail` matrix.
    fn model_with_detail() -> UiModel {
        use infinity_core::events::{ClusterFocusInfo, Event};
        let mut m = UiModel::default();
        for i in 0..4u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(0.0, crate::state::clustered_event(&[(0, 0, &[0, 1, 2, 3])]));
        m.apply(
            0.0,
            Event::BackboneTree {
                newick: "g0;".into(),
                pivots: vec!["g0".into()],
            },
        );
        m.apply(
            0.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "((g2:1,g3:2):10,(g0:1,g1:1):10);".into(),
                members: vec![],
            }),
        );
        let d = [
            [0.0, 2.0, 21.0, 22.0],
            [2.0, 0.0, 21.0, 22.0],
            [21.0, 21.0, 0.0, 3.0],
            [22.0, 22.0, 3.0, 0.0],
        ];
        m.apply(
            0.0,
            Event::ClusterDetail(Box::new(ClusterDetailInfo {
                cluster: 0,
                focus: 2,
                ids: vec![0, 1, 2, 3],
                names: (0..4).map(|i| format!("g{i}")).collect(),
                matrix: Some(d.iter().map(|r| r.to_vec()).collect()),
                raw_matrix: None,
                core_bp: 0,
                pass_snps: vec![],
                masked_by_rule: vec![],
                snp_thresholds: vec![],
                snp_groups: vec![],
            })),
        );
        m
    }

    #[test]
    fn the_cache_groups_from_the_detail_event_and_follows_the_slider() {
        let m = model_with_detail();
        let mut c = GroupCache::default();
        assert!(c.refresh(&m, None, 5));
        assert!(
            !c.refresh(&m, None, 5),
            "nothing changed, nothing recomputed"
        );
        assert_eq!(c.set.groups.len(), 2);
        // The tree shows g2, g3 first: that pair is G1.
        assert_eq!(c.set.groups[0].members, vec![2, 3]);
        assert_eq!(c.set.groups[1].members, vec![0, 1]);
        assert!(c.set.missing.is_empty());
        assert!(c.refresh(&m, None, 2));
        assert_eq!(c.set.groups.len(), 1, "only the 2-SNP pair at 2");
        assert!(c.refresh(&m, None, 25));
        assert_eq!(c.set.groups.len(), 1);
        assert_eq!(c.set.groups[0].members.len(), 4);
        assert_eq!(c.set.groups[0].diameter, 22.0);
    }

    #[test]
    fn without_a_detail_the_run_folders_matrix_is_read() {
        let mut m = model_with_detail();
        m.details.clear();
        m.detail_version += 1;
        let dir = std::env::temp_dir().join(format!("infinity-groups-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut c = GroupCache::default();
        c.refresh(&m, Some(&dir), 5);
        assert_eq!(c.set.missing, vec![0], "no matrix anywhere yet");
        std::fs::write(
            tsv_path(&dir, 0),
            "genome\tg0\tg1\tg2\tg3\ng0\t0\t2\t21\t22\ng1\t2\t0\t21\t22\n\
             g2\t21\t21\t0\t3\ng3\t22\t22\t3\t0\n",
        )
        .unwrap();
        // The run finishing is what makes the file worth trying again.
        m.finished = true;
        m.version += 1;
        c.refresh(&m, Some(&dir), 5);
        assert!(c.set.missing.is_empty());
        assert_eq!(c.set.groups.len(), 2);
        // A file for other genomes (an older run in the same folder) is refused.
        std::fs::write(tsv_path(&dir, 0), "genome\tg0\tg1\ng0\t0\t1\ng1\t1\t0\n").unwrap();
        let mut fresh = GroupCache::default();
        fresh.refresh(&m, Some(&dir), 5);
        assert_eq!(fresh.set.missing, vec![0]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn the_slider_is_fine_at_the_low_end_and_round_trips() {
        assert_eq!(slider_to_snps(0.0), 0);
        assert_eq!(slider_to_snps(1.0), 100);
        assert_eq!(slider_to_snps(0.6), 20);
        // One SNP per 3 % of the track at the low end.
        assert_eq!(slider_to_snps(0.03), 1);
        for n in 0..=100u32 {
            assert_eq!(slider_to_snps(snps_to_slider(n)), n, "{n}");
        }
        assert_eq!(slider_to_snps(f64::NAN), 0);
    }
}
