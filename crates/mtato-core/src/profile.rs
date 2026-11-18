//! Profiles, pairwise distances and cluster trees (6.9).
//!
//! A member's *profile* is everything that distinguishes it from its pivot: its
//! variants plus the pivot intervals it does not cover. Distances between two
//! members are then derived from their profiles alone, without looking at either
//! sequence again — which is what makes the cost of a cluster quadratic in the
//! number of *differences* rather than in genome size.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::compare::{DiffPiece, MemberComparison, PieceKey};
use crate::interval::{merge, total_len, union_len, Iv};
use crate::pivot::PivotIndex;
use crate::tree::{midpoint_root, neighbor_joining, Tree};
use crate::variant::{Variant, VariantKind};

/// The pivot sequence with contigs concatenated, so that a profile can be
/// indexed by a single flat position. Built once per cluster.
pub fn concat_pivot(pivot: &PivotIndex) -> Vec<u8> {
    let mut out = Vec::with_capacity(pivot.total_len as usize);
    for c in &pivot.contigs {
        out.extend_from_slice(&c.seq);
    }
    out
}

/// Start offset of each pivot contig in flat coordinates.
pub fn contig_offsets(pivot: &PivotIndex) -> Vec<u64> {
    let mut out = Vec::with_capacity(pivot.contig_lens.len());
    let mut acc = 0u64;
    for &l in &pivot.contig_lens {
        out.push(acc);
        acc += l;
    }
    out
}

/// A key that identifies an indel for set comparison between two profiles.
///
/// Two members share an indel when the keys are equal, which is why variants
/// are left-aligned in 6.7: without canonical coordinates the same event seen
/// from two members would look like two events.
pub type IndelKey = (u64, u64, VariantKind, String, String);

/// Everything that distinguishes a member from its pivot (section 5).
#[derive(Debug, Clone, PartialEq)]
pub struct Profile {
    pub genome_id: u32,
    /// Sorted variants in pivot coordinates. Empty for the pivot itself.
    pub variants: Vec<Variant>,
    /// SNP positions, in flat pivot coordinates, mapped to the member's base.
    pub snp_map: HashMap<u64, u8>,
    /// Sorted indel keys, in flat pivot coordinates.
    pub indels: Vec<IndelKey>,
    /// Pivot intervals this member does not cover, merged, flat coordinates.
    pub uncovered: Vec<Iv>,
    /// Focus 1 records: the difference pieces, before they were resolved.
    pub pieces: Vec<PieceKey>,
    /// Best focus this profile represents: 1 or 2.
    pub focus: u8,
}

impl Profile {
    /// The pivot's own profile, which is empty (6.9).
    pub fn pivot(genome_id: u32) -> Self {
        Self {
            genome_id,
            variants: Vec::new(),
            snp_map: HashMap::new(),
            indels: Vec::new(),
            uncovered: Vec::new(),
            pieces: Vec::new(),
            focus: 2,
        }
    }

    /// Build a Focus 2 profile from a finished comparison.
    pub fn from_comparison(cmp: &MemberComparison, offsets: &[u64]) -> Self {
        let mut snp_map = HashMap::new();
        let mut indels = Vec::new();
        for v in &cmp.variants {
            let base = offsets.get(v.pivot_contig as usize).copied().unwrap_or(0);
            match v.kind {
                VariantKind::Snp => {
                    if let Some(b) = v.alt_allele.as_bytes().first() {
                        snp_map.insert(base + v.pivot_start, *b);
                    }
                }
                VariantKind::Insertion | VariantKind::Deletion => {
                    indels.push((
                        base + v.pivot_start,
                        base + v.pivot_end,
                        v.kind,
                        v.ref_allele.clone(),
                        v.alt_allele.clone(),
                    ));
                }
                VariantKind::Complex => {
                    // A COMPLEX interval is not covered, so it is excluded from
                    // the comparable denominator rather than counted as a
                    // difference of a specific kind.
                }
            }
        }
        indels.sort();
        indels.dedup();
        Self {
            genome_id: cmp.genome_id,
            variants: cmp.variants.clone(),
            snp_map,
            indels,
            uncovered: merge(cmp.uncovered.clone()),
            pieces: cmp.focus1.pieces.iter().map(|p| p.cache_key()).collect(),
            focus: 2,
        }
    }

    /// Build a Focus 1 profile: difference pieces only, no variants yet.
    pub fn from_focus1(genome_id: u32, pieces: &[DiffPiece], uncovered: Vec<Iv>) -> Self {
        let mut keys: Vec<PieceKey> = pieces.iter().map(|p| p.cache_key()).collect();
        keys.sort();
        keys.dedup();
        Self {
            genome_id,
            variants: Vec::new(),
            snp_map: HashMap::new(),
            indels: Vec::new(),
            uncovered: merge(uncovered),
            pieces: keys,
            focus: 1,
        }
    }

    /// A hash of everything that defines this profile, used for the grouping of
    /// 6.9 so that identical genomes are handled once.
    pub fn identity(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        for v in &self.variants {
            v.pivot_contig.hash(&mut h);
            v.pivot_start.hash(&mut h);
            v.pivot_end.hash(&mut h);
            v.ref_allele.hash(&mut h);
            v.alt_allele.hash(&mut h);
            (v.kind as u8).hash(&mut h);
        }
        for iv in &self.uncovered {
            iv.hash(&mut h);
        }
        h.finish()
    }
}

/// The distance between two members of a cluster (6.9).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PairDistance {
    pub snps: u64,
    pub indels: u64,
    /// Pivot positions where the two can actually be compared.
    pub comparable: u64,
    /// `snps / comparable`, substitutions per site. `None` when nothing is
    /// comparable, so a zero is never invented (section 0 rule 8).
    pub dist: Option<f64>,
}

/// Pairwise distance inside a cluster.
///
/// ```text
/// comparable(A,B) = pivot length - | uncovered(A) UNION uncovered(B) |
/// snps(A,B)       = comparable pivot positions where A's base differs from B's
/// indels(A,B)     = insertions/deletions inside comparable intervals present
///                   in one profile and not identically in the other
/// dist(A,B)       = snps(A,B) / comparable(A,B)
/// ```
///
/// A member without a SNP at a position has the pivot's base there.
pub fn pair_distance(
    a: &Profile,
    b: &Profile,
    pivot_len: u64,
    pivot_concat: &[u8],
) -> PairDistance {
    let uncovered_union = union_len(&a.uncovered, &b.uncovered);
    let comparable = pivot_len.saturating_sub(uncovered_union);

    // Union of the two SNP position sets; everything else is the pivot base on
    // both sides and therefore cannot differ.
    let mut positions: Vec<u64> = Vec::with_capacity(a.snp_map.len() + b.snp_map.len());
    positions.extend(a.snp_map.keys().copied());
    positions.extend(b.snp_map.keys().copied());
    positions.sort_unstable();
    positions.dedup();

    let mut snps = 0u64;
    for p in positions {
        if p >= pivot_len {
            continue;
        }
        if uncovered_union > 0 && in_any(&[&a.uncovered, &b.uncovered], p) {
            continue;
        }
        let base = pivot_concat.get(p as usize).copied().unwrap_or(b'N');
        let ba = a.snp_map.get(&p).copied().unwrap_or(base);
        let bb = b.snp_map.get(&p).copied().unwrap_or(base);
        if ba != bb {
            snps += 1;
        }
    }

    // Indels: symmetric difference of the key sets, restricted to comparable
    // positions.
    let mut indels = 0u64;
    let (mut i, mut j) = (0, 0);
    while i < a.indels.len() && j < b.indels.len() {
        if a.indels[i] == b.indels[j] {
            i += 1;
            j += 1;
        } else if a.indels[i] < b.indels[j] {
            if comparable_position(a.indels[i].0, pivot_len, a, b) {
                indels += 1;
            }
            i += 1;
        } else {
            if comparable_position(b.indels[j].0, pivot_len, a, b) {
                indels += 1;
            }
            j += 1;
        }
    }
    for k in &a.indels[i..] {
        if comparable_position(k.0, pivot_len, a, b) {
            indels += 1;
        }
    }
    for k in &b.indels[j..] {
        if comparable_position(k.0, pivot_len, a, b) {
            indels += 1;
        }
    }

    PairDistance {
        snps,
        indels,
        comparable,
        dist: if comparable == 0 {
            None
        } else {
            Some(snps as f64 / comparable as f64)
        },
    }
}

fn comparable_position(pos: u64, pivot_len: u64, a: &Profile, b: &Profile) -> bool {
    pos < pivot_len && !in_any(&[&a.uncovered, &b.uncovered], pos)
}

fn in_any(lists: &[&[Iv]], pos: u64) -> bool {
    lists.iter().any(|l| crate::interval::contains(l, pos))
}

/// The Focus 1 estimate of 6.9: the number of difference-piece records that are
/// in one list and not in the other.
///
/// This is an *approximation* and must always be shown with the "≈" marker
/// (section 0 rule 8). The true SNP count is at least half of it, rounded up,
/// because one differing piece can hold anything from one SNP to a large
/// rearrangement.
pub fn focus1_estimate(a: &Profile, b: &Profile) -> u64 {
    let (mut i, mut j) = (0, 0);
    let mut diff = 0u64;
    while i < a.pieces.len() && j < b.pieces.len() {
        if a.pieces[i] == b.pieces[j] {
            i += 1;
            j += 1;
        } else if a.pieces[i] < b.pieces[j] {
            diff += 1;
            i += 1;
        } else {
            diff += 1;
            j += 1;
        }
    }
    diff += (a.pieces.len() - i) as u64 + (b.pieces.len() - j) as u64;
    diff
}

/// The lower bound on the true SNP count implied by a Focus 1 estimate.
pub fn focus1_snp_lower_bound(estimate: u64) -> u64 {
    estimate.div_ceil(2)
}

/// Which quantity a tree is built from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TreeBasis {
    /// Sketch distance: the Focus 0 star.
    Sketch,
    /// Number of differing difference-piece records: Focus 1, approximate.
    Focus1Estimate,
    /// SNP counts: Focus 2, exact.
    #[default]
    Snps,
    /// Substitutions per site: what the full tree uses inside clusters so that
    /// units match the backbone (6.9).
    SubsPerSite,
}

impl TreeBasis {
    /// True when the distances behind this basis are estimates rather than
    /// counts, which is what the GUI shows with an "≈" marker.
    pub fn is_approximate(&self) -> bool {
        matches!(self, TreeBasis::Sketch | TreeBasis::Focus1Estimate)
    }

    pub fn focus(&self) -> u8 {
        match self {
            TreeBasis::Sketch => 0,
            TreeBasis::Focus1Estimate => 1,
            TreeBasis::Snps | TreeBasis::SubsPerSite => 2,
        }
    }
}

/// A full pairwise matrix over a set of profiles.
#[derive(Debug, Clone, Default)]
pub struct DistanceMatrix {
    pub ids: Vec<u32>,
    pub snps: Vec<Vec<u64>>,
    pub indels: Vec<Vec<u64>>,
    pub comparable: Vec<Vec<u64>>,
    pub dist: Vec<Vec<f64>>,
    /// Focus 1 estimates, filled only when the profiles are Focus 1.
    pub estimates: Vec<Vec<u64>>,
    pub basis: TreeBasis,
}

impl DistanceMatrix {
    /// The matrix actually used to build a tree, given the basis.
    fn matrix_for(&self, basis: TreeBasis) -> Vec<Vec<f64>> {
        let n = self.ids.len();
        match basis {
            TreeBasis::Snps => (0..n)
                .map(|i| (0..n).map(|j| self.snps[i][j] as f64).collect())
                .collect(),
            TreeBasis::Focus1Estimate => (0..n)
                .map(|i| (0..n).map(|j| self.estimates[i][j] as f64).collect())
                .collect(),
            TreeBasis::SubsPerSite | TreeBasis::Sketch => self.dist.clone(),
        }
    }
}

/// Build the pairwise matrix for the profiles of one cluster.
pub fn distance_matrix(
    profiles: &[&Profile],
    pivot_len: u64,
    pivot_concat: &[u8],
) -> DistanceMatrix {
    let n = profiles.len();
    let mut snps = vec![vec![0u64; n]; n];
    let mut indels = vec![vec![0u64; n]; n];
    let mut comparable = vec![vec![0u64; n]; n];
    let mut dist = vec![vec![0f64; n]; n];
    let mut estimates = vec![vec![0u64; n]; n];

    for i in 0..n {
        comparable[i][i] = pivot_len.saturating_sub(total_len(&profiles[i].uncovered));
        for j in (i + 1)..n {
            let d = pair_distance(profiles[i], profiles[j], pivot_len, pivot_concat);
            snps[i][j] = d.snps;
            snps[j][i] = d.snps;
            indels[i][j] = d.indels;
            indels[j][i] = d.indels;
            comparable[i][j] = d.comparable;
            comparable[j][i] = d.comparable;
            // `dist` is `None` when nothing is comparable, and neighbour
            // joining still has to be given a number. Zero would assert "these
            // two are identical", inventing a measurement out of the absence of
            // one (section 0 rule 8) and pulling unrelated genomes together in
            // the tree. 1.0 is the largest distance the sketch metric of 6.4
            // can produce, so an incomparable pair is placed as far apart as
            // the scale allows rather than as close as it can be.
            // `cluster_<id>.dist.tsv` is unaffected: it tests `comparable`
            // itself and writes `NA`.
            let dv = d.dist.unwrap_or(1.0);
            dist[i][j] = dv;
            dist[j][i] = dv;
            let e = focus1_estimate(profiles[i], profiles[j]);
            estimates[i][j] = e;
            estimates[j][i] = e;
        }
    }
    DistanceMatrix {
        ids: profiles.iter().map(|p| p.genome_id).collect(),
        snps,
        indels,
        comparable,
        dist,
        estimates,
        basis: TreeBasis::Snps,
    }
}

/// The star tree of Focus 0: until a cluster is refined, its members hang from
/// the pivot, each at its sketch distance (6.4).
pub fn star_tree(pivot_name: &str, members: &[(String, f64)]) -> Tree {
    let mut children = vec![Tree::Leaf {
        name: pivot_name.to_string(),
        len: 0.0,
    }];
    for (name, d) in members {
        children.push(Tree::Leaf {
            name: name.clone(),
            len: (*d).max(0.0),
        });
    }
    if children.len() == 1 {
        return children.pop().unwrap();
    }
    Tree::Internal { len: 0.0, children }
}

/// Neighbour-joining tree of one cluster, midpoint rooted, with negative branch
/// lengths set to zero (6.9).
pub fn cluster_tree(names: &[String], matrix: &DistanceMatrix, basis: TreeBasis) -> Tree {
    let n = names.len();
    if n == 0 {
        return Tree::Internal {
            len: 0.0,
            children: Vec::new(),
        };
    }
    if n == 1 {
        return Tree::Leaf {
            name: names[0].clone(),
            len: 0.0,
        };
    }
    let d = matrix.matrix_for(basis);
    let mut t = neighbor_joining(names, &d);
    t.clamp_negative_lengths();
    midpoint_root(t)
}

/// Group profiles that are identical, so distances are computed between groups
/// only and the members of a group attach to the tree as zero-length siblings
/// (6.9, 6.10).
///
/// Returns `(group representative indices, members per group)`. Groups are
/// ordered by the genome id of their first member, so the result does not
/// depend on hash iteration order.
pub fn group_identical(profiles: &[&Profile]) -> (Vec<usize>, Vec<Vec<usize>>) {
    let mut by_hash: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, p) in profiles.iter().enumerate() {
        by_hash.entry(p.identity()).or_default().push(i);
    }
    // Hash equality is not identity equality: verify before merging, so a hash
    // collision can never make two different genomes look the same.
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut keys: Vec<u64> = by_hash.keys().copied().collect();
    keys.sort_unstable();
    for h in keys {
        let idxs = &by_hash[&h];
        for &i in idxs {
            let target = groups.iter_mut().find(|g| {
                let rep = profiles[g[0]];
                profiles[i].variants == rep.variants && profiles[i].uncovered == rep.uncovered
            });
            match target {
                Some(g) => g.push(i),
                None => groups.push(vec![i]),
            }
        }
    }
    for g in &mut groups {
        g.sort_unstable();
    }
    groups.sort_by_key(|g| profiles[g[0]].genome_id);
    let reps = groups.iter().map(|g| g[0]).collect();
    (reps, groups)
}

/// Expand a tree whose leaves are group representatives into one whose leaves
/// are genomes, attaching the other members of each group as zero-length
/// siblings (6.9).
pub fn expand_groups(tree: Tree, groups: &[Vec<usize>], names: &[String]) -> Tree {
    let mut t = tree;
    for g in groups {
        if g.len() <= 1 {
            continue;
        }
        let rep = &names[g[0]];
        let siblings: Vec<Tree> = g
            .iter()
            .map(|&i| Tree::Leaf {
                name: names[i].clone(),
                len: 0.0,
            })
            .collect();
        // The representative keeps the branch length that leads to the group.
        let subtree = Tree::Internal {
            len: 0.0,
            children: siblings,
        };
        crate::tree::graft(&mut t, rep, subtree);
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::Variant;

    fn var(kind: VariantKind, pos: u64, ref_a: &str, alt: &str) -> Variant {
        Variant {
            kind,
            pivot_contig: 0,
            pivot_start: pos,
            pivot_end: pos + ref_a.len() as u64,
            ref_allele: ref_a.into(),
            alt_allele: alt.into(),
            genome_contig: 0,
            genome_start: pos,
            genome_end: pos + ref_a.len() as u64,
            strand: '+',
            flags: 0,
        }
    }

    fn profile(id: u32, variants: Vec<Variant>, uncovered: Vec<Iv>) -> Profile {
        let cmp = MemberComparison {
            genome_id: id,
            focus1: crate::compare::Focus1 {
                pieces: Vec::new(),
                breakpoints: Vec::new(),
                matched: Vec::new(),
                unplaced: Vec::new(),
                covered: Vec::new(),
                not_called: Vec::new(),
                variants: Vec::new(),
                novel_bp: 0,
                counters: Default::default(),
            },
            variants,
            uncovered,
            duplicated: Vec::new(),
            novel_bp: 0,
            counters: Default::default(),
        };
        Profile::from_comparison(&cmp, &[0])
    }

    #[test]
    fn two_identical_profiles_are_at_distance_zero() {
        let pivot = b"ACGTACGTACGTACGTACGT".to_vec();
        let a = profile(0, vec![], vec![]);
        let b = profile(1, vec![], vec![]);
        let d = pair_distance(&a, &b, pivot.len() as u64, &pivot);
        assert_eq!(d.snps, 0);
        assert_eq!(d.indels, 0);
        assert_eq!(d.comparable, pivot.len() as u64);
        assert_eq!(d.dist, Some(0.0));
    }

    #[test]
    fn one_snp_in_one_member_gives_a_distance_of_one() {
        let pivot = b"ACGTACGTACGTACGTACGT".to_vec();
        let a = profile(0, vec![], vec![]);
        let b = profile(1, vec![var(VariantKind::Snp, 5, "A", "T")], vec![]);
        let d = pair_distance(&a, &b, pivot.len() as u64, &pivot);
        assert_eq!(d.snps, 1);
        assert_eq!(d.comparable, 20);
        assert!((d.dist.unwrap() - 1.0 / 20.0).abs() < 1e-12);
    }

    #[test]
    fn the_same_snp_in_both_members_is_not_a_difference() {
        let pivot = b"ACGTACGTACGTACGTACGT".to_vec();
        let a = profile(0, vec![var(VariantKind::Snp, 5, "A", "T")], vec![]);
        let b = profile(1, vec![var(VariantKind::Snp, 5, "A", "T")], vec![]);
        assert_eq!(pair_distance(&a, &b, 20, &pivot).snps, 0);
    }

    #[test]
    fn two_different_snps_at_one_position_count_once() {
        let pivot = b"ACGTACGTACGTACGTACGT".to_vec();
        let a = profile(0, vec![var(VariantKind::Snp, 5, "A", "T")], vec![]);
        let b = profile(1, vec![var(VariantKind::Snp, 5, "A", "G")], vec![]);
        assert_eq!(pair_distance(&a, &b, 20, &pivot).snps, 1);
    }

    #[test]
    fn uncovered_regions_shrink_the_comparable_denominator() {
        let pivot = vec![b'A'; 100];
        let a = profile(0, vec![], vec![(0, 10)]);
        let b = profile(1, vec![], vec![(5, 20)]);
        let d = pair_distance(&a, &b, 100, &pivot);
        // union of [0,10) and [5,20) is [0,20) = 20 bp
        assert_eq!(d.comparable, 80);
    }

    #[test]
    fn a_snp_inside_an_uncovered_region_is_not_counted() {
        let pivot = vec![b'A'; 100];
        let a = profile(0, vec![var(VariantKind::Snp, 5, "A", "T")], vec![]);
        let b = profile(1, vec![], vec![(0, 10)]);
        // Position 5 is uncovered in b, so it is not comparable.
        assert_eq!(pair_distance(&a, &b, 100, &pivot).snps, 0);
    }

    #[test]
    fn indels_present_in_one_profile_only_are_counted() {
        let pivot = b"ACGTACGTACGTACGTACGT".to_vec();
        let a = profile(0, vec![], vec![]);
        let b = profile(1, vec![var(VariantKind::Insertion, 5, "", "TTT")], vec![]);
        let d = pair_distance(&a, &b, 20, &pivot);
        assert_eq!(d.indels, 1);
        assert_eq!(d.snps, 0);
    }

    #[test]
    fn identical_indels_are_not_a_difference() {
        let pivot = b"ACGTACGTACGTACGTACGT".to_vec();
        let v = || var(VariantKind::Deletion, 5, "AC", "");
        let a = profile(0, vec![v()], vec![]);
        let b = profile(1, vec![v()], vec![]);
        assert_eq!(pair_distance(&a, &b, 20, &pivot).indels, 0);
    }

    #[test]
    fn nothing_comparable_gives_no_distance_rather_than_zero() {
        let pivot = vec![b'A'; 10];
        let a = profile(0, vec![], vec![(0, 10)]);
        let b = profile(1, vec![], vec![(0, 10)]);
        let d = pair_distance(&a, &b, 10, &pivot);
        assert_eq!(d.comparable, 0);
        assert_eq!(d.dist, None, "0/0 must not be reported as 0");
    }

    #[test]
    fn focus1_estimate_counts_the_symmetric_difference() {
        let key = |h: u64| (h, 0u32, 0u64, 100u64);
        let a = Profile {
            pieces: vec![key(1), key(2)],
            ..Profile::pivot(0)
        };
        let b = Profile {
            pieces: vec![key(2), key(3)],
            ..Profile::pivot(1)
        };
        assert_eq!(focus1_estimate(&a, &b), 2);
        assert_eq!(focus1_snp_lower_bound(2), 1);
        assert_eq!(focus1_snp_lower_bound(5), 3);
        assert_eq!(focus1_estimate(&a, &a), 0);
    }

    #[test]
    fn identical_profiles_group_together_and_others_do_not() {
        let v = var(VariantKind::Snp, 5, "A", "T");
        let a = profile(0, vec![v.clone()], vec![]);
        let b = profile(1, vec![v.clone()], vec![]);
        let c = profile(2, vec![var(VariantKind::Snp, 9, "C", "G")], vec![]);
        let refs: Vec<&Profile> = vec![&a, &b, &c];
        let (reps, groups) = group_identical(&refs);
        assert_eq!(groups.len(), 2);
        assert_eq!(reps.len(), 2);
        assert!(groups.iter().any(|g| g == &vec![0, 1]));
        assert!(groups.iter().any(|g| g == &vec![2]));
    }

    #[test]
    fn grouping_is_deterministic() {
        let v = var(VariantKind::Snp, 5, "A", "T");
        let a = profile(0, vec![v.clone()], vec![]);
        let b = profile(1, vec![v.clone()], vec![]);
        let c = profile(2, vec![], vec![]);
        let refs: Vec<&Profile> = vec![&a, &b, &c];
        let g1 = group_identical(&refs).1;
        let g2 = group_identical(&refs).1;
        assert_eq!(g1, g2);
    }

    #[test]
    fn star_tree_hangs_every_member_off_the_pivot() {
        let t = star_tree(
            "pivot",
            &[("m1".to_string(), 0.001), ("m2".to_string(), 0.002)],
        );
        assert_eq!(t.leaf_count(), 3);
        let s = crate::tree::to_newick(&t);
        assert!(s.contains("pivot"), "{s}");
    }

    #[test]
    fn star_tree_of_a_singleton_is_a_bare_leaf() {
        let t = star_tree("only", &[]);
        assert_eq!(t.leaf_count(), 1);
    }

    #[test]
    fn cluster_tree_has_one_leaf_per_profile_and_no_negative_lengths() {
        let pivot = vec![b'A'; 1000];
        let profiles: Vec<Profile> = (0..6u32)
            .map(|i| {
                let vs: Vec<Variant> = (0..i)
                    .map(|j| var(VariantKind::Snp, (j * 37 + 11) as u64, "A", "T"))
                    .collect();
                profile(i, vs, vec![])
            })
            .collect();
        let refs: Vec<&Profile> = profiles.iter().collect();
        let m = distance_matrix(&refs, 1000, &pivot);
        let names: Vec<String> = (0..6).map(|i| format!("g{i}")).collect();
        let t = cluster_tree(&names, &m, TreeBasis::Snps);
        assert_eq!(t.leaf_count(), 6);
        assert!(t.total_length() >= 0.0);
        // Distances must be symmetric and zero on the diagonal.
        for i in 0..6 {
            assert_eq!(m.snps[i][i], 0);
            for j in 0..6 {
                assert_eq!(m.snps[i][j], m.snps[j][i]);
            }
        }
        // A nested series of SNP counts must come out monotone along the chain.
        assert!(m.snps[0][5] >= m.snps[0][1]);
    }

    #[test]
    fn expanding_a_group_adds_zero_length_siblings() {
        let names: Vec<String> = vec!["rep".into(), "twin".into(), "other".into()];
        let t = crate::tree::parse_newick("(rep:1,other:1);").unwrap();
        let groups = vec![vec![0usize, 1usize], vec![2usize]];
        let out = expand_groups(t, &groups, &names);
        assert_eq!(out.leaf_count(), 3);
        let s = crate::tree::to_newick(&out);
        assert!(s.contains("twin"), "{s}");
    }

    #[test]
    fn tree_basis_reports_approximation_honestly() {
        assert!(TreeBasis::Sketch.is_approximate());
        assert!(TreeBasis::Focus1Estimate.is_approximate());
        assert!(!TreeBasis::Snps.is_approximate());
        assert!(!TreeBasis::SubsPerSite.is_approximate());
        assert_eq!(TreeBasis::Sketch.focus(), 0);
        assert_eq!(TreeBasis::Focus1Estimate.focus(), 1);
        assert_eq!(TreeBasis::Snps.focus(), 2);
    }
}
