//! Pairwise numbers for two selected genomes (section 8: "Selecting two leaves
//! shows their pairwise distance, SNP count, and focus").
//!
//! The event stream carries each member's variants against the *pivot* and
//! never a pairwise matrix, so this module recomputes the definition of 6.9 from
//! those variant lists. That is exact rather than an estimate, for two reasons:
//! a member with no SNP at a position has the pivot's base there, and a variant's
//! alt allele always differs from the pivot's ref. Two members therefore differ
//! at a position exactly when their alt alleles differ, which is what
//! [`pair_from_maps`] counts — the pivot's own sequence is not needed.
//!
//! When the run has finished, `RunResults` carries the engine's own matrix and
//! [`pair_from_matrix`] uses it instead; that is the number written to
//! `cluster_<id>.dist.tsv`.

use std::collections::HashMap;

use infinity_core::events::MapUpdateInfo;
use infinity_core::interval::{merge, union_len, Iv};
use infinity_core::profile::{DistanceMatrix, IndelKey};
use infinity_core::variant::{Variant, VariantKind};

use crate::fmt;

/// Where a pairwise number came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairSource {
    /// The engine's own distance matrix (a finished run, or an opened folder).
    Matrix,
    /// Recomputed here from the two members' variant lists.
    Variants,
    /// Not computable yet: at least one of the two is below Focus 2, so the
    /// variant lists do not exist.
    Unavailable,
}

/// The numbers shown when exactly two leaves are selected.
#[derive(Debug, Clone, PartialEq)]
pub struct PairView {
    pub a: u32,
    pub b: u32,
    pub snps: u64,
    pub indels: u64,
    /// Pivot positions the two can actually be compared at (6.9).
    pub comparable: u64,
    /// `snps / comparable`, substitutions per site. `None` when nothing is
    /// comparable, so a zero is never invented (section 0 rule 8).
    pub dist: Option<f64>,
    pub focus: u8,
    /// True when the numbers are estimates and must carry the "≈" marker.
    pub approximate: bool,
    pub source: PairSource,
}

impl PairView {
    /// A pair that cannot be measured yet says so rather than showing zeros.
    pub fn unavailable(a: u32, b: u32, focus: u8) -> Self {
        Self {
            a,
            b,
            snps: 0,
            indels: 0,
            comparable: 0,
            dist: None,
            focus,
            approximate: true,
            source: PairSource::Unavailable,
        }
    }

    pub fn available(&self) -> bool {
        self.source != PairSource::Unavailable
    }

    /// The distance as text, with the "≈" marker when it is an estimate.
    pub fn dist_text(&self) -> String {
        match self.dist {
            Some(d) => fmt::distance(d, self.approximate),
            None => "n/a".to_string(),
        }
    }

    pub fn snps_text(&self) -> String {
        if !self.available() {
            return "not measured yet".to_string();
        }
        format!(
            "{}{}",
            fmt::approx_mark(self.approximate),
            fmt::count(self.snps)
        )
    }
}

/// Start offset of each pivot contig in flat coordinates. Mirrors
/// `profile::contig_offsets`, but from the lengths the event carries.
pub fn contig_offsets(lens: &[u64]) -> Vec<u64> {
    let mut out = Vec::with_capacity(lens.len());
    let mut acc = 0u64;
    for &l in lens {
        out.push(acc);
        acc = acc.saturating_add(l);
    }
    out
}

/// Total pivot length, the denominator of 6.9.
pub fn pivot_len(m: &MapUpdateInfo) -> u64 {
    m.pivot_contig_lens.iter().copied().sum()
}

/// Flat pivot position of a variant, or `None` when its contig is unknown.
pub fn flat_pos(offsets: &[u64], contig: u32, pos: u64) -> Option<u64> {
    offsets.get(contig as usize).map(|&o| o.saturating_add(pos))
}

/// SNP alt bases keyed by flat pivot position, exactly as
/// `Profile::from_comparison` builds them.
pub fn snp_map(m: &MapUpdateInfo) -> HashMap<u64, u8> {
    let offsets = contig_offsets(&m.pivot_contig_lens);
    let mut out = HashMap::new();
    for v in &m.variants {
        if v.kind != VariantKind::Snp {
            continue;
        }
        if let (Some(p), Some(b)) = (
            flat_pos(&offsets, v.pivot_contig, v.pivot_start),
            v.alt_allele.as_bytes().first().copied(),
        ) {
            out.insert(p, b);
        }
    }
    out
}

/// Sorted, de-duplicated indel keys in flat pivot coordinates.
pub fn indel_keys(m: &MapUpdateInfo) -> Vec<IndelKey> {
    let offsets = contig_offsets(&m.pivot_contig_lens);
    let mut out: Vec<IndelKey> = Vec::new();
    for v in &m.variants {
        if !v.kind.is_indel() {
            continue;
        }
        if let (Some(s), Some(e)) = (
            flat_pos(&offsets, v.pivot_contig, v.pivot_start),
            flat_pos(&offsets, v.pivot_contig, v.pivot_end),
        ) {
            out.push((s, e, v.kind, v.ref_allele.clone(), v.alt_allele.clone()));
        }
    }
    out.sort();
    out.dedup();
    out
}

/// True when a flat position falls inside one of the merged intervals.
pub fn in_any(merged: &[Iv], pos: u64) -> bool {
    let i = merged.partition_point(|&(s, _)| s <= pos);
    i > 0 && pos < merged[i - 1].1
}

/// One member's contribution to a pairwise comparison, precomputed so that
/// comparing every pair of a cluster is quadratic in *pairs*, not in variants.
#[derive(Debug, Clone, Default)]
pub struct MemberProfile {
    pub id: u32,
    /// SNP alt bases by flat pivot position.
    pub snps: HashMap<u64, u8>,
    /// Sorted, de-duplicated indel keys in flat pivot coordinates.
    pub indels: Vec<IndelKey>,
    /// Merged pivot intervals this member does not cover.
    pub uncovered: Vec<Iv>,
    pub pivot_len: u64,
}

/// Reduce one member's map data to a [`MemberProfile`].
pub fn profile(m: &MapUpdateInfo, id: u32) -> MemberProfile {
    MemberProfile {
        id,
        snps: snp_map(m),
        indels: indel_keys(m),
        uncovered: merge(m.uncovered.clone()),
        pivot_len: pivot_len(m),
    }
}

/// Recompute 6.9 from two members' variant lists.
///
/// `focus` is the weaker of the two members' focus levels; below Focus 2 the
/// variant lists do not exist and the result is [`PairView::unavailable`].
pub fn pair_from_maps(
    a: &MapUpdateInfo,
    b: &MapUpdateInfo,
    a_id: u32,
    b_id: u32,
    focus: u8,
) -> PairView {
    let plen = pivot_len(a).min(pivot_len(b));
    let pa = MemberProfile {
        id: a_id,
        snps: snp_map(a),
        indels: indel_keys(a),
        uncovered: merge(a.uncovered.clone()),
        pivot_len: plen,
    };
    let pb = MemberProfile {
        id: b_id,
        snps: snp_map(b),
        indels: indel_keys(b),
        uncovered: merge(b.uncovered.clone()),
        pivot_len: plen,
    };
    pair_from_profiles(&pa, &pb, focus)
}

/// The same computation over two precomputed profiles.
pub fn pair_from_profiles(a: &MemberProfile, b: &MemberProfile, focus: u8) -> PairView {
    if focus < fmt::FOCUS_EXACT {
        return PairView::unavailable(a.id, b.id, focus);
    }
    let plen = a.pivot_len.min(b.pivot_len);
    let uncovered_union = union_len(&a.uncovered, &b.uncovered);
    let comparable = plen.saturating_sub(uncovered_union);

    let mut positions: Vec<u64> = Vec::with_capacity(a.snps.len() + b.snps.len());
    positions.extend(a.snps.keys().copied());
    positions.extend(b.snps.keys().copied());
    positions.sort_unstable();
    positions.dedup();

    let mut snps = 0u64;
    for p in positions {
        if p >= plen {
            continue;
        }
        if uncovered_union > 0 && (in_any(&a.uncovered, p) || in_any(&b.uncovered, p)) {
            continue;
        }
        // A member with no entry has the pivot's base, and a variant's alt
        // always differs from that base, so the two differ exactly when the
        // entries differ. The pivot's own sequence is not needed.
        if a.snps.get(&p) != b.snps.get(&p) {
            snps += 1;
        }
    }

    let indels = symmetric_difference(&a.indels, &b.indels)
        .into_iter()
        .filter(|k| k.0 < plen && !(in_any(&a.uncovered, k.0) || in_any(&b.uncovered, k.0)))
        .count() as u64;

    let dist = if comparable > 0 {
        Some(snps as f64 / comparable as f64)
    } else {
        None
    };
    PairView {
        a: a.id,
        b: b.id,
        snps,
        indels,
        comparable,
        dist,
        focus,
        approximate: false,
        source: PairSource::Variants,
    }
}

/// Elements of two sorted lists that are in exactly one of them.
pub fn symmetric_difference<T: Ord + Clone>(a: &[T], b: &[T]) -> Vec<T> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            out.push(a[i].clone());
            i += 1;
        } else {
            out.push(b[j].clone());
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Read a pair straight out of the engine's matrix, when the run has finished.
pub fn pair_from_matrix(
    m: &DistanceMatrix,
    ia: usize,
    ib: usize,
    a_id: u32,
    b_id: u32,
    focus: u8,
) -> PairView {
    let dist = if m.comparable[ia][ib] > 0 {
        Some(m.dist[ia][ib])
    } else {
        None
    };
    PairView {
        a: a_id,
        b: b_id,
        snps: m.snps[ia][ib],
        indels: m.indels[ia][ib],
        comparable: m.comparable[ia][ib],
        dist,
        focus,
        approximate: focus < fmt::FOCUS_EXACT,
        source: PairSource::Matrix,
    }
}

/// The variants of one member, in the order the map panel and the export use.
pub fn sorted_variants(m: &MapUpdateInfo) -> Vec<&Variant> {
    let mut v: Vec<&Variant> = m.variants.iter().collect();
    v.sort_by_key(|x| (x.pivot_contig, x.pivot_start, x.pivot_end));
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::PieceView;

    fn variant(kind: VariantKind, contig: u32, start: u64, end: u64, r: &str, a: &str) -> Variant {
        Variant {
            kind,
            pivot_contig: contig,
            pivot_start: start,
            pivot_end: end,
            ref_allele: r.to_string(),
            alt_allele: a.to_string(),
            genome_contig: 0,
            genome_start: start,
            genome_end: end,
            strand: '+',
            flags: 0,
        }
    }

    /// A pivot of two contigs, 100 + 50 bases.
    fn member(name: &str, variants: Vec<Variant>, uncovered: Vec<Iv>) -> MapUpdateInfo {
        MapUpdateInfo {
            cluster: 0,
            genome: 1,
            genome_name: name.to_string(),
            pivot_name: "pivot".to_string(),
            pieces: vec![PieceView {
                pivot_contig: 0,
                pivot_contig_name: "c0".to_string(),
                start: 10,
                end: 20,
                resolved: true,
            }],
            variants,
            uncovered,
            pivot_contigs: vec!["c0".to_string(), "c1".to_string()],
            pivot_contig_lens: vec![100, 50],
        }
    }

    #[test]
    fn contig_offsets_are_a_prefix_sum() {
        assert_eq!(contig_offsets(&[100, 50, 7]), vec![0, 100, 150]);
        assert_eq!(contig_offsets(&[]), Vec::<u64>::new());
        assert_eq!(pivot_len(&member("a", vec![], vec![])), 150);
    }

    #[test]
    fn identical_members_are_zero_apart() {
        let v = vec![variant(VariantKind::Snp, 0, 10, 11, "A", "G")];
        let a = member("a", v.clone(), vec![]);
        let b = member("b", v, vec![]);
        let p = pair_from_maps(&a, &b, 1, 2, 2);
        assert_eq!(p.snps, 0);
        assert_eq!(p.indels, 0);
        assert_eq!(p.comparable, 150);
        assert_eq!(p.dist, Some(0.0));
        assert!(!p.approximate);
        assert_eq!(p.source, PairSource::Variants);
    }

    #[test]
    fn one_shared_snp_and_one_private_snp() {
        let shared = variant(VariantKind::Snp, 0, 10, 11, "A", "G");
        let only_a = variant(VariantKind::Snp, 0, 20, 21, "C", "T");
        let a = member("a", vec![shared.clone(), only_a], vec![]);
        let b = member("b", vec![shared], vec![]);
        let p = pair_from_maps(&a, &b, 1, 2, 2);
        // The shared SNP is the same base on both sides, so it is not a
        // difference; the private one is.
        assert_eq!(p.snps, 1);
        assert!((p.dist.unwrap() - 1.0 / 150.0).abs() < 1e-12);
    }

    #[test]
    fn the_same_position_with_different_alts_counts_once() {
        let a = member(
            "a",
            vec![variant(VariantKind::Snp, 0, 10, 11, "A", "G")],
            vec![],
        );
        let b = member(
            "b",
            vec![variant(VariantKind::Snp, 0, 10, 11, "A", "T")],
            vec![],
        );
        assert_eq!(pair_from_maps(&a, &b, 1, 2, 2).snps, 1);
    }

    #[test]
    fn positions_span_contigs_through_the_offsets() {
        // A SNP at contig 1 position 0 and one at contig 0 position 100 would
        // collide without the offsets; they must not.
        let a = member(
            "a",
            vec![variant(VariantKind::Snp, 1, 0, 1, "A", "G")],
            vec![],
        );
        let b = member(
            "b",
            vec![variant(VariantKind::Snp, 0, 99, 100, "A", "G")],
            vec![],
        );
        assert_eq!(pair_from_maps(&a, &b, 1, 2, 2).snps, 2);
        assert_eq!(snp_map(&a).keys().copied().collect::<Vec<_>>(), vec![100]);
        assert_eq!(snp_map(&b).keys().copied().collect::<Vec<_>>(), vec![99]);
    }

    #[test]
    fn uncovered_positions_are_excluded_from_both_numbers() {
        let a = member(
            "a",
            vec![
                variant(VariantKind::Snp, 0, 10, 11, "A", "G"),
                variant(VariantKind::Deletion, 0, 40, 45, "ACGTG", ""),
            ],
            vec![(50, 60)],
        );
        // b has a private SNP inside a's uncovered stretch, and a private
        // deletion inside its own.
        let b = member(
            "b",
            vec![
                variant(VariantKind::Snp, 0, 55, 56, "C", "T"),
                variant(VariantKind::Deletion, 0, 70, 72, "GG", ""),
            ],
            vec![(70, 72)],
        );
        let p = pair_from_maps(&a, &b, 1, 2, 2);
        assert_eq!(
            p.comparable,
            150 - (60 - 50) - (72 - 70),
            "union of the gaps"
        );
        // The SNP at 10 differs, the one at 55 is not comparable.
        assert_eq!(p.snps, 1);
        // The deletion at 40 is only in a; the one at 70 is inside b's gap.
        assert_eq!(p.indels, 1);
    }

    #[test]
    fn shared_indels_cancel_and_private_ones_count() {
        let shared = variant(VariantKind::Insertion, 0, 30, 30, "", "AAA");
        let only_a = variant(VariantKind::Deletion, 0, 40, 43, "TTT", "");
        let a = member("a", vec![shared.clone(), only_a], vec![]);
        let b = member("b", vec![shared], vec![]);
        assert_eq!(pair_from_maps(&a, &b, 1, 2, 2).indels, 1);
    }

    #[test]
    fn a_pair_below_focus_two_says_it_was_not_measured() {
        let a = member("a", vec![], vec![]);
        let b = member("b", vec![], vec![]);
        let p = pair_from_maps(&a, &b, 1, 2, 1);
        assert!(!p.available());
        assert_eq!(p.source, PairSource::Unavailable);
        assert!(p.approximate, "an unmeasured pair must never look exact");
        assert_eq!(p.snps_text(), "not measured yet");
        assert_eq!(p.dist_text(), "n/a");
    }

    #[test]
    fn nothing_comparable_leaves_the_distance_unknown() {
        let a = member("a", vec![], vec![(0, 150)]);
        let b = member("b", vec![], vec![]);
        let p = pair_from_maps(&a, &b, 1, 2, 2);
        assert_eq!(p.comparable, 0);
        assert_eq!(p.dist, None);
        assert_eq!(p.dist_text(), "n/a");
    }

    #[test]
    fn in_any_handles_empty_and_boundary_positions() {
        let m = merge(vec![(10, 20), (30, 35)]);
        assert!(!in_any(&m, 9));
        assert!(in_any(&m, 10));
        assert!(in_any(&m, 19));
        assert!(!in_any(&m, 20), "half-open");
        assert!(in_any(&m, 34));
        assert!(!in_any(&[], 5));
    }

    #[test]
    fn symmetric_difference_is_order_independent() {
        let a = vec![1, 3, 5];
        let b = vec![3, 4];
        assert_eq!(symmetric_difference(&a, &b), vec![1, 4, 5]);
        assert_eq!(symmetric_difference(&b, &a), vec![1, 4, 5]);
        assert!(symmetric_difference(&a, &a).is_empty());
        assert_eq!(symmetric_difference(&a, &[]), vec![1, 3, 5]);
    }

    #[test]
    fn the_matrix_path_reports_its_own_source() {
        let m = DistanceMatrix {
            ids: vec![7, 9],
            snps: vec![vec![0, 4], vec![4, 0]],
            indels: vec![vec![0, 1], vec![1, 0]],
            comparable: vec![vec![1000, 990], vec![990, 1000]],
            dist: vec![vec![0.0, 4.0 / 990.0], vec![4.0 / 990.0, 0.0]],
            estimates: vec![vec![0, 0], vec![0, 0]],
            basis: infinity_core::profile::TreeBasis::Snps,
        };
        let p = pair_from_matrix(&m, 0, 1, 7, 9, 2);
        assert_eq!(p.source, PairSource::Matrix);
        assert_eq!(p.snps, 4);
        assert_eq!(p.indels, 1);
        assert!(!p.approximate);
        let p1 = pair_from_matrix(&m, 0, 1, 7, 9, 1);
        assert!(p1.approximate);
        assert!(p1.dist_text().starts_with('≈'));
    }

    #[test]
    fn precomputed_profiles_agree_with_the_direct_path() {
        let a = member(
            "a",
            vec![
                variant(VariantKind::Snp, 0, 10, 11, "A", "G"),
                variant(VariantKind::Deletion, 0, 40, 45, "ACGTG", ""),
            ],
            vec![(50, 60)],
        );
        let b = member(
            "b",
            vec![
                variant(VariantKind::Snp, 0, 10, 11, "A", "G"),
                variant(VariantKind::Insertion, 0, 41, 41, "", "TT"),
            ],
            vec![(70, 72)],
        );
        let direct = pair_from_maps(&a, &b, 1, 2, 2);
        let plen = pivot_len(&a).min(pivot_len(&b));
        let mut pa = profile(&a, 1);
        let mut pb = profile(&b, 2);
        pa.pivot_len = plen;
        pb.pivot_len = plen;
        assert_eq!(pair_from_profiles(&pa, &pb, 2), direct);
        assert!(!pair_from_profiles(&pa, &pb, 1).available());
    }

    #[test]
    fn variants_come_out_in_pivot_order_for_the_map_and_the_export() {
        let m = member(
            "a",
            vec![
                variant(VariantKind::Snp, 1, 0, 1, "A", "G"),
                variant(VariantKind::Snp, 0, 50, 51, "A", "G"),
                variant(VariantKind::Snp, 0, 10, 11, "A", "G"),
            ],
            vec![],
        );
        let s = sorted_variants(&m);
        let keys: Vec<(u32, u64)> = s.iter().map(|v| (v.pivot_contig, v.pivot_start)).collect();
        assert_eq!(keys, vec![(0, 10), (0, 50), (1, 0)]);
    }
}
