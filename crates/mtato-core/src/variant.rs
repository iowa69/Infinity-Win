//! Variants: the `diff` procedure of 6.7, plus left-alignment.
//!
//! Coordinates inside this module are piece-local, 0-based and half-open. The
//! caller translates them into pivot coordinates. Output files use 1-based
//! positions (section 5), and that conversion happens once, in `output.rs`.

use serde::{Deserialize, Serialize};

use crate::align::{global_align, AlignScores, Alignment, Op};

/// What kind of difference a variant is.
///
/// `Ord` is derived so that variants can be sorted into a canonical order,
/// which is what makes output files byte-identical across runs (section 0
/// rule 3). The order of the variants themselves is not meaningful.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum VariantKind {
    Snp,
    Insertion,
    Deletion,
    Complex,
}

impl VariantKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            VariantKind::Snp => "snp",
            VariantKind::Insertion => "ins",
            VariantKind::Deletion => "del",
            VariantKind::Complex => "complex",
        }
    }

    /// True for the kinds that count towards `indels(A,B)` in 6.9.
    pub fn is_indel(&self) -> bool {
        matches!(self, VariantKind::Insertion | VariantKind::Deletion)
    }
}

/// Variant flags (sections 6.7, 6.8 and 15).
pub mod flags {
    /// The interval is covered by two or more places of the member.
    pub const DUP: u32 = 1 << 0;
    /// The block or anchor matched a repeat in the pivot.
    pub const REPEAT: u32 = 1 << 1;
    /// The piece was reported as COMPLEX because identity was too low.
    pub const LOW_IDENTITY: u32 = 1 << 2;
    /// Found by end-of-contig extension rather than between two anchors.
    pub const END_EXT: u32 = 1 << 3;
    /// The alignment had to give up (band too narrow); reported as COMPLEX.
    pub const ALIGN_FAILED: u32 = 1 << 4;

    pub fn names(f: u32) -> Vec<&'static str> {
        let mut v = Vec::new();
        if f & DUP != 0 {
            v.push("dup");
        }
        if f & REPEAT != 0 {
            v.push("repeat");
        }
        if f & LOW_IDENTITY != 0 {
            v.push("low_identity");
        }
        if f & END_EXT != 0 {
            v.push("end_ext");
        }
        if f & ALIGN_FAILED != 0 {
            v.push("align_failed");
        }
        v
    }
}

/// A difference in piece-local coordinates, before translation to the pivot.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RawVariant {
    pub kind: VariantKind,
    /// Half-open interval of the piece reference that the variant replaces.
    /// Empty (`start == end`) for an insertion.
    pub start: usize,
    pub end: usize,
    pub ref_allele: Vec<u8>,
    pub alt_allele: Vec<u8>,
    pub flags: u32,
}

/// A variant in pivot coordinates: one difference between a member and its
/// pivot (section 5).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Variant {
    pub kind: VariantKind,
    pub pivot_contig: u32,
    /// 0-based, half-open, on the pivot's forward strand.
    pub pivot_start: u64,
    pub pivot_end: u64,
    pub ref_allele: String,
    pub alt_allele: String,
    /// Member contig and position, in member forward-strand coordinates.
    pub genome_contig: u32,
    pub genome_start: u64,
    pub genome_end: u64,
    /// Orientation of the member relative to the pivot.
    pub strand: char,
    pub flags: u32,
}

impl Variant {
    /// Length on the pivot. Zero for an insertion.
    pub fn pivot_len(&self) -> u64 {
        self.pivot_end.saturating_sub(self.pivot_start)
    }

    /// A key that identifies this variant for set comparison in 6.9. Two
    /// members share an indel when the keys are equal.
    pub fn identity_key(&self) -> (u32, u64, u64, VariantKind, String, String) {
        (
            self.pivot_contig,
            self.pivot_start,
            self.pivot_end,
            self.kind,
            self.ref_allele.clone(),
            self.alt_allele.clone(),
        )
    }

    pub fn flags_string(&self) -> String {
        flags::names(self.flags).join(",")
    }
}

/// The subset of parameters that `diff` needs (6.7).
#[derive(Debug, Clone, Copy)]
pub struct DiffParams {
    pub dp_max_len: usize,
    pub min_identity: f64,
    pub scores: AlignScores,
}

impl From<&crate::params::Parameters> for DiffParams {
    fn from(p: &crate::params::Parameters) -> Self {
        Self {
            dp_max_len: p.dp_max_len,
            min_identity: p.min_identity,
            scores: AlignScores::default(),
        }
    }
}

/// The outcome of resolving one difference piece.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct DiffResult {
    pub variants: Vec<RawVariant>,
    /// Piece-reference offsets where a base was N on either side. Those
    /// positions are never SNPs and count as not covered (6.7).
    pub n_positions: Vec<usize>,
    /// Set when the piece could not be aligned and became one COMPLEX variant.
    pub align_failed: bool,
    /// True when a dynamic-programming alignment was actually run, as opposed
    /// to the cheap trim / equal-length paths. Section 11 counts these, so the
    /// number must come from what happened rather than from an estimate.
    pub used_dp: bool,
}

/// Resolve a difference piece exactly, as specified in 6.7.
///
/// `ref_s` is the pivot piece and `alt_s` the member piece (already
/// reverse-complemented when the member runs the other way). The two begin and
/// end with the same `k` bases, which the prefix/suffix trim removes.
pub fn diff(ref_s: &[u8], alt_s: &[u8], p: &DiffParams) -> DiffResult {
    let lr = ref_s.len();
    let la = alt_s.len();

    // Longest common prefix.
    let mut pre = 0usize;
    while pre < lr && pre < la && ref_s[pre] == alt_s[pre] {
        pre += 1;
    }
    // Longest common suffix of what remains.
    let mut suf = 0usize;
    while suf < lr - pre && suf < la - pre && ref_s[lr - 1 - suf] == alt_s[la - 1 - suf] {
        suf += 1;
    }
    let r = &ref_s[pre..lr - suf];
    let a = &alt_s[pre..la - suf];

    if r.is_empty() && a.is_empty() {
        return DiffResult::default();
    }

    // Pure insertion or deletion after trimming.
    if r.is_empty() {
        let mut alt = a.to_vec();
        let pos = left_align_insertion(ref_s, pre, &mut alt);
        return DiffResult {
            variants: vec![RawVariant {
                kind: VariantKind::Insertion,
                start: pos,
                end: pos,
                ref_allele: Vec::new(),
                alt_allele: alt,
                flags: 0,
            }],
            n_positions: n_offsets(r, a, pre),
            align_failed: false,
            used_dp: false,
        };
    }
    if a.is_empty() {
        let pos = left_align_deletion(ref_s, pre, r.len());
        return DiffResult {
            variants: vec![RawVariant {
                kind: VariantKind::Deletion,
                start: pos,
                end: pos + r.len(),
                ref_allele: r.to_vec(),
                alt_allele: Vec::new(),
                flags: 0,
            }],
            n_positions: n_offsets(r, a, pre),
            align_failed: false,
            used_dp: false,
        };
    }

    // Same length: cheap path for the overwhelmingly common case.
    if r.len() == a.len() {
        let mismatches: Vec<usize> = (0..r.len()).filter(|&i| r[i] != a[i]).collect();
        let ratio = mismatches.len() as f64 / r.len() as f64;
        if mismatches.len() <= 3 || ratio <= 0.10 {
            let (variants, n_positions) = snps_from_positions(ref_s, alt_s, r, a, pre, &mismatches);
            return DiffResult {
                variants,
                n_positions,
                align_failed: false,
                used_dp: false,
            };
        }
    }

    if r.len().max(a.len()) <= p.dp_max_len {
        match global_align(r, a, &p.scores) {
            Some(al) if al.identity() >= p.min_identity => {
                let (variants, n_positions) = from_alignment(ref_s, alt_s, &al, pre);
                return DiffResult {
                    variants,
                    n_positions,
                    align_failed: false,
                    used_dp: true,
                };
            }
            Some(_) => {
                // Identity below the floor: one COMPLEX variant for r (6.7).
                return DiffResult {
                    variants: vec![complex(r, a, pre, flags::LOW_IDENTITY)],
                    n_positions: n_offsets(r, a, pre),
                    align_failed: false,
                    used_dp: true,
                };
            }
            None => {
                return DiffResult {
                    variants: vec![complex(r, a, pre, flags::ALIGN_FAILED)],
                    n_positions: n_offsets(r, a, pre),
                    align_failed: true,
                    used_dp: true,
                };
            }
        }
    }

    DiffResult {
        variants: vec![complex(r, a, pre, flags::LOW_IDENTITY)],
        n_positions: n_offsets(r, a, pre),
        align_failed: false,
        used_dp: false,
    }
}

fn complex(r: &[u8], a: &[u8], pre: usize, flag: u32) -> RawVariant {
    RawVariant {
        kind: VariantKind::Complex,
        start: pre,
        end: pre + r.len(),
        ref_allele: r.to_vec(),
        alt_allele: a.to_vec(),
        flags: flag,
    }
}

/// Offsets (in full-piece reference coordinates) where either side is N.
fn n_offsets(r: &[u8], a: &[u8], pre: usize) -> Vec<usize> {
    let mut v = Vec::new();
    for (i, &b) in r.iter().enumerate() {
        if b == b'N' {
            v.push(pre + i);
        }
    }
    // Alt-side Ns map onto the reference position they sit opposite, which is
    // only well defined when the two sides are the same length; otherwise the
    // whole trimmed region is already reported through the reference scan.
    if a.len() == r.len() {
        for (i, &b) in a.iter().enumerate() {
            if b == b'N' && r[i] != b'N' {
                v.push(pre + i);
            }
        }
    }
    v.sort_unstable();
    v.dedup();
    v
}

fn snps_from_positions(
    ref_s: &[u8],
    alt_s: &[u8],
    r: &[u8],
    a: &[u8],
    pre: usize,
    mismatches: &[usize],
) -> (Vec<RawVariant>, Vec<usize>) {
    let mut variants = Vec::with_capacity(mismatches.len());
    let mut n_positions = Vec::new();
    for &i in mismatches {
        let (rb, ab) = (r[i], a[i]);
        if rb == b'N' || ab == b'N' {
            n_positions.push(pre + i);
            continue;
        }
        variants.push(RawVariant {
            kind: VariantKind::Snp,
            start: pre + i,
            end: pre + i + 1,
            ref_allele: vec![rb],
            alt_allele: vec![ab],
            flags: 0,
        });
    }
    // Ns that are equal on both sides are not mismatches, so scan for them too.
    for i in 0..r.len() {
        if r[i] == b'N' && a[i] == b'N' {
            n_positions.push(pre + i);
        }
    }
    let _ = (ref_s, alt_s);
    n_positions.sort_unstable();
    n_positions.dedup();
    (variants, n_positions)
}

/// Read SNPs, insertions and deletions out of an alignment, left-aligning
/// every indel. Coordinates are returned in full-piece reference space.
fn from_alignment(
    ref_s: &[u8],
    alt_s: &[u8],
    al: &Alignment,
    pre: usize,
) -> (Vec<RawVariant>, Vec<usize>) {
    let mut variants: Vec<RawVariant> = Vec::new();
    let mut n_positions: Vec<usize> = Vec::new();
    // `ri`/`ai` walk the trimmed remainders; the piece ref/alt start at `pre`.
    let r = &ref_s[pre..];
    let a = &alt_s[pre..];
    let mut ri = 0usize;
    let mut ai = 0usize;
    let mut ops = al.ops.iter().peekable();

    while let Some(op) = ops.next() {
        match op {
            Op::Match => {
                ri += 1;
                ai += 1;
            }
            Op::Mismatch => {
                let (rb, ab) = (r[ri], a[ai]);
                if rb == b'N' || ab == b'N' {
                    n_positions.push(pre + ri);
                } else {
                    variants.push(RawVariant {
                        kind: VariantKind::Snp,
                        start: pre + ri,
                        end: pre + ri + 1,
                        ref_allele: vec![rb],
                        alt_allele: vec![ab],
                        flags: 0,
                    });
                }
                ri += 1;
                ai += 1;
            }
            Op::Insertion => {
                let start_ai = ai;
                let start_ri = ri;
                while ops.peek() == Some(&&Op::Insertion) {
                    ops.next();
                    ai += 1;
                }
                ai += 1;
                let mut ins = a[start_ai..ai].to_vec();
                let pos = left_align_insertion(ref_s, pre + start_ri, &mut ins);
                variants.push(RawVariant {
                    kind: VariantKind::Insertion,
                    start: pos,
                    end: pos,
                    ref_allele: Vec::new(),
                    alt_allele: ins,
                    flags: 0,
                });
            }
            Op::Deletion => {
                let start_ri = ri;
                while ops.peek() == Some(&&Op::Deletion) {
                    ops.next();
                    ri += 1;
                }
                ri += 1;
                let len = ri - start_ri;
                let pos = left_align_deletion(ref_s, pre + start_ri, len);
                variants.push(RawVariant {
                    kind: VariantKind::Deletion,
                    start: pos,
                    end: pos + len,
                    ref_allele: ref_s[pos..pos + len].to_vec(),
                    alt_allele: Vec::new(),
                    flags: 0,
                });
            }
        }
    }
    variants.sort_by_key(|v| (v.start, v.end));
    (variants, n_positions)
}

/// Shift an insertion as far left as it will go (6.7, standard left-alignment).
///
/// `pos` is the piece-reference offset the insertion sits before; `ins` is
/// rotated in place so that it still describes the same sequence change.
pub fn left_align_insertion(ref_s: &[u8], mut pos: usize, ins: &mut Vec<u8>) -> usize {
    while pos > 0 {
        let last = *ins.last().unwrap_or(&0);
        if ref_s[pos - 1] == last {
            ins.pop();
            ins.insert(0, ref_s[pos - 1]);
            pos -= 1;
        } else {
            break;
        }
    }
    pos
}

/// Shift a deletion as far left as it will go.
pub fn left_align_deletion(ref_s: &[u8], mut pos: usize, len: usize) -> usize {
    if len == 0 {
        return pos;
    }
    while pos > 0 && pos + len <= ref_s.len() && ref_s[pos - 1] == ref_s[pos + len - 1] {
        pos -= 1;
    }
    pos
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;

    fn dp() -> DiffParams {
        DiffParams::from(&Parameters::default())
    }

    fn v(x: &[u8]) -> String {
        String::from_utf8_lossy(x).into_owned()
    }

    #[test]
    fn identical_pieces_produce_nothing() {
        let r = b"ACGTACGTACGT";
        assert!(diff(r, r, &dp()).variants.is_empty());
    }

    #[test]
    fn one_snp_is_one_snp_with_the_right_bases_and_position() {
        // Property 3, at the piece level: the two pieces differ at exactly one
        // position, index 6, T -> A. Everything else is identical, so the
        // prefix/suffix trim must leave a single base on each side.
        let r = b"ACGTTTTTGCA";
        let a = b"ACGTTTATGCA";
        assert_eq!(r.len(), a.len());
        let d = diff(r, a, &dp());
        assert_eq!(d.variants.len(), 1);
        let x = &d.variants[0];
        assert_eq!(x.kind, VariantKind::Snp);
        assert_eq!(x.start, 6);
        assert_eq!(x.end, 7);
        assert_eq!(v(&x.ref_allele), "T");
        assert_eq!(v(&x.alt_allele), "A");
    }

    #[test]
    fn up_to_three_mismatches_stay_snps_whatever_the_length() {
        let mut r = vec![b'A'; 100];
        let mut a = r.clone();
        r[10] = b'C';
        a[10] = b'G';
        r[50] = b'C';
        a[50] = b'G';
        r[90] = b'C';
        a[90] = b'G';
        let d = diff(&r, &a, &dp());
        assert_eq!(d.variants.len(), 3);
        assert!(d.variants.iter().all(|x| x.kind == VariantKind::Snp));
        assert_eq!(
            d.variants.iter().map(|x| x.start).collect::<Vec<_>>(),
            vec![10, 50, 90]
        );
    }

    #[test]
    fn a_pure_insertion_is_one_insertion() {
        let r = b"ACGTACGT";
        let a = b"ACGTTTTACGT";
        let d = diff(r, a, &dp());
        assert_eq!(d.variants.len(), 1);
        assert_eq!(d.variants[0].kind, VariantKind::Insertion);
        assert_eq!(v(&d.variants[0].alt_allele), "TTT");
        assert!(d.variants[0].ref_allele.is_empty());
    }

    #[test]
    fn a_pure_deletion_is_one_deletion() {
        let r = b"ACGTTTTACGT";
        let a = b"ACGTACGT";
        let d = diff(r, a, &dp());
        assert_eq!(d.variants.len(), 1);
        assert_eq!(d.variants[0].kind, VariantKind::Deletion);
        assert_eq!(v(&d.variants[0].ref_allele), "TTT");
        assert!(d.variants[0].alt_allele.is_empty());
    }

    #[test]
    fn indels_are_left_aligned_in_a_homopolymer() {
        // ref A A G  vs  alt A A A G: the extra A must land at the leftmost
        // position of the run, so the same event always gets the same
        // coordinates no matter which side it is seen from.
        let r = b"TTTTAAGTTTT";
        let a = b"TTTTAAAGTTTT";
        let d1 = diff(r, a, &dp());
        assert_eq!(d1.variants.len(), 1);
        assert_eq!(d1.variants[0].kind, VariantKind::Insertion);
        assert_eq!(d1.variants[0].start, 4, "shifted to the start of the A run");

        // Symmetrically for the deletion.
        let d2 = diff(a, r, &dp());
        assert_eq!(d2.variants.len(), 1);
        assert_eq!(d2.variants[0].kind, VariantKind::Deletion);
        assert_eq!(d2.variants[0].start, 4);
    }

    #[test]
    fn left_alignment_helpers_are_exact() {
        let r = b"CAAAAT";
        let mut ins = b"A".to_vec();
        assert_eq!(left_align_insertion(r, 4, &mut ins), 1);
        assert_eq!(ins, b"A".to_vec());
        let r2 = b"ACGTACGT";
        assert_eq!(left_align_deletion(r2, 3, 2), 3, "no shift possible");
        let r3 = b"ATATATX";
        assert_eq!(
            left_align_deletion(r3, 4, 2),
            0,
            "shifts through the repeat"
        );
    }

    #[test]
    fn an_n_position_is_never_a_snp() {
        let r = b"ACGTNCGT";
        let a = b"ACGTACGT";
        let d = diff(r, a, &dp());
        assert!(d.variants.is_empty(), "N vs A is not a SNP");
        assert_eq!(d.n_positions, vec![4]);
    }

    #[test]
    fn a_divergent_piece_becomes_one_complex_variant() {
        let r = b"AAAAACCCCCGGGGGTTTTTAAAA";
        let a = b"AAAAGGGGGCCCCCTTTTTGGGGG";
        let p = DiffParams {
            min_identity: 0.99,
            ..dp()
        };
        let d = diff(r, a, &p);
        assert_eq!(d.variants.len(), 1);
        assert_eq!(d.variants[0].kind, VariantKind::Complex);
    }

    #[test]
    fn a_piece_too_long_for_dp_becomes_complex() {
        let r = vec![b'A'; 200];
        let mut a = vec![b'C'; 200];
        a[..5].copy_from_slice(&r[..5]);
        a[195..].copy_from_slice(&r[195..]);
        let p = DiffParams {
            dp_max_len: 50,
            ..dp()
        };
        let d = diff(&r, &a, &p);
        assert_eq!(d.variants.len(), 1);
        assert_eq!(d.variants[0].kind, VariantKind::Complex);
    }

    #[test]
    fn mixed_snp_and_indel_are_both_found() {
        let r = b"TTTTACGTACGTACGTTTTT";
        let a = b"TTTTACGAACGTACCGTTTTT";
        let d = diff(r, a, &dp());
        let kinds: Vec<_> = d.variants.iter().map(|x| x.kind).collect();
        assert!(kinds.contains(&VariantKind::Snp), "{kinds:?}");
        assert!(kinds.iter().any(|k| k.is_indel()), "{kinds:?}");
        // variants come out in position order
        let starts: Vec<_> = d.variants.iter().map(|x| x.start).collect();
        let mut sorted = starts.clone();
        sorted.sort_unstable();
        assert_eq!(starts, sorted);
    }

    #[test]
    fn identity_key_separates_different_indels_at_the_same_place() {
        let a = Variant {
            kind: VariantKind::Insertion,
            pivot_contig: 0,
            pivot_start: 10,
            pivot_end: 10,
            ref_allele: String::new(),
            alt_allele: "A".into(),
            genome_contig: 0,
            genome_start: 10,
            genome_end: 10,
            strand: '+',
            flags: 0,
        };
        let mut b = a.clone();
        b.alt_allele = "C".into();
        assert_ne!(a.identity_key(), b.identity_key());
        assert_eq!(a.identity_key(), a.clone().identity_key());
    }

    #[test]
    fn flag_names_round_trip() {
        let f = flags::DUP | flags::REPEAT;
        let names = flags::names(f);
        assert_eq!(names, vec!["dup", "repeat"]);
        assert!(flags::names(0).is_empty());
    }
}
