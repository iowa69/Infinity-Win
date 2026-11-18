//! Alignment primitives: banded affine-gap global alignment and x-drop
//! ungapped extension (6.6 step D and 6.7).
//!
//! Written in-house rather than pulled from a crate, for three reasons: the
//! scoring scheme is fixed by the outline and must match exactly, the result
//! must be deterministic and independent of thread count, and a pure-Rust
//! aligner keeps the Windows build free of extra toolchain requirements
//! (section 12).

use serde::{Deserialize, Serialize};

/// Affine gap scores. The outline fixes them: match +1, mismatch -3, gap open
/// -5, gap extend -1.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignScores {
    pub r#match: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for AlignScores {
    fn default() -> Self {
        Self {
            r#match: 1,
            mismatch: -3,
            gap_open: -5,
            gap_extend: -1,
        }
    }
}

/// One step of an alignment, always described relative to the reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Op {
    /// Aligned pair of equal bases.
    Match,
    /// Aligned pair of different bases (a substitution).
    Mismatch,
    /// Base present in `alt` and absent from `ref`.
    Insertion,
    /// Base present in `ref` and absent from `alt`.
    Deletion,
}

/// A finished alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub ops: Vec<Op>,
    pub score: i32,
    /// Number of `Match` ops.
    pub matches: usize,
    /// Number of aligned columns, i.e. `ops.len()`.
    pub columns: usize,
}

impl Alignment {
    /// Identity over aligned columns, in `[0, 1]`.
    pub fn identity(&self) -> f64 {
        if self.columns == 0 {
            return 1.0;
        }
        self.matches as f64 / self.columns as f64
    }
}

const NEG: i32 = i32::MIN / 4;

/// Half-width of the band, in cells, around the diagonal.
///
/// Pieces are stretches between two shared anchors of two nearly identical
/// genomes, so the optimal path stays close to the diagonal; the band is what
/// keeps a 5 kb piece (`dp_max_len`) at a few megabytes instead of the 100 MB
/// a full Gotoh matrix would need. The width grows with the length difference,
/// which is the minimum any correct alignment must be able to drift.
fn band_width(len_r: usize, len_a: usize) -> usize {
    let diff = len_r.abs_diff(len_a);
    let longest = len_r.max(len_a);
    let pad = 96 + longest / 32;
    diff + pad
}

/// Global affine-gap alignment (Gotoh) inside a band.
///
/// Returns `None` when the band is too narrow to carry a path from the top
/// left to the bottom right corner. The caller turns that into a COMPLEX
/// variant, which is the honest answer: an approximate alignment would be
/// reported as exact, and section 0 rule 8 forbids that.
pub fn global_align(ref_s: &[u8], alt: &[u8], s: &AlignScores) -> Option<Alignment> {
    let n = ref_s.len();
    let m = alt_s_len(alt);
    if n == 0 && m == 0 {
        return Some(Alignment {
            ops: Vec::new(),
            score: 0,
            matches: 0,
            columns: 0,
        });
    }
    if n == 0 {
        return Some(Alignment {
            ops: vec![Op::Insertion; m],
            score: s.gap_open + s.gap_extend * (m as i32 - 1),
            matches: 0,
            columns: m,
        });
    }
    if m == 0 {
        return Some(Alignment {
            ops: vec![Op::Deletion; n],
            score: s.gap_open + s.gap_extend * (n as i32 - 1),
            matches: 0,
            columns: n,
        });
    }

    // Band limits: for row i, columns j in [lo_i, hi_i].
    let w = band_width(n, m) as i64;
    let mut lo = vec![0i64; n + 1];
    let mut hi = vec![0i64; n + 1];
    // The diagonal that the corner (n, m) sits on, expressed as j - i*(m/n).
    for i in 0..=n {
        let centre = (i as i64) * (m as i64) / (n as i64);
        let l = (centre - w).max(0);
        let r = (centre + w).min(m as i64);
        lo[i] = l;
        hi[i] = r;
    }
    // Row 0 must start at column 0 and row n must end at column m.
    lo[0] = 0;
    hi[n] = m as i64;

    // Three score matrices, stored banded. Directions are packed into a u8.
    let mut width: Vec<usize> = Vec::with_capacity(n + 1);
    for i in 0..=n {
        width.push((hi[i] - lo[i] + 1).max(0) as usize);
    }
    let mut off = vec![0usize; n + 2];
    for i in 0..=n {
        off[i + 1] = off[i] + width[i];
    }
    let cells = off[n + 1];
    if cells == 0 || cells > 200_000_000 {
        return None;
    }

    // M: alignment ending with a match/mismatch; I: ending with an insertion
    // (base from alt); D: ending with a deletion (base from ref).
    let mut mm = vec![NEG; cells];
    let mut mi = vec![NEG; cells];
    let mut md = vec![NEG; cells];
    // Traceback: 0 = none, 1 = from M, 2 = from I, 3 = from D, per matrix.
    let mut tm = vec![0u8; cells];
    let mut ti = vec![0u8; cells];
    let mut td = vec![0u8; cells];

    let idx = |i: usize, j: i64| -> Option<usize> {
        if j < lo[i] || j > hi[i] {
            return None;
        }
        Some(off[i] + (j - lo[i]) as usize)
    };

    // Initialisation.
    mm[0] = 0;
    if let Some(c) = idx(0, 1) {
        // first insertion
        let mut prev = 0i32;
        for j in 1..=hi[0] {
            if let Some(cj) = idx(0, j) {
                let v = if j == 1 {
                    s.gap_open
                } else {
                    prev.saturating_add(s.gap_extend)
                };
                mi[cj] = v;
                ti[cj] = 2;
                prev = v;
                let _ = c;
            }
        }
    }

    for i in 1..=n {
        let rb = ref_s[i - 1];
        for j in lo[i].max(1)..=hi[i] {
            let ju = j as usize;
            let cur = idx(i, j).unwrap();
            let diag = idx(i - 1, j - 1);
            let up = idx(i - 1, j);
            let left = idx(i, j - 1);

            // Match / mismatch from any of the three previous states.
            if let Some(d) = diag {
                let sc = if rb == alt[ju - 1] {
                    s.r#match
                } else {
                    s.mismatch
                };
                let best = mm[d].max(mi[d]).max(md[d]);
                if best > NEG / 2 {
                    let v = best + sc;
                    if v > mm[cur] {
                        mm[cur] = v;
                        tm[cur] = if mm[d] >= mi[d] && mm[d] >= md[d] {
                            1
                        } else if mi[d] >= md[d] {
                            2
                        } else {
                            3
                        };
                    }
                }
            }
            // Deletion: consume a ref base (move down).
            if let Some(u) = up {
                let cand_d = md[u].saturating_add(s.gap_extend);
                let cand_o = mm[u].max(mi[u]).saturating_add(s.gap_open);
                let (v, from) = if cand_d >= cand_o {
                    (cand_d, 3)
                } else {
                    (cand_o, if mm[u] >= mi[u] { 1 } else { 2 })
                };
                if v > NEG / 2 && v > md[cur] {
                    md[cur] = v;
                    td[cur] = from;
                }
            }
            // Insertion: consume an alt base (move right).
            if let Some(l) = left {
                let cand_i = mi[l].saturating_add(s.gap_extend);
                let cand_o = mm[l].max(md[l]).saturating_add(s.gap_open);
                let (v, from) = if cand_i >= cand_o {
                    (cand_i, 2)
                } else {
                    (cand_o, if mm[l] >= md[l] { 1 } else { 3 })
                };
                if v > NEG / 2 && v > mi[cur] {
                    mi[cur] = v;
                    ti[cur] = from;
                }
            }
            // Suppress impossible first-column cells.
            if j == 0 {
                mi[cur] = NEG;
                md[cur] = NEG;
            }
        }
        // Column 0 of row i: only deletions can reach it.
        if lo[i] == 0 {
            let cur = off[i];
            if let Some(u) = idx(i - 1, 0) {
                let cand_d = md[u].saturating_add(s.gap_extend);
                let cand_o = mm[u].saturating_add(s.gap_open);
                let (v, from) = if cand_d >= cand_o {
                    (cand_d, 3)
                } else {
                    (cand_o, 1)
                };
                if v > NEG / 2 {
                    md[cur] = v;
                    td[cur] = from;
                }
                mm[cur] = NEG;
                mi[cur] = NEG;
            }
        }
    }

    let last = idx(n, m as i64)?;
    // Pick the best ending state.
    let (mut state, score) = if mm[last] >= mi[last] && mm[last] >= md[last] {
        (1u8, mm[last])
    } else if mi[last] >= md[last] {
        (2u8, mi[last])
    } else {
        (3u8, md[last])
    };
    if score <= NEG / 2 {
        return None; // no path inside the band
    }

    // Traceback.
    let mut ops: Vec<Op> = Vec::new();
    let mut i = n as i64;
    let mut j = m as i64;
    let mut guard = 0usize;
    let limit = (n + m) * 4 + 16;
    while i > 0 || j > 0 {
        guard += 1;
        if guard > limit {
            return None;
        }
        let cur = idx(i as usize, j)?;
        let (from, op): (u8, Op) = match state {
            1 => {
                let op = if ref_s[i as usize - 1] == alt[j as usize - 1] {
                    Op::Match
                } else {
                    Op::Mismatch
                };
                (tm[cur], op)
            }
            2 => (ti[cur], Op::Insertion),
            3 => (td[cur], Op::Deletion),
            _ => return None,
        };
        ops.push(op);
        state = from;
        match op {
            Op::Match | Op::Mismatch => {
                i -= 1;
                j -= 1;
            }
            Op::Insertion => j -= 1,
            Op::Deletion => i -= 1,
        }
        if state == 0 && (i > 0 || j > 0) {
            return None;
        }
    }
    ops.reverse();
    let matches = ops.iter().filter(|o| **o == Op::Match).count();
    Some(Alignment {
        columns: ops.len(),
        score,
        matches,
        ops,
    })
}

fn alt_s_len(alt: &[u8]) -> usize {
    alt.len()
}

/// Extend in any combination of directions, ungapped.
///
/// `r_dir` and `a_dir` are `+1` or `-1`: the walk compares
/// `ref_s[r_start + i*r_dir]` with `alt[a_start + i*a_dir]` for `i = 0, 1, ...`
/// and stops at the first index out of range, or when the running score has
/// fallen `xdrop` below its best. `mismatches` holds the step indices, not the
/// coordinates, so the caller can convert them whichever way it needs.
///
/// This generalises [`xdrop_forward`] and [`xdrop_backward`], and it is what
/// 6.6 step D needs when the member runs opposite to the pivot: there the
/// member walks backwards while the pivot walks forwards.
///
/// `complement_alt` compares the complement of the `alt` base rather than the
/// base itself. That is required when the two strands are antiparallel: a
/// *matching* position then holds complementary bases, and comparing them
/// directly would report a mismatch at every single position and stop the
/// extension after five steps. Only the comparison is complemented — the caller
/// still reports the member's actual base as the variant allele.
pub fn xdrop_walk(
    ref_s: &[u8],
    alt: &[u8],
    r_start: i64,
    a_start: i64,
    r_dir: i64,
    a_dir: i64,
    complement_alt: bool,
    scores: &AlignScores,
    xdrop: i32,
) -> XdropResult {
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_len = 0usize;
    let mut best_mismatches = 0usize;
    let mut mismatches: Vec<usize> = Vec::new();
    let mut hit_end = true;
    let mut i = 0i64;
    loop {
        let rp = r_start + i * r_dir;
        let ap = a_start + i * a_dir;
        if rp < 0 || ap < 0 || rp >= ref_s.len() as i64 || ap >= alt.len() as i64 {
            break;
        }
        let ab = alt[ap as usize];
        let ab = if complement_alt {
            crate::seq::complement_byte(ab)
        } else {
            ab
        };
        let eq = ref_s[rp as usize] == ab;
        score += if eq { scores.r#match } else { scores.mismatch };
        if !eq {
            mismatches.push(i as usize);
        }
        i += 1;
        if score > best {
            best = score;
            best_len = i as usize;
            best_mismatches = mismatches.len();
        }
        if best - score >= xdrop {
            hit_end = false;
            break;
        }
    }
    mismatches.truncate(best_mismatches);
    XdropResult {
        len: best_len,
        mismatches,
        hit_end,
    }
}

/// The pivot interval covered by a walk of `kept` steps from `start` in
/// direction `dir`.
pub fn walked_range(start: i64, dir: i64, kept: usize) -> (i64, i64) {
    if kept == 0 {
        return (start, start);
    }
    if dir > 0 {
        (start, start + kept as i64)
    } else {
        (start - kept as i64 + 1, start + 1)
    }
}

/// Result of an ungapped x-drop extension.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct XdropResult {
    /// Bases kept, counted from the start position outward.
    pub len: usize,
    /// Offsets within the kept part where the two sequences differ.
    pub mismatches: Vec<usize>,
    /// True when the extension stopped because a sequence ran out.
    pub hit_end: bool,
}

/// Extend forward from `r_start`/`a_start`, ungapped, with +1 for a match and
/// `mismatch` for a difference. Stop when the running score has fallen `xdrop`
/// below its best value, and keep only up to that best (6.6 step D).
pub fn xdrop_forward(
    ref_s: &[u8],
    alt: &[u8],
    r_start: usize,
    a_start: usize,
    scores: &AlignScores,
    xdrop: i32,
) -> XdropResult {
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_len = 0usize;
    let mut best_mismatches = 0usize;
    let mut mismatches: Vec<usize> = Vec::new();
    let mut hit_end = true;
    let mut i = 0usize;
    loop {
        let rp = r_start + i;
        let ap = a_start + i;
        if rp >= ref_s.len() || ap >= alt.len() {
            break;
        }
        let eq = ref_s[rp] == alt[ap];
        score += if eq { scores.r#match } else { scores.mismatch };
        if !eq {
            mismatches.push(i);
        }
        i += 1;
        if score > best {
            best = score;
            best_len = i;
            best_mismatches = mismatches.len();
        }
        if best - score >= xdrop {
            hit_end = false;
            break;
        }
    }
    mismatches.truncate(best_mismatches);
    XdropResult {
        len: best_len,
        mismatches,
        hit_end,
    }
}

/// The same, extending backwards from `r_start`/`a_start` (both exclusive).
pub fn xdrop_backward(
    ref_s: &[u8],
    alt: &[u8],
    r_start: usize,
    a_start: usize,
    scores: &AlignScores,
    xdrop: i32,
) -> XdropResult {
    let mut score = 0i32;
    let mut best = 0i32;
    let mut best_len = 0usize;
    let mut best_mismatches = 0usize;
    let mut mismatches: Vec<usize> = Vec::new();
    let mut hit_end = true;
    let mut i = 0usize;
    loop {
        if r_start == 0 || a_start == 0 || i >= r_start || i >= a_start {
            break;
        }
        let rp = r_start - i - 1;
        let ap = a_start - i - 1;
        let eq = ref_s[rp] == alt[ap];
        score += if eq { scores.r#match } else { scores.mismatch };
        if !eq {
            // Offsets are reported as distance back from the start position.
            mismatches.push(i);
        }
        i += 1;
        if score > best {
            best = score;
            best_len = i;
            best_mismatches = mismatches.len();
        }
        if best - score >= xdrop {
            hit_end = false;
            break;
        }
    }
    mismatches.truncate(best_mismatches);
    XdropResult {
        len: best_len,
        mismatches,
        hit_end,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(x: &str) -> Vec<u8> {
        x.as_bytes().to_vec()
    }

    #[test]
    fn identical_sequences_align_with_no_gaps() {
        let a = s("ACGTACGTACGT");
        let al = global_align(&a, &a, &AlignScores::default()).unwrap();
        assert_eq!(al.matches, a.len());
        assert_eq!(al.columns, a.len());
        assert_eq!(al.identity(), 1.0);
        assert!(al.ops.iter().all(|o| *o == Op::Match));
    }

    #[test]
    fn a_single_substitution_is_one_mismatch() {
        let r = s("ACGTACGTACGT");
        let a = s("ACGAACGTACGT");
        let al = global_align(&r, &a, &AlignScores::default()).unwrap();
        assert_eq!(al.ops.iter().filter(|o| **o == Op::Mismatch).count(), 1);
        assert_eq!(al.matches, 11);
    }

    #[test]
    fn a_single_insertion_is_one_gap() {
        let r = s("ACGTACGTACGT");
        let a = s("ACGTTACGTACGT");
        let al = global_align(&r, &a, &AlignScores::default()).unwrap();
        let ins = al.ops.iter().filter(|o| **o == Op::Insertion).count();
        let del = al.ops.iter().filter(|o| **o == Op::Deletion).count();
        assert_eq!(ins, 1);
        assert_eq!(del, 0);
        assert_eq!(al.matches, 12);
    }

    #[test]
    fn a_single_deletion_is_one_gap() {
        let r = s("ACGTTACGTACGT");
        let a = s("ACGTACGTACGT");
        let al = global_align(&r, &a, &AlignScores::default()).unwrap();
        assert_eq!(al.ops.iter().filter(|o| **o == Op::Deletion).count(), 1);
        assert_eq!(al.matches, 12);
    }

    #[test]
    fn empty_sides_are_all_gaps() {
        let r = s("ACGT");
        let al = global_align(&r, &[], &AlignScores::default()).unwrap();
        assert_eq!(al.ops, vec![Op::Deletion; 4]);
        let al2 = global_align(&[], &r, &AlignScores::default()).unwrap();
        assert_eq!(al2.ops, vec![Op::Insertion; 4]);
        let al3 = global_align(&[], &[], &AlignScores::default()).unwrap();
        assert_eq!(al3.columns, 0);
        assert_eq!(al3.identity(), 1.0);
    }

    #[test]
    fn alignment_consumes_both_sequences_exactly() {
        // Whatever the band does, the ops must account for every base of both
        // inputs: matches+mismatches+deletions == len(ref),
        // matches+mismatches+insertions == len(alt).
        let r = s("ACGTACGTTTGGGGCAAACCGGTT");
        let a = s("ACGAACGGTTGGGCAAACCAGGTT");
        let al = global_align(&r, &a, &AlignScores::default()).unwrap();
        let (mut m, mut mm, mut i, mut d) = (0, 0, 0, 0);
        for o in &al.ops {
            match o {
                Op::Match => m += 1,
                Op::Mismatch => mm += 1,
                Op::Insertion => i += 1,
                Op::Deletion => d += 1,
            }
        }
        assert_eq!(m + mm + d, r.len());
        assert_eq!(m + mm + i, a.len());
    }

    #[test]
    fn long_pieces_stay_inside_the_band_and_stay_correct() {
        // 4 kb with a 30 bp insertion in the middle: the band must carry the
        // path and the op counts must still balance.
        let mut r: Vec<u8> = (0..4000).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let mut a = r.clone();
        let ins: Vec<u8> = (0..30).map(|_| b'T').collect();
        a.splice(2000..2000, ins.iter().copied());
        let al = global_align(&r, &a, &AlignScores::default()).expect("band wide enough");
        let ins_count = al.ops.iter().filter(|o| **o == Op::Insertion).count();
        assert!(ins_count >= 30, "got {ins_count}");
        let (mut m, mut mm, mut i, mut d) = (0, 0, 0, 0);
        for o in &al.ops {
            match o {
                Op::Match => m += 1,
                Op::Mismatch => mm += 1,
                Op::Insertion => i += 1,
                Op::Deletion => d += 1,
            }
        }
        assert_eq!(m + mm + d, r.len());
        assert_eq!(m + mm + i, a.len());
        r.clear();
    }

    #[test]
    fn xdrop_keeps_the_best_prefix_and_reports_mismatches() {
        // One mismatch early, then an unbroken run of matches to the end.
        let r = s("AAAAAAAAAACCCCCCCCCC");
        let a = s("AAAATAAAAACCCCCCCCCC");
        let res = xdrop_forward(&r, &a, 0, 0, &AlignScores::default(), 15);
        assert_eq!(res.len, 20, "the whole run is kept");
        assert_eq!(res.mismatches, vec![4]);
        assert!(res.hit_end, "ran out of sequence rather than dropping");
    }

    #[test]
    fn xdrop_walk_handles_opposite_directions() {
        // The member runs backwards while the reference runs forwards, which is
        // what 6.6 step D meets when the member is reverse-complemented.
        // r forward from index 6 is ten A's; a backward from index 9 is the
        // same ten A's read the other way.
        let r = s("GGGGGGAAAAAAAAAA");
        let a = s("AAAAAAAAAAGGGGGG");
        let res = xdrop_walk(&r, &a, 6, 9, 1, -1, false, &AlignScores::default(), 15);
        assert_eq!(res.len, 10);
        assert!(res.mismatches.is_empty());
        assert!(res.hit_end, "both sides ran out together");
    }

    #[test]
    fn xdrop_walk_complements_the_alt_strand_when_asked() {
        // The antiparallel case of 6.6 step D: a member contig that runs
        // opposite to the pivot. Compared directly every position mismatches;
        // compared as complements every position matches.
        let r = s("AAAAAAAAAAGGGGGGGGGG");
        let a = s("CCCCCCCCCCTTTTTTTTTT"); // reverse complement of r
        let scores = AlignScores::default();
        // Walk r forward from 0 and a backward from its last base.
        let plain = xdrop_walk(&r, &a, 0, 19, 1, -1, false, &scores, 15);
        assert_eq!(plain.len, 0, "without complementing nothing matches");
        let cmp = xdrop_walk(&r, &a, 0, 19, 1, -1, true, &scores, 15);
        assert_eq!(cmp.len, 20, "with complementing the whole run matches");
        assert!(cmp.mismatches.is_empty());
    }

    #[test]
    fn xdrop_walk_complement_still_reports_real_mismatches() {
        // One base of the member is not the complement of the pivot: that must
        // show up as exactly one mismatch, not as a collapse of the walk.
        let r = s("AAAAAAAAAAGGGGGGGGGG");
        let mut a = s("CCCCCCCCCCTTTTTTTTTT");
        a[4] = b'A'; // should be C, to complement the G it faces
        let scores = AlignScores::default();
        let res = xdrop_walk(&r, &a, 0, 19, 1, -1, true, &scores, 15);
        assert_eq!(res.len, 20);
        // a[4] is reached at step 19 - 4 = 15.
        assert_eq!(res.mismatches, vec![15]);
    }

    #[test]
    fn xdrop_walk_forward_matches_xdrop_forward() {
        let r = s("AAAAAAAAAACCCCCCCCCC");
        let a = s("AAAATAAAAACCCCCCCCCC");
        let w = xdrop_walk(&r, &a, 0, 0, 1, 1, false, &AlignScores::default(), 15);
        let f = xdrop_forward(&r, &a, 0, 0, &AlignScores::default(), 15);
        assert_eq!(w.len, f.len);
        assert_eq!(w.mismatches, f.mismatches);
        assert_eq!(w.hit_end, f.hit_end);
    }

    #[test]
    fn walked_range_covers_exactly_the_steps_taken() {
        assert_eq!(walked_range(100, 1, 10), (100, 110));
        assert_eq!(walked_range(100, -1, 10), (91, 101));
        assert_eq!(walked_range(100, 1, 0), (100, 100));
    }

    #[test]
    fn xdrop_stops_on_divergence() {
        let r = s("AAAAAAAAAAGGGGGGGGGG");
        let a = s("AAAAAAAAAATTTTTTTTTT");
        let res = xdrop_forward(&r, &a, 0, 0, &AlignScores::default(), 15);
        assert_eq!(res.len, 10, "kept up to the best score");
        assert!(!res.hit_end);
    }

    #[test]
    fn xdrop_backward_mirrors_forward() {
        let r = s("CCCCCCCCCCAAAAAAAAAA");
        let a = s("CCCCCCCCCAAAAAAAAAAA");
        let back = xdrop_backward(&r, &a, 20, 20, &AlignScores::default(), 15);
        assert_eq!(back.len, 20);
        assert_eq!(back.mismatches, vec![10]);
    }

    #[test]
    fn xdrop_reports_running_out_of_sequence() {
        let r = s("AAAA");
        let a = s("AAAA");
        let res = xdrop_forward(&r, &a, 0, 0, &AlignScores::default(), 15);
        assert_eq!(res.len, 4);
        assert!(res.hit_end);
    }
}
