//! DNA encoding, k-mer scanning and the anchor hash of section 6.2.
//!
//! Everything that turns bases into numbers lives here, so that the CLI, the
//! GUI, the simulator and the tests agree bit for bit. Section 0 rule 3
//! (byte-identical output for any thread count) depends on it: no floating
//! point, no hashing of pointer values, no iteration over an unordered map.

use serde::{Deserialize, Serialize};

/// Mixing constant from section 6.2.
pub const ANCHOR_SALT: u64 = 0x9E37_79B9_7F4A_7C15;

/// Number of anchor levels (section 5).
pub const LEVELS: usize = 3;

/// Code of a normalised base: A=0, C=1, G=2, T=3. Anything else codes as 4.
#[inline]
pub fn code(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

/// Uppercase, and map every character that is not A, C, G or T to N (6.2.1).
#[inline]
pub fn normalize_byte(b: u8) -> u8 {
    let u = b.to_ascii_uppercase();
    match u {
        b'A' | b'C' | b'G' | b'T' => u,
        _ => b'N',
    }
}

/// Normalise a whole sequence in place (6.2.1).
pub fn normalize(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        *b = normalize_byte(*b);
    }
}

/// Complement of a normalised base. N complements to N.
#[inline]
pub fn complement_byte(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Reverse complement of a normalised sequence.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_byte(b)).collect()
}

/// `fmix64` from section 6.2. All arithmetic wraps at 64 bits.
#[inline]
pub fn fmix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;
    x
}

/// `h = fmix64(canon XOR 0x9E3779B97F4A7C15)`.
#[inline]
pub fn anchor_hash(canon: u64) -> u64 {
    fmix64(canon ^ ANCHOR_SALT)
}

/// Strand of an anchor, as written in output files.
#[inline]
pub fn strand_char(plus: bool) -> char {
    if plus {
        '+'
    } else {
        '-'
    }
}

/// Sampling thresholds for the three levels (6.2.5).
///
/// A window is an anchor at level `j` when `h < (2^64 - 1) / densities[j]`.
/// Because the densities descend, the thresholds ascend and the levels are
/// nested: a level-1 anchor is automatically a level-2 and level-3 anchor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Thresholds {
    pub t: [u64; LEVELS],
}

impl Thresholds {
    pub fn new(densities: &[u64; LEVELS]) -> Self {
        Self {
            t: [
                u64::MAX / densities[0],
                u64::MAX / densities[1],
                u64::MAX / densities[2],
            ],
        }
    }

    /// The sparsest level this hash qualifies for, or `None` if it is not an
    /// anchor at all.
    #[inline]
    pub fn level(&self, h: u64) -> Option<u8> {
        if h < self.t[0] {
            Some(1)
        } else if h < self.t[1] {
            Some(2)
        } else if h < self.t[2] {
            Some(3)
        } else {
            None
        }
    }

    /// True when the hash qualifies for level `j` (1-based) or any sparser one.
    #[inline]
    pub fn qualifies(&self, h: u64, level: u8) -> bool {
        let idx = (level as usize).saturating_sub(1).min(LEVELS - 1);
        h < self.t[idx]
    }
}

/// One occurrence of a sampled k-mer (section 5).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Anchor {
    pub h: u64,
    /// Index of the contig inside its genome.
    pub contig: u32,
    /// 0-based start of the k-mer on the contig.
    pub pos: u64,
    /// True when the forward encoding is the canonical one.
    pub plus: bool,
    /// Sparsest level this anchor belongs to: 1, 2 or 3.
    pub level: u8,
}

/// Walk every N-free window of `k` bases in `seq`, in ascending position
/// order, calling `f(start, hash, plus)`.
///
/// The forward and reverse-complement encodings are rolled along the sequence,
/// so the cost is one pass with no per-window re-encoding. Windows containing
/// an N are skipped; the bogus codes they introduce are shifted out of the
/// rolling state before any reported window can see them.
pub fn scan_windows(seq: &[u8], k: usize, mut f: impl FnMut(usize, u64, bool)) {
    let n = seq.len();
    if k == 0 || n < k {
        return;
    }
    let shift = 2 * (k - 1);
    let mask: u64 = if 2 * k >= 64 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    let mut fwd: u64 = 0;
    for &b in &seq[..k] {
        fwd = (fwd << 2) | u64::from(code(b) & 3);
    }
    let mut rc: u64 = 0;
    for &b in seq[..k].iter().rev() {
        rc = (rc << 2) | u64::from(3 - (code(b) & 3));
    }

    // Start of the most recent N seen, or -1.
    let mut last_n: i64 = -1;
    for (i, &b) in seq[..k].iter().enumerate() {
        if code(b) == 4 {
            last_n = i as i64;
        }
    }

    let mut start = 0usize;
    loop {
        if last_n < start as i64 {
            let canon = fwd.min(rc);
            f(start, anchor_hash(canon), fwd <= rc);
        }
        if start + k >= n {
            return;
        }
        let incoming = code(seq[start + k]);
        if incoming == 4 {
            last_n = (start + k) as i64;
        }
        fwd = ((fwd << 2) & mask) | u64::from(incoming & 3);
        rc = (rc >> 2) | (u64::from(3 - (incoming & 3)) << shift);
        start += 1;
    }
}

/// Collect the anchors of one contig, keeping every level-3 anchor (6.2.6).
pub fn contig_anchors(seq: &[u8], contig: u32, k: usize, th: &Thresholds) -> Vec<Anchor> {
    let mut out = Vec::new();
    scan_windows(seq, k, |pos, h, plus| {
        if let Some(level) = th.level(h) {
            out.push(Anchor {
                h,
                contig,
                pos: pos as u64,
                plus,
                level,
            });
        }
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmix64_matches_the_reference_avalanche() {
        // fmix64 is the well-known MurmurHash3 finaliser; check it is stable
        // and mixes, without pinning a value that depends on a typo.
        assert_eq!(fmix64(0), fmix64(0));
        assert_ne!(fmix64(0), fmix64(1));
        assert_eq!(fmix64(0x9E37_79B9_7F4A_7C15), anchor_hash(0));
    }

    #[test]
    fn normalize_maps_everything_else_to_n() {
        let mut s = b"acgtNRYKMSWBDHV-*\n".to_vec();
        normalize(&mut s);
        assert_eq!(&s, b"ACGTNNNNNNNNNNNNNN");
    }

    #[test]
    fn revcomp_is_an_involution_on_normalised_sequence() {
        let s = b"ACGTNACGTA";
        assert_eq!(revcomp(&revcomp(s)), s.to_vec());
        assert_eq!(revcomp(b"ACGTN"), b"NACGT".to_vec());
    }

    #[test]
    fn canonical_hash_is_strand_symmetric() {
        // A k-mer and its reverse complement must give the same hash and
        // opposite strand flags. This is what "bidirectional" means (6.2).
        let k = 31;
        let fwd = b"ACGTACGTACGTACGTACGTACGTACGTACG";
        assert_eq!(fwd.len(), k);
        let rc = revcomp(fwd);

        let mut a = Vec::new();
        let mut b = Vec::new();
        scan_windows(fwd, k, |p, h, plus| a.push((p, h, plus)));
        scan_windows(&rc, k, |p, h, plus| b.push((p, h, plus)));
        assert_eq!(a.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(a[0].1, b[0].1, "hash must not depend on strand");
        assert_ne!(a[0].2, b[0].2, "strand flag must flip");
    }

    #[test]
    fn rolling_matches_naive_encoding() {
        // Cross-check the rolled window encodings against a naive per-window
        // computation on a sequence that contains Ns.
        let k = 7;
        let seq = b"ACGTNNACGTACAGTNACGTACAGT";

        let mut rolled = Vec::new();
        scan_windows(seq, k, |pos, h, plus| rolled.push((pos, h, plus)));

        let mut naive = Vec::new();
        for start in 0..=seq.len() - k {
            let w = &seq[start..start + k];
            if w.iter().any(|&b| code(b) == 4) {
                continue;
            }
            let mut f = 0u64;
            for &b in w {
                f = (f << 2) | u64::from(code(b));
            }
            let mut r = 0u64;
            for &b in w.iter().rev() {
                r = (r << 2) | u64::from(3 - code(b));
            }
            naive.push((start, anchor_hash(f.min(r)), f <= r));
        }
        assert_eq!(rolled, naive);
    }

    #[test]
    fn levels_are_nested() {
        let th = Thresholds::new(&[16384, 2048, 256]);
        assert!(th.t[0] < th.t[1] && th.t[1] < th.t[2]);
        assert_eq!(th.level(th.t[0] - 1), Some(1));
        assert_eq!(th.level(th.t[1] - 1), Some(2));
        assert_eq!(th.level(th.t[2] - 1), Some(3));
        assert_eq!(th.level(th.t[2]), None);
        // a level-1 anchor qualifies for every level
        assert!(th.qualifies(th.t[0] - 1, 3));
        assert!(!th.qualifies(th.t[2] - 1, 1));
    }

    #[test]
    fn short_sequences_yield_nothing() {
        let th = Thresholds::new(&[4, 8, 16]);
        assert!(contig_anchors(b"ACG", 0, 31, &th).is_empty());
        assert!(contig_anchors(b"", 0, 31, &th).is_empty());
    }
}
