//! Nested content-defined blocks and their strand-independent hashes (6.3).
//!
//! A block is the stretch of a contig from the start of one anchor to the end
//! of the next anchor of the same level, so neighbouring blocks overlap by
//! exactly `k` bases. That overlap is what makes a block self-delimiting: the
//! cut points come from the sequence itself, so the same region of two genomes
//! produces the same block boundaries even when one of them has an indel
//! upstream.

use serde::{Deserialize, Serialize};
use xxhash_rust::xxh3::xxh3_64;

use crate::seq::revcomp;

/// One block of one contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Block {
    /// 0-based, half-open, in contig coordinates.
    pub start: u64,
    pub end: u64,
    /// `min(xxh3(bytes), xxh3(revcomp(bytes)))`.
    pub hash: u64,
    /// True when the reverse-complement encoding was the smaller one.
    pub flipped: bool,
}

impl Block {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }

    /// Position of the anchor that opens this block, if any. Head blocks that
    /// start at 0 with the first anchor later on still report that anchor.
    pub fn anchor_start(&self) -> u64 {
        self.start
    }
}

/// The block boundaries of one contig at one level, given the ascending start
/// positions of that level's anchors.
///
/// Layout, exactly as 6.3:
/// - interior block `i`: `[p_i, p_{i+1} + k)`
/// - head block: `[0, p_1 + k)`, omitted when `p_1 == 0`
/// - tail block: `[p_n, L)`, omitted when `p_n + k == L`
/// - no anchors at all: the whole contig is one block
pub fn block_ranges(anchor_positions: &[u64], contig_len: u64, k: u64) -> Vec<(u64, u64)> {
    if contig_len == 0 {
        return Vec::new();
    }
    let n = anchor_positions.len();
    if n == 0 {
        return vec![(0, contig_len)];
    }
    let mut out = Vec::with_capacity(n + 1);
    let first = anchor_positions[0];
    let last = anchor_positions[n - 1];

    // Head.
    if first != 0 {
        out.push((0, first + k));
    }
    // Interior.
    for w in anchor_positions.windows(2) {
        let end = (w[1] + k).min(contig_len);
        out.push((w[0], end));
    }
    // Tail.
    if last + k != contig_len {
        // The tail overlaps the last interior block by k bases.
        if last < contig_len {
            out.push((last, contig_len));
        }
    }
    out
}

/// Strand-independent hash of a byte slice (6.3).
pub fn hash_block(seq: &[u8]) -> (u64, bool) {
    let hf = xxh3_64(seq);
    let rc = revcomp(seq);
    let hr = xxh3_64(&rc);
    if hr < hf {
        (hr, true)
    } else {
        (hf, false)
    }
}

/// All blocks of one contig at one level, hashed.
pub fn contig_blocks(seq: &[u8], anchor_positions: &[u64], k: usize) -> Vec<Block> {
    let ranges = block_ranges(anchor_positions, seq.len() as u64, k as u64);
    let mut out = Vec::with_capacity(ranges.len());
    for (s, e) in ranges {
        let slice = &seq[s as usize..e as usize];
        let (hash, flipped) = hash_block(slice);
        out.push(Block {
            start: s,
            end: e,
            hash,
            flipped,
        });
    }
    out
}

/// All blocks of all contigs of a genome at one level.
///
/// `anchors_by_contig` comes from [`crate::genome::anchors_by_contig`].
pub fn genome_blocks(
    contigs: &[crate::fasta::Contig],
    anchors_by_contig: &[Vec<u64>],
    k: usize,
) -> Vec<Vec<Block>> {
    contigs
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let pos = anchors_by_contig
                .get(i)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            contig_blocks(&c.seq, pos, k)
        })
        .collect()
}

/// Map a member position onto the pivot through a matched block pair (6.6 A).
///
/// `rel` is `flipped_member XOR flipped_pivot`. Positions are 0-based.
#[inline]
pub fn map_through_block(member_x: u64, gs: u64, ps: u64, pe: u64, rel: bool) -> u64 {
    if rel {
        pe - 1 - (member_x - gs)
    } else {
        ps + (member_x - gs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_anchors_means_one_block_for_the_whole_contig() {
        assert_eq!(block_ranges(&[], 100, 5), vec![(0, 100)]);
        assert_eq!(block_ranges(&[], 0, 5), Vec::<(u64, u64)>::new());
    }

    #[test]
    fn blocks_tiling_matches_section_6_3() {
        // L = 50, k = 5, anchors at 10, 25, 40.
        let r = block_ranges(&[10, 25, 40], 50, 5);
        assert_eq!(
            r,
            vec![
                (0, 15),  // head [0, p1+k)
                (10, 30), // interior [p1, p2+k)
                (25, 45), // interior [p2, p3+k)
                (40, 50), // tail [p3, L) because 40+5 != 50
            ]
        );
    }

    #[test]
    fn head_omitted_when_first_anchor_is_at_zero() {
        let r = block_ranges(&[0, 20], 30, 5);
        assert_eq!(r, vec![(0, 25), (20, 30)]);
    }

    #[test]
    fn tail_omitted_when_last_anchor_ends_the_contig() {
        let r = block_ranges(&[0, 20], 25, 5);
        assert_eq!(r, vec![(0, 25)], "20 + 5 == L so there is no tail");
    }

    #[test]
    fn neighbours_overlap_by_exactly_k() {
        let k = 7u64;
        let r = block_ranges(&[3, 40, 90, 150], 200, k);
        for w in r.windows(2) {
            // Head -> first interior also overlaps by k.
            assert!(w[0].1 > w[1].0, "{:?}", w);
        }
        // interior blocks overlap by exactly k
        assert_eq!(r[1].1 - r[2].0, k);
        assert_eq!(r[2].1 - r[3].0, k);
    }

    #[test]
    fn block_hash_ignores_strand() {
        let s = b"ACGTACGTGGGTCAAACGTACGT";
        let rc = revcomp(s);
        let (h1, f1) = hash_block(s);
        let (h2, f2) = hash_block(&rc);
        assert_eq!(h1, h2, "the hash must not depend on orientation");
        assert_ne!(f1, f2, "the flipped flag must");
    }

    #[test]
    fn different_sequence_gives_a_different_hash() {
        let a = b"ACGTACGTGGGTCAAACGTACGT";
        let mut b = a.to_vec();
        b[10] = b'A';
        assert_ne!(hash_block(a).0, hash_block(&b).0);
    }

    #[test]
    fn mapping_through_a_block_is_reversible() {
        // rel = false: forward mapping.
        assert_eq!(map_through_block(15, 10, 100, 200, false), 105);
        // rel = true: the member block is the reverse complement.
        // member [10,20) maps to pivot [100,110) reversed.
        assert_eq!(map_through_block(10, 10, 100, 110, true), 109);
        assert_eq!(map_through_block(19, 10, 100, 110, true), 100);
        // round trip
        for x in 10..20u64 {
            let p = map_through_block(x, 10, 100, 110, true);
            let back = map_through_block(p, 100, 10, 20, true);
            assert_eq!(back, x);
        }
    }
}
