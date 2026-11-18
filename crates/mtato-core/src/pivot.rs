//! The pivot index (6.5).
//!
//! For the pivot of a cluster being refined, hold in memory: its sequence, its
//! level-3 anchors sorted by position, an anchor map of canonical k-mer hash to
//! position for anchors that occur exactly once, and one block dictionary per
//! level. Clusters with a single genome need no index.

use std::collections::HashMap;

use crate::blocks::{block_ranges, hash_block, Block};
use crate::fasta::{Contig, GenomeData};
use crate::params::Parameters;
use crate::seq::{scan_windows, Anchor, Thresholds, LEVELS};

/// Where a block sits in the pivot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockLoc {
    pub contig: u32,
    pub start: u64,
    pub end: u64,
    /// True when the block's reverse-complement encoding was the smaller one.
    pub flipped: bool,
}

/// Where an anchor sits in the pivot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnchorLoc {
    pub contig: u32,
    pub pos: u64,
    pub plus: bool,
}

/// Everything needed to compare members of one cluster against its pivot.
///
/// Built once per cluster and shared by all members, which is what makes the
/// per-member cost proportional to what differs rather than to genome size.
#[derive(Debug)]
pub struct PivotIndex {
    pub name: String,
    pub contigs: Vec<Contig>,
    pub contig_lens: Vec<u64>,
    pub total_len: u64,
    pub k: usize,
    pub th: Thresholds,
    /// All level-3 anchors, ordered by contig then position.
    pub anchors: Vec<Anchor>,
    /// Anchor start positions per level, per contig. Index `[level-1][contig]`.
    pub anchor_pos: Vec<Vec<Vec<u64>>>,
    /// Blocks per level, per contig. Index `[level-1][contig]`.
    pub blocks: Vec<Vec<Vec<Block>>>,
    /// Block hash -> occurrences, per level. Occurrences are stored in contig
    /// then position order, so "take the first" is deterministic.
    pub block_dict: Vec<HashMap<u64, Vec<BlockLoc>>>,
    /// Canonical anchor hash -> position, for anchors occurring exactly once
    /// among the pivot's anchors. Repeats are absent on purpose: mapping a
    /// repeated k-mer to one arbitrary copy would invent a position (section
    /// 15).
    pub anchor_map: HashMap<u64, AnchorLoc>,
    /// Number of distinct hashes that were dropped from `anchor_map` because
    /// they occur more than once. Reported as a warning: it measures how much
    /// of the pivot is repetitive.
    pub repeated_anchor_hashes: u64,
}

impl PivotIndex {
    /// Build the index for a pivot genome.
    pub fn new(name: &str, data: &GenomeData, params: &Parameters) -> Self {
        let k = params.k;
        let th = Thresholds::new(&params.densities);
        let n_contigs = data.contigs.len();

        // One scan produces every level, because the level thresholds are
        // nested (6.2.5).
        let mut anchors: Vec<Anchor> = Vec::new();
        let mut anchor_pos = vec![vec![Vec::new(); n_contigs]; LEVELS];
        for (ci, contig) in data.contigs.iter().enumerate() {
            let cid = ci as u32;
            scan_windows(&contig.seq, k, |pos, h, plus| {
                if let Some(level) = th.level(h) {
                    anchors.push(Anchor {
                        h,
                        contig: cid,
                        pos: pos as u64,
                        plus,
                        level,
                    });
                    // Nested thresholds: a level-`level` anchor is an anchor at
                    // every finer level too (6.2.5).
                    for l in level as usize..=LEVELS {
                        anchor_pos[l - 1][ci].push(pos as u64);
                    }
                }
            });
        }

        // Block dictionaries, one per level.
        let mut blocks = vec![Vec::new(); LEVELS];
        let mut block_dict: Vec<HashMap<u64, Vec<BlockLoc>>> =
            (0..LEVELS).map(|_| HashMap::new()).collect();
        for level in 1..=LEVELS {
            let mut per_contig = Vec::with_capacity(n_contigs);
            for (ci, contig) in data.contigs.iter().enumerate() {
                let pos = &anchor_pos[level - 1][ci];
                let ranges = block_ranges(pos, contig.len(), k as u64);
                let mut v = Vec::with_capacity(ranges.len());
                for (s, e) in ranges {
                    let (hash, flipped) = hash_block(&contig.seq[s as usize..e as usize]);
                    let b = Block {
                        start: s,
                        end: e,
                        hash,
                        flipped,
                    };
                    block_dict[level - 1]
                        .entry(hash)
                        .or_default()
                        .push(BlockLoc {
                            contig: ci as u32,
                            start: s,
                            end: e,
                            flipped,
                        });
                    v.push(b);
                }
                per_contig.push(v);
            }
            blocks[level - 1] = per_contig;
        }

        // Anchor map: keep only hashes that occur exactly once.
        let mut counts: HashMap<u64, u32> = HashMap::with_capacity(anchors.len());
        for a in &anchors {
            *counts.entry(a.h).or_insert(0) += 1;
        }
        let mut anchor_map: HashMap<u64, AnchorLoc> = HashMap::with_capacity(anchors.len());
        for a in &anchors {
            if counts.get(&a.h) == Some(&1) {
                anchor_map.insert(
                    a.h,
                    AnchorLoc {
                        contig: a.contig,
                        pos: a.pos,
                        plus: a.plus,
                    },
                );
            }
        }
        let distinct_repeated = counts.values().filter(|&&c| c > 1).count() as u64;

        let contig_lens: Vec<u64> = data.contigs.iter().map(|c| c.len()).collect();
        let total_len = contig_lens.iter().sum();

        Self {
            name: name.to_string(),
            contigs: data.contigs.clone(),
            contig_lens,
            total_len,
            k,
            th,
            anchors,
            anchor_pos,
            blocks,
            block_dict,
            anchor_map,
            repeated_anchor_hashes: distinct_repeated,
        }
    }

    /// Concatenated pivot sequence of one contig.
    pub fn contig_seq(&self, c: u32) -> &[u8] {
        &self.contigs[c as usize].seq
    }

    pub fn contig_name(&self, c: u32) -> &str {
        &self.contigs[c as usize].name
    }

    /// Look a block hash up at one level (1-based).
    pub fn lookup_block(&self, level: usize, hash: u64) -> Option<&[BlockLoc]> {
        self.block_dict
            .get(level - 1)?
            .get(&hash)
            .map(|v| v.as_slice())
    }

    /// Blocks of one contig at one level, in position order.
    pub fn contig_blocks(&self, level: usize, contig: usize) -> &[Block] {
        self.blocks[level - 1]
            .get(contig)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Anchor start positions of one contig at one level.
    pub fn contig_anchor_pos(&self, level: usize, contig: usize) -> &[u64] {
        self.anchor_pos[level - 1]
            .get(contig)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Level-3 blocks of one contig whose start lies in `[s, e.saturating_sub(k)]`.
    ///
    /// This is how the zoom descends (6.10): the finer blocks of a level `j`
    /// block `[s, e)` are exactly those whose anchor start falls in that range,
    /// because a level-`j` block spanning anchors `p_i .. p_{i+1}` is the union
    /// of the level-`j+1` blocks starting at every anchor from `p_i` to
    /// `p_{i+1}` inclusive.
    pub fn inner_blocks(&self, level: usize, contig: usize, s: u64, e: u64) -> Vec<Block> {
        let lo = s;
        let hi = e.saturating_sub(self.k as u64);
        self.contig_blocks(level, contig)
            .iter()
            .filter(|b| b.start >= lo && b.start <= hi)
            .copied()
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::{anchor_hash, code};

    fn params(k: usize, densities: [u64; LEVELS]) -> Parameters {
        Parameters {
            k,
            densities,
            ..Default::default()
        }
    }

    fn genome(seqs: &[&[u8]]) -> GenomeData {
        GenomeData {
            name: "pivot".into(),
            contigs: seqs
                .iter()
                .enumerate()
                .map(|(i, s)| Contig {
                    name: format!("c{i}"),
                    seq: s.to_vec(),
                })
                .collect(),
        }
    }

    fn randomish(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            out.push(b"ACGT"[(x >> 33) as usize & 3]);
        }
        out
    }

    #[test]
    fn index_covers_the_whole_contig_at_every_level() {
        let seq = randomish(20_000, 7);
        let g = genome(&[&seq]);
        let p = params(21, [256, 64, 16]);
        let idx = PivotIndex::new("pivot", &g, &p);

        for level in 1..=LEVELS {
            let blocks = idx.contig_blocks(level, 0);
            assert!(!blocks.is_empty(), "level {level} has no blocks");
            assert_eq!(blocks[0].start, 0, "level {level} must start at 0");
            assert_eq!(
                blocks.last().unwrap().end,
                seq.len() as u64,
                "level {level} must reach the end"
            );
            // Consecutive blocks overlap by exactly k and never go backwards.
            for w in blocks.windows(2) {
                assert!(w[1].start > w[0].start);
                assert!(w[1].start < w[0].end);
            }
        }
    }

    #[test]
    fn coarser_levels_have_fewer_blocks_and_nesting_holds() {
        let seq = randomish(40_000, 11);
        let g = genome(&[&seq]);
        let p = params(21, [256, 64, 16]);
        let idx = PivotIndex::new("pivot", &g, &p);
        let n1 = idx.contig_blocks(1, 0).len();
        let n2 = idx.contig_blocks(2, 0).len();
        let n3 = idx.contig_blocks(3, 0).len();
        assert!(n1 < n2 && n2 < n3, "{n1} {n2} {n3}");
        assert!(
            idx.contig_anchor_pos(1, 0).len() <= idx.contig_anchor_pos(2, 0).len()
                && idx.contig_anchor_pos(2, 0).len() <= idx.contig_anchor_pos(3, 0).len()
        );
    }

    #[test]
    fn inner_blocks_cover_their_parent() {
        // The zoom of 6.10 is only correct if the finer blocks of a coarse block
        // cover at least the same span. The last inner block legitimately runs
        // past the parent's end, because blocks overlap by k and the parent ends
        // at an anchor that is also the start of an inner block.
        let seq = randomish(30_000, 3);
        let g = genome(&[&seq]);
        let p = params(21, [256, 64, 16]);
        let idx = PivotIndex::new("pivot", &g, &p);
        for b1 in idx.contig_blocks(1, 0) {
            let inner = idx.inner_blocks(3, 0, b1.start, b1.end);
            assert!(!inner.is_empty(), "empty interior for {:?}", b1);
            assert_eq!(inner[0].start, b1.start);
            assert!(
                inner.last().unwrap().end >= b1.end,
                "inner blocks must cover the parent: {:?} vs {:?}",
                inner.last().unwrap(),
                b1
            );
            // Contiguous, overlapping by exactly k.
            for w in inner.windows(2) {
                assert_eq!(w[1].start + p.k as u64, w[0].end);
            }
        }
    }

    #[test]
    fn every_leaf_block_is_reachable_by_descent() {
        // The set of leaf blocks the hierarchical descent can reach must be
        // exactly the set flat mode looks at, otherwise the two modes could
        // match different bases and disagree (property 6, section 10.2).
        let seq = randomish(30_000, 5);
        let g = genome(&[&seq]);
        let p = params(21, [256, 64, 16]);
        let idx = PivotIndex::new("pivot", &g, &p);
        let k = p.k as u64;

        let mut reached: Vec<(u64, u64)> = Vec::new();
        for b1 in idx.contig_blocks(1, 0) {
            for b2 in idx.inner_blocks(2, 0, b1.start, b1.end) {
                for b3 in idx.inner_blocks(3, 0, b2.start, b2.end) {
                    reached.push((b3.start, b3.end));
                }
            }
        }
        reached.sort_unstable();
        reached.dedup();

        let flat: Vec<(u64, u64)> = idx
            .contig_blocks(3, 0)
            .iter()
            .map(|b| (b.start, b.end))
            .collect();
        assert_eq!(reached, flat, "descent must reach every leaf block");
        let _ = k;
    }

    #[test]
    fn repeated_kmers_are_excluded_from_the_anchor_map() {
        // A tandem repeat: the same 21-mer occurs twice, so it must not appear
        // in the anchor map at all.
        let unit = randomish(60, 5);
        let mut seq = randomish(2000, 9);
        seq.extend_from_slice(&unit);
        seq.extend_from_slice(&randomish(500, 13));
        seq.extend_from_slice(&unit);
        seq.extend_from_slice(&randomish(2000, 17));
        let g = genome(&[&seq]);
        let p = params(21, [16, 8, 4]);
        let idx = PivotIndex::new("pivot", &g, &p);

        // Count duplicates directly and compare with what the index reports.
        let mut counts: HashMap<u64, u32> = HashMap::new();
        for a in &idx.anchors {
            *counts.entry(a.h).or_insert(0) += 1;
        }
        let expected_repeated = counts.values().filter(|&&c| c > 1).count() as u64;
        assert_eq!(idx.repeated_anchor_hashes, expected_repeated);
        assert!(expected_repeated > 0, "the tandem repeat must be detected");
        for (h, loc) in &idx.anchor_map {
            assert_eq!(counts[h], 1, "a repeated hash leaked into the anchor map");
            // The recorded position must really hold a k-mer that hashes to `h`,
            // with the strand flag that section 6.2 prescribes.
            let s = &idx.contig_seq(loc.contig)[loc.pos as usize..loc.pos as usize + p.k];
            let mut fwd = 0u64;
            let mut rc = 0u64;
            for &b in s {
                fwd = (fwd << 2) | u64::from(code(b) & 3);
            }
            for &b in s.iter().rev() {
                rc = (rc << 2) | u64::from(3 - (code(b) & 3));
            }
            assert_eq!(anchor_hash(fwd.min(rc)), *h);
            assert_eq!(loc.plus, fwd <= rc);
        }
    }

    #[test]
    fn block_dictionary_finds_a_block_by_hash() {
        let seq = randomish(15_000, 21);
        let g = genome(&[&seq]);
        let p = params(21, [256, 64, 16]);
        let idx = PivotIndex::new("pivot", &g, &p);
        let b = &idx.contig_blocks(3, 0)[1];
        let hits = idx.lookup_block(3, b.hash).expect("must be present");
        assert!(hits.iter().any(|l| l.start == b.start && l.end == b.end));
        assert_eq!(
            hits.iter().find(|l| l.start == b.start).unwrap().flipped,
            b.flipped
        );
        // A hash that is not in the dictionary must not be inventable.
        assert!(idx.lookup_block(3, b.hash ^ 0xDEAD_BEEF).is_none());
    }

    #[test]
    fn a_reverse_complemented_pivot_gives_the_same_block_hashes() {
        // Property 1 of section 10.2, checked at the index level: the block
        // dictionaries must be identical for a contig and its reverse
        // complement, because block hashes are strand-independent (6.3).
        let seq = randomish(20_000, 31);
        let rc = crate::seq::revcomp(&seq);
        let p = params(21, [256, 64, 16]);
        let a = PivotIndex::new("p", &genome(&[&seq]), &p);
        let b = PivotIndex::new("p", &genome(&[&rc]), &p);
        for level in 1..=LEVELS {
            let mut ha: Vec<u64> = a.contig_blocks(level, 0).iter().map(|x| x.hash).collect();
            let mut hb: Vec<u64> = b.contig_blocks(level, 0).iter().map(|x| x.hash).collect();
            ha.sort_unstable();
            hb.sort_unstable();
            assert_eq!(ha, hb, "level {level} block hashes differ under reversal");
        }
        // The anchor hash multiset must match too.
        let mut aa: Vec<u64> = a.anchors.iter().map(|x| x.h).collect();
        let mut bb: Vec<u64> = b.anchors.iter().map(|x| x.h).collect();
        aa.sort_unstable();
        bb.sort_unstable();
        assert_eq!(aa, bb);
    }

    #[test]
    fn empty_and_tiny_contigs_do_not_panic() {
        let g = genome(&[b"", b"AC", b"AAAAAAAAAAAAAAAAAAAAAAAA"]);
        let p = params(21, [256, 64, 16]);
        let idx = PivotIndex::new("p", &g, &p);
        assert_eq!(idx.contig_blocks(1, 0).len(), 0);
        // A contig shorter than k has no anchors, so it is one block.
        assert_eq!(idx.contig_blocks(3, 1).len(), 1);
        assert_eq!(idx.total_len, 0 + 2 + 24);
    }
}
