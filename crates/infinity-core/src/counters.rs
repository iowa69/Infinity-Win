//! The live counters of sections 7 and 11.
//!
//! Every counter is a plain integer and is only ever incremented, so a run can
//! be reproduced bit for bit and the counters can be summed across threads in
//! any order.

use serde::{Deserialize, Serialize};

use crate::seq::LEVELS;

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Counters {
    /// Leaf or interior blocks looked up in a pivot dictionary (6.6 step A).
    pub blocks_looked_up: u64,
    /// Blocks that hit, split by level. Index 0 is level 1.
    pub blocks_matched: [u64; LEVELS],
    /// Pieces that had to be compared base by base (6.6 step C).
    pub pieces_compared: u64,
    /// Pieces resolved from the cache instead of being aligned (6.10).
    pub cache_hits: u64,
    /// Dynamic-programming alignments actually run (6.7).
    pub alignments_run: u64,
    /// Bases that went through those alignments.
    pub bases_aligned: u64,
    /// k-mers inspected while scanning.
    pub kmers_scanned: u64,
    /// Genomes scanned.
    pub genomes_scanned: u64,
    /// Ungapped extensions at contig ends (6.6 step D).
    pub end_extensions: u64,
    /// Shared points dropped as non-collinear (6.6 step B).
    pub stray_points_dropped: u64,
    /// Breakpoints recorded (6.6 step C).
    pub breakpoints: u64,

    // --- the number the whole idea is judged on ---
    /// Blocks that matched the pivot and were skipped, summed over all levels.
    /// Always exactly `blocks_matched.iter().sum()` — what
    /// [`Counters::total_matched`] recomputes — because `compare::match_blocks`
    /// increments the two together. It is *not* `blocks_looked_up`, which counts
    /// every block examined whether it matched or not (measured on a 50-genome
    /// run: 38155 against 42175). Redundant, but kept because an external
    /// report script reads this key out of `run.json`.
    pub blocks_total: u64,
}

impl Counters {
    pub fn merge(&mut self, other: &Counters) {
        self.blocks_looked_up += other.blocks_looked_up;
        for i in 0..LEVELS {
            self.blocks_matched[i] += other.blocks_matched[i];
        }
        self.pieces_compared += other.pieces_compared;
        self.cache_hits += other.cache_hits;
        self.alignments_run += other.alignments_run;
        self.bases_aligned += other.bases_aligned;
        self.kmers_scanned += other.kmers_scanned;
        self.genomes_scanned += other.genomes_scanned;
        self.end_extensions += other.end_extensions;
        self.stray_points_dropped += other.stray_points_dropped;
        self.breakpoints += other.breakpoints;
        self.blocks_total += other.blocks_total;
    }

    /// Fraction of block lookups that were identical, in `[0, 1]`.
    /// `None` when nothing was looked up, so a ratio is never invented.
    pub fn block_identity(&self) -> Option<f64> {
        if self.blocks_looked_up == 0 {
            return None;
        }
        let matched: u64 = self.blocks_matched.iter().sum();
        Some(matched as f64 / self.blocks_looked_up as f64)
    }

    /// True when the run was exact and used every part of the idea. Used by the
    /// benchmark to label a configuration.
    pub fn total_matched(&self) -> u64 {
        self.blocks_matched.iter().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_is_additive_and_order_independent() {
        let a = Counters {
            blocks_looked_up: 10,
            blocks_matched: [1, 2, 3],
            pieces_compared: 4,
            ..Default::default()
        };
        let b = Counters {
            blocks_looked_up: 5,
            blocks_matched: [0, 1, 1],
            pieces_compared: 2,
            ..Default::default()
        };
        let mut x = a.clone();
        x.merge(&b);
        let mut y = b.clone();
        y.merge(&a);
        assert_eq!(x, y);
        assert_eq!(x.blocks_looked_up, 15);
        assert_eq!(x.blocks_matched, [1, 3, 4]);
    }

    #[test]
    fn identity_is_none_when_nothing_was_looked_up() {
        assert_eq!(Counters::default().block_identity(), None);
        let c = Counters {
            blocks_looked_up: 200,
            blocks_matched: [0, 0, 199],
            ..Default::default()
        };
        assert_eq!(c.block_identity(), Some(0.995));
    }
}
