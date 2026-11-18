//! Occupancy of reference positions inside one cluster's pedigree.
//!
//! The whole accuracy argument of section 10.1 rests on one invariant: **no two
//! events in a cluster ever touch the same reference position.** When it holds,
//! the true distance between two genomes is exactly the sum of the SNP counts on
//! the path between them in the truth tree, because nothing can overwrite,
//! revert or delete another event's work.
//!
//! Two bitsets enforce it:
//!
//! * `snp_used` - every position substituted by any SNP, including the ones that
//!   turned the reference into the cluster root. A branch SNP may not land here,
//!   which keeps a leaf's mutation list disjoint from its root's.
//! * `event_used` - every position touched by an event on the cluster's tree: a
//!   branch SNP, an insertion anchor, or a deleted range. A deletion window may
//!   not contain one of these, so a loss can never swallow a branch SNP or an
//!   inserted segment belonging to any member of the cluster.
//!
//! Root divergence positions deliberately stay out of `event_used`: at 2%
//! divergence they sit every 50 bases, and forbidding them would make segment
//! loss impossible. They are shared by the whole cluster, so they cancel out of
//! every within-cluster pairwise distance and a loss covering them only costs
//! coverage, never a SNP.

use mtato_core::seq::code;

use crate::rng::Rng;

/// A growable bitset over `[0, size)`.
#[derive(Debug, Clone)]
struct PositionSet {
    words: Vec<u64>,
    size: usize,
}

impl PositionSet {
    fn new(size: usize) -> Self {
        Self {
            words: vec![0; size.div_ceil(64)],
            size,
        }
    }

    fn insert(&mut self, pos: usize) {
        if pos < self.size {
            self.words[pos >> 6] |= 1u64 << (pos & 63);
        }
    }

    fn insert_range(&mut self, start: usize, end: usize) {
        let end = end.min(self.size);
        if start >= end {
            return;
        }
        let first = start >> 6;
        let last = (end - 1) >> 6;
        let low_mask = u64::MAX << (start & 63);
        // A shift by 64 is undefined, so build the high mask from the count of
        // bits to keep instead.
        let kept = end & 63;
        let high_mask = if kept == 0 {
            u64::MAX
        } else {
            (1u64 << kept) - 1
        };
        if first == last {
            self.words[first] |= low_mask & high_mask;
            return;
        }
        self.words[first] |= low_mask;
        for w in &mut self.words[first + 1..last] {
            *w = u64::MAX;
        }
        self.words[last] |= high_mask;
    }

    fn contains(&self, pos: usize) -> bool {
        pos < self.size && (self.words[pos >> 6] >> (pos & 63)) & 1 == 1
    }

    /// True when any position of `[start, end)` is set.
    fn any_in_range(&self, start: usize, end: usize) -> bool {
        let end = end.min(self.size);
        if start >= end {
            return false;
        }
        let first = start >> 6;
        let last = (end - 1) >> 6;
        let low_mask = u64::MAX << (start & 63);
        let kept = end & 63;
        let high_mask = if kept == 0 {
            u64::MAX
        } else {
            (1u64 << kept) - 1
        };
        if first == last {
            return self.words[first] & low_mask & high_mask != 0;
        }
        if self.words[first] & low_mask != 0 || self.words[last] & high_mask != 0 {
            return true;
        }
        self.words[first + 1..last].iter().any(|&w| w != 0)
    }

    /// Maximal runs of clear positions, as half-open ranges.
    ///
    /// Only used as the fallback when rejection sampling has run out of luck,
    /// so it favours clarity over speed: fully set and fully clear words are
    /// skipped wholesale and only the rest is walked bit by bit.
    fn free_runs(&self) -> Vec<(usize, usize)> {
        let mut runs: Vec<(usize, usize)> = Vec::new();
        let mut open: Option<usize> = None;
        let last = self.words.len().saturating_sub(1);
        for (i, &w) in self.words.iter().enumerate() {
            let mut inv = !w;
            if i == last {
                let spare = self.words.len() * 64 - self.size;
                if spare > 0 {
                    inv &= u64::MAX >> spare;
                }
            }
            let base = i * 64;
            if inv == 0 {
                if let Some(s) = open.take() {
                    runs.push((s, base));
                }
                continue;
            }
            let filled = if i == last { self.size - base } else { 64 };
            if inv == u64::MAX >> (64 - filled) {
                if open.is_none() {
                    open = Some(base);
                }
                continue;
            }
            for b in 0..filled {
                if (inv >> b) & 1 == 1 {
                    if open.is_none() {
                        open = Some(base + b);
                    }
                } else if let Some(s) = open.take() {
                    runs.push((s, base + b));
                }
            }
        }
        if let Some(s) = open {
            runs.push((s, self.size));
        }
        runs
    }
}

/// The positions a cluster's pedigree has already committed to.
#[derive(Debug, Clone)]
pub(crate) struct FreeSpace {
    size: usize,
    snp_used: PositionSet,
    event_used: PositionSet,
}

impl FreeSpace {
    /// An empty occupancy map over `[0, size)`.
    pub(crate) fn new(size: usize) -> Self {
        Self {
            size,
            snp_used: PositionSet::new(size),
            event_used: PositionSet::new(size),
        }
    }

    /// Pick `count` distinct free positions, returned in ascending order.
    ///
    /// `require_acgt` drops positions whose reference base is N: the engine
    /// never calls an N a SNP, so counting one as truth would be a lie. When
    /// `block_events` is set the positions also go into `event_used`, which is
    /// what keeps later deletion windows away from them.
    ///
    /// Positions are claimed as they are drawn, so a single call can never
    /// return the same one twice - that is the whole exclusivity invariant.
    ///
    /// Rejection sampling handles the normal case in a couple of draws. If it
    /// stalls 256 times in a row - a small reference, a high divergence, most of
    /// it already used - the free positions are enumerated once and drawn from
    /// exactly. Failing that there is genuinely no room and the caller gets an
    /// error rather than a silently shortened branch.
    pub(crate) fn take_positions(
        &mut self,
        rng: &mut Rng,
        ref_seq: &[u8],
        count: usize,
        require_acgt: bool,
        block_events: bool,
    ) -> Result<Vec<usize>, String> {
        let mut out: Vec<usize> = Vec::with_capacity(count);
        let mut misses = 0usize;
        while out.len() < count {
            if misses >= 256 {
                let mut pool = self.free_pool(ref_seq, require_acgt);
                if pool.len() < count - out.len() {
                    return Err(format!(
                        "not enough free reference positions: {} needed, {} available",
                        count - out.len(),
                        pool.len()
                    ));
                }
                while out.len() < count {
                    let p = pool.swap_remove(rng.below_usize(pool.len()));
                    self.mark(p, block_events);
                    out.push(p);
                }
                break;
            }
            let p = rng.below_usize(self.size);
            if self.is_free(p) && (!require_acgt || code(ref_seq[p]) < 4) {
                self.mark(p, block_events);
                out.push(p);
                misses = 0;
            } else {
                misses += 1;
            }
        }
        out.sort_unstable();
        Ok(out)
    }

    /// Pick the start of a window of `len` consecutive positions that no event
    /// of this cluster has touched, and claim it. `None` when there is no room.
    ///
    /// Because every deleted range and every insertion anchor is in
    /// `event_used`, a window that passes is guaranteed to be entirely present
    /// in the genome and to contain no inserted segment.
    pub(crate) fn take_window(&mut self, rng: &mut Rng, len: usize) -> Option<usize> {
        if len == 0 || len > self.size {
            return None;
        }
        let limit = self.size - len;
        for _ in 0..256 {
            let start = rng.below_usize(limit + 1);
            if !self.event_used.any_in_range(start, start + len) {
                self.event_used.insert_range(start, start + len);
                return Some(start);
            }
        }
        // Rejection has stalled: walk the free runs and take a window inside one.
        let runs: Vec<(usize, usize)> = self
            .event_used
            .free_runs()
            .into_iter()
            .filter(|&(s, e)| e - s >= len)
            .collect();
        if runs.is_empty() {
            return None;
        }
        let total: usize = runs.iter().map(|&(s, e)| e - s - len + 1).sum();
        let mut pick = rng.below_usize(total);
        for &(s, e) in &runs {
            let span = e - s - len + 1;
            if pick < span {
                let start = s + pick;
                self.event_used.insert_range(start, start + len);
                return Some(start);
            }
            pick -= span;
        }
        None
    }

    fn is_free(&self, pos: usize) -> bool {
        !self.snp_used.contains(pos) && !self.event_used.contains(pos)
    }

    fn mark(&mut self, pos: usize, block_events: bool) {
        self.snp_used.insert(pos);
        if block_events {
            self.event_used.insert(pos);
        }
    }

    fn free_pool(&self, ref_seq: &[u8], require_acgt: bool) -> Vec<usize> {
        let mut pool = Vec::new();
        for &(s, e) in &self.snp_used.free_runs() {
            for (offset, &base) in ref_seq[s..e].iter().enumerate() {
                let p = s + offset;
                if !self.event_used.contains(p) && (!require_acgt || code(base) < 4) {
                    pool.push(p);
                }
            }
        }
        pool
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    fn set_with(size: usize, ranges: &[(usize, usize)]) -> PositionSet {
        let mut s = PositionSet::new(size);
        for &(a, b) in ranges {
            s.insert_range(a, b);
        }
        s
    }

    fn naive(s: &PositionSet) -> Vec<bool> {
        (0..s.size).map(|p| s.contains(p)).collect()
    }

    #[test]
    fn range_insert_matches_single_inserts() {
        let cases: Vec<(usize, &[(usize, usize)])> = vec![
            (200, &[(0, 1)]),
            (200, &[(0, 64)]),
            (200, &[(0, 65)]),
            (200, &[(63, 129)]),
            (200, &[(10, 20), (100, 200)]),
            (200, &[(199, 200)]),
            (200, &[(0, 200)]),
            (64, &[(0, 64)]),
            (65, &[(1, 65)]),
            (1, &[(0, 1)]),
        ];
        for (size, ranges) in cases {
            let fast = set_with(size, ranges);
            let mut slow = PositionSet::new(size);
            for &(a, b) in ranges {
                for p in a..b.min(size) {
                    slow.insert(p);
                }
            }
            assert_eq!(naive(&fast), naive(&slow), "size {size} ranges {ranges:?}");
        }
    }

    #[test]
    fn any_in_range_agrees_with_contains() {
        let s = set_with(300, &[(5, 9), (64, 66), (128, 200)]);
        for start in 0..300 {
            for end in start..300 {
                let want = (start..end).any(|p| s.contains(p));
                assert_eq!(s.any_in_range(start, end), want, "[{start},{end})");
            }
        }
    }

    #[test]
    fn free_runs_partition_the_clear_positions() {
        let s = set_with(300, &[(0, 3), (62, 70), (295, 300)]);
        let runs = s.free_runs();
        assert_eq!(runs, vec![(3, 62), (70, 295)]);
        let total: usize = runs.iter().map(|&(a, b)| b - a).sum();
        assert_eq!(total, 300 - 3 - 8 - 5);
    }

    #[test]
    fn free_runs_handles_all_clear_and_all_set() {
        assert_eq!(PositionSet::new(130).free_runs(), vec![(0, 130)]);
        assert!(set_with(130, &[(0, 130)]).free_runs().is_empty());
        assert_eq!(PositionSet::new(0).free_runs(), vec![]);
    }

    #[test]
    fn take_positions_returns_distinct_free_values() {
        let mut rng = Rng::from_seed(1);
        let ref_seq = vec![b'A'; 1000];
        let mut space = FreeSpace::new(1000);
        let got = space
            .take_positions(&mut rng, &ref_seq, 100, true, true)
            .unwrap();
        assert_eq!(got.len(), 100);
        let mut sorted = got.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), 100);
        for &p in &got {
            assert!(space.snp_used.contains(p));
            assert!(space.event_used.contains(p));
        }
        let again = space
            .take_positions(&mut rng, &ref_seq, 50, true, false)
            .unwrap();
        assert!(again.iter().all(|p| !got.contains(p)));
        // block_events = false leaves event_used alone.
        assert!(again.iter().all(|&p| !space.event_used.contains(p)));
    }

    #[test]
    fn take_positions_skips_n_bases() {
        let mut rng = Rng::from_seed(2);
        let mut ref_seq = vec![b'A'; 500];
        ref_seq.iter_mut().take(250).for_each(|b| *b = b'N');
        let mut space = FreeSpace::new(500);
        let got = space
            .take_positions(&mut rng, &ref_seq, 200, true, true)
            .unwrap();
        assert!(got.iter().all(|&p| p >= 250));
        assert_eq!(got.len(), 200);
    }

    #[test]
    fn take_positions_falls_back_to_enumeration_when_room_runs_out() {
        let mut rng = Rng::from_seed(3);
        let ref_seq = vec![b'C'; 300];
        let mut space = FreeSpace::new(300);
        // Fill all but 10 positions; rejection sampling cannot find them.
        let first = space
            .take_positions(&mut rng, &ref_seq, 290, true, false)
            .unwrap();
        assert_eq!(first.len(), 290);
        let last = space
            .take_positions(&mut rng, &ref_seq, 10, true, false)
            .unwrap();
        assert_eq!(last.len(), 10);
        let mut all = first;
        all.extend(last);
        all.sort_unstable();
        assert_eq!(all, (0..300).collect::<Vec<_>>());
        assert!(space
            .take_positions(&mut rng, &ref_seq, 1, true, false)
            .is_err());
    }

    #[test]
    fn take_window_avoids_claimed_positions() {
        let mut rng = Rng::from_seed(4);
        let mut space = FreeSpace::new(100_000);
        space.event_used.insert_range(1000, 1100);
        let mut taken = Vec::new();
        for _ in 0..50 {
            let start = space.take_window(&mut rng, 500).unwrap();
            assert!(start + 500 <= 100_000);
            let overlaps_forbidden = start < 1100 && 1000 < start + 500;
            assert!(
                !overlaps_forbidden,
                "window at {start} overlaps the seeded band"
            );
            taken.push(start);
        }
        // Claiming is what makes the next window disjoint from this one.
        taken.sort_unstable();
        assert!(taken.windows(2).all(|w| w[0] + 500 <= w[1]));
        assert!(FreeSpace::new(10).take_window(&mut rng, 11).is_none());
    }

    #[test]
    fn take_window_finds_a_tight_gap_by_enumeration() {
        let mut rng = Rng::from_seed(5);
        let mut space = FreeSpace::new(1000);
        space.event_used.insert_range(0, 400);
        space.event_used.insert_range(500, 1000);
        let start = space.take_window(&mut rng, 100).unwrap();
        assert_eq!(start, 400);
        // The only gap is now taken.
        assert!(space.take_window(&mut rng, 100).is_none());
    }
}
