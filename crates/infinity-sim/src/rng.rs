//! Deterministic pseudo-random generation (section 10.1).
//!
//! The simulator must be reproducible from a seed alone, so it carries its own
//! generator instead of depending on `rand`: no external version bump can
//! silently change a benchmark's truth files, and no platform-dependent
//! entropy source can reach the output.
//!
//! The generator is xoshiro256\*\* seeded from SplitMix64, the standard pairing.
//! Every consumer works from a *derived stream* ([`stream`]), so the reference,
//! each cluster root, each cluster pedigree and each genome's fragmentation are
//! independent: changing `contigs` does not change a single base of the
//! underlying sequences, which the fragmentation test relies on.

use infinity_core::seq::fmix64;

/// Stream tag: the random reference sequence.
pub(crate) const TAG_REFERENCE: u64 = 1;
/// Stream tag: the substitutions that turn the reference into a cluster root.
pub(crate) const TAG_ROOT: u64 = 2;
/// Stream tag: a cluster's tree topology and branch lengths.
pub(crate) const TAG_TREE: u64 = 3;
/// Stream tag: the mutations placed on a cluster's branches.
pub(crate) const TAG_EVENTS: u64 = 4;
/// Stream tag: cutting one genome into contigs.
pub(crate) const TAG_FRAGMENT: u64 = 5;

/// SplitMix64, used only to expand one seed word into the xoshiro state.
fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// xoshiro256\*\* generator.
#[derive(Debug, Clone)]
pub(crate) struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Build a generator from one 64-bit seed.
    pub(crate) fn from_seed(seed: u64) -> Self {
        let mut x = seed;
        let mut s = [0u64; 4];
        for v in s.iter_mut() {
            *v = splitmix64(&mut x);
        }
        // xoshiro has one absorbing state; SplitMix64 can reach it from
        // exactly one seed, so patch it rather than lose reproducibility.
        if s.iter().all(|&v| v == 0) {
            s[0] = 0x9E37_79B9_7F4A_7C15;
        }
        Self { s }
    }

    /// Next 64-bit word.
    pub(crate) fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// A uniform value in `[0, n)`, without modulo bias (Lemire's method).
    /// Panics when `n` is zero, which is always a caller bug.
    pub(crate) fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "Rng::below(0)");
        let mut x = self.next_u64();
        let mut prod = u128::from(x) * u128::from(n);
        let mut low = prod as u64;
        if low < n {
            let threshold = n.wrapping_neg() % n;
            while low < threshold {
                x = self.next_u64();
                prod = u128::from(x) * u128::from(n);
                low = prod as u64;
            }
        }
        (prod >> 64) as u64
    }

    /// A uniform value in `[0, n)`.
    pub(crate) fn below_usize(&mut self, n: usize) -> usize {
        self.below(n as u64) as usize
    }

    /// A uniform value in `[lo, hi]`.
    pub(crate) fn range_inclusive(&mut self, lo: usize, hi: usize) -> usize {
        debug_assert!(lo <= hi);
        lo + self.below_usize(hi - lo + 1)
    }

    /// A fair coin flip, taken from the high bit.
    pub(crate) fn coin(&mut self) -> bool {
        (self.next_u64() >> 63) != 0
    }

    /// A uniform value in `[0, 1)`, with 53 bits of resolution.
    pub(crate) fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0f64 / (1u64 << 53) as f64)
    }

    /// A uniform value in `(0, 1]`, the form the Poisson samplers need.
    fn unit_open(&mut self) -> f64 {
        1.0 - self.unit()
    }

    /// Fisher-Yates shuffle.
    pub(crate) fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.below_usize(i + 1);
            items.swap(i, j);
        }
    }

    /// A Poisson draw with the given mean.
    ///
    /// Small means use Knuth's product-of-uniforms method; larger means use the
    /// equivalent exponential inter-arrival count, which stays exact but costs
    /// `O(mean)` draws. `SimConfig::validate` keeps the mean far below the
    /// point where that would matter.
    pub(crate) fn poisson(&mut self, mean: f64) -> usize {
        if !mean.is_finite() || mean <= 0.0 {
            return 0;
        }
        if mean < 30.0 {
            let limit = (-mean).exp();
            let mut k = 0usize;
            let mut prod = 1.0f64;
            loop {
                prod *= self.unit_open();
                if prod <= limit {
                    return k;
                }
                k += 1;
            }
        }
        let mut k = 0usize;
        let mut acc = 0.0f64;
        loop {
            acc -= self.unit_open().ln();
            if acc > mean {
                return k;
            }
            k += 1;
        }
    }
}

/// An independent generator for one purpose, derived from the run seed.
///
/// `tag` says what the stream is for and `index` which object of that kind it
/// belongs to (a cluster number, a genome number). Both are mixed through
/// `fmix64` so that neighbouring tags and indices do not produce correlated
/// streams.
pub(crate) fn stream(seed: u64, tag: u64, index: u64) -> Rng {
    let mixed = tag.wrapping_mul(0x0100_0000_01b3).wrapping_add(index);
    Rng::from_seed(fmix64(seed ^ fmix64(mixed)))
}

/// `count` distinct values from `[lo, hi)`, sorted ascending (Floyd's
/// algorithm). Exact and `O(count)`, unlike rejection sampling, which degrades
/// once `count` approaches the size of the range - the simulator needs it for
/// cut points, where a genome of 100 contigs leaves little slack.
pub(crate) fn distinct_sorted(rng: &mut Rng, lo: usize, hi: usize, count: usize) -> Vec<usize> {
    let n = hi.saturating_sub(lo);
    let count = count.min(n);
    let mut out: Vec<usize> = Vec::with_capacity(count);
    for j in (n - count)..n {
        let candidate = lo + rng.below_usize(j + 1);
        match out.binary_search(&candidate) {
            // Already present: Floyd's rule takes `lo + j` instead, which keeps
            // the vector sorted for free because every earlier draw was below it.
            Ok(_) => out.push(lo + j),
            Err(at) => out.insert(at, candidate),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_stream() {
        let mut a = stream(7, TAG_TREE, 3);
        let mut b = stream(7, TAG_TREE, 3);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_tag_or_index_diverges() {
        assert_ne!(
            stream(7, TAG_TREE, 3).next_u64(),
            stream(7, TAG_TREE, 4).next_u64()
        );
        assert_ne!(
            stream(7, TAG_TREE, 3).next_u64(),
            stream(7, TAG_EVENTS, 3).next_u64()
        );
        assert_ne!(
            stream(7, TAG_TREE, 3).next_u64(),
            stream(8, TAG_TREE, 3).next_u64()
        );
    }

    #[test]
    fn below_stays_in_range_and_spreads() {
        let mut rng = Rng::from_seed(11);
        let mut seen = [0usize; 7];
        for _ in 0..7000 {
            let v = rng.below(7) as usize;
            assert!(v < 7);
            seen[v] += 1;
        }
        // A fair draw over 7 values, 1000 expected each; allow a wide margin.
        for (i, &c) in seen.iter().enumerate() {
            assert!((700..1300).contains(&c), "bucket {i} got {c}");
        }
    }

    #[test]
    fn poisson_matches_its_mean() {
        let mut rng = Rng::from_seed(5);
        for &mean in &[0.5f64, 3.0, 12.0, 40.0, 500.0] {
            let draws = 4000;
            let total: usize = (0..draws).map(|_| rng.poisson(mean)).sum();
            let got = total as f64 / draws as f64;
            let tol = (6.0 * mean.sqrt() / (draws as f64).sqrt()).max(0.05);
            assert!((got - mean).abs() < tol, "mean {mean}: got {got}");
        }
        let mut rng = Rng::from_seed(5);
        assert_eq!(rng.poisson(0.0), 0);
        assert_eq!(rng.poisson(-1.0), 0);
    }

    #[test]
    fn poisson_is_exactly_zero_for_tiny_means_often_enough() {
        let mut rng = Rng::from_seed(3);
        let zeros = (0..10_000).filter(|_| rng.poisson(0.5) == 0).count();
        // exp(-0.5) = 0.6065
        assert!((5900..6250).contains(&zeros), "zeros = {zeros}");
    }

    #[test]
    fn distinct_sorted_covers_the_whole_range_when_asked() {
        let mut rng = Rng::from_seed(2);
        let v = distinct_sorted(&mut rng, 0, 5, 5);
        assert_eq!(v, vec![0, 1, 2, 3, 4]);
        let v = distinct_sorted(&mut rng, 10, 20, 4);
        assert_eq!(v.len(), 4);
        assert!(v.windows(2).all(|w| w[0] < w[1]));
        assert!(v.iter().all(|&x| (10..20).contains(&x)));
        // Asking for more than exists saturates instead of panicking.
        assert_eq!(distinct_sorted(&mut rng, 0, 3, 99).len(), 3);
    }

    #[test]
    fn distinct_sorted_works_when_the_range_is_tight() {
        let mut rng = Rng::from_seed(9);
        for _ in 0..50 {
            let v = distinct_sorted(&mut rng, 1, 100, 99);
            assert_eq!(v, (1..100).collect::<Vec<_>>());
        }
    }

    #[test]
    fn shuffle_keeps_the_multiset() {
        let mut rng = Rng::from_seed(4);
        let mut v: Vec<usize> = (0..20).collect();
        let original = v.clone();
        rng.shuffle(&mut v);
        let mut sorted = v.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, original);
    }

    #[test]
    fn unit_is_in_the_half_open_interval() {
        let mut rng = Rng::from_seed(6);
        for _ in 0..1000 {
            let u = rng.unit();
            assert!((0.0..1.0).contains(&u));
        }
    }
}
