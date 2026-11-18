//! Sketches and the sketch distance of section 6.4.
//!
//! The sketch of a genome is the set of hashes of its anchors at
//! `sketch_level`. The distance is the Mash formula applied to that set, which
//! is why every distance the program reports carries its focus level: a sketch
//! distance at Focus 0 is a *sample* estimate, not a count of differences.

/// Number of hashes shared by two sorted, de-duplicated sketches.
pub fn shared_count(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j) = (0, 0);
    let mut shared = 0;
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            shared += 1;
            i += 1;
            j += 1;
        } else if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    shared
}

/// Jaccard containment `J = s / (|A| + |B| - s)`.
pub fn jaccard(shared: usize, len_a: usize, len_b: usize) -> f64 {
    let denom = len_a + len_b - shared;
    if denom == 0 {
        return 0.0;
    }
    shared as f64 / denom as f64
}

/// Sketch distance of 6.4:
///
/// ```text
/// D = 1                              if J = 0
/// D = -(1/k) * ln( 2J / (1 + J) )    otherwise
/// ```
///
/// The `2J / (1 + J)` term is the containment-corrected form: it estimates the
/// fraction of shared k-mers as if both genomes had the same accessory content,
/// which is what makes the distance usable across genomes of different size.
/// The result is clamped to `[0, 1]`; `ln` of a number above 1 would be
/// negative and a negative distance is meaningless.
pub fn sketch_distance(shared: usize, len_a: usize, len_b: usize, k: usize) -> f64 {
    let j = jaccard(shared, len_a, len_b);
    if j <= 0.0 {
        return 1.0;
    }
    let p = 2.0 * j / (1.0 + j);
    if p >= 1.0 {
        return 0.0;
    }
    let d = -(p.ln()) / k as f64;
    d.clamp(0.0, 1.0)
}

/// Distance between two sketches.
pub fn distance(a: &[u64], b: &[u64], k: usize) -> f64 {
    sketch_distance(shared_count(a, b), a.len(), b.len(), k)
}

/// Distances from one genome to many, given the shared counts already computed
/// through the pivot index. Kept separate so the clustering loop can reuse it.
pub fn distance_from_shared(shared: usize, len_a: usize, len_b: usize, k: usize) -> f64 {
    sketch_distance(shared, len_a, len_b, k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_count_on_sorted_deduped_lists() {
        let a = vec![1, 3, 5, 7, 9];
        let b = vec![2, 3, 5, 8];
        assert_eq!(shared_count(&a, &b), 2);
        assert_eq!(shared_count(&a, &a), 5);
        assert_eq!(shared_count(&a, &[]), 0);
    }

    #[test]
    fn identical_sketches_are_at_distance_zero() {
        let a = vec![10, 20, 30];
        assert_eq!(distance(&a, &a, 31), 0.0);
    }

    #[test]
    fn disjoint_sketches_are_at_distance_one() {
        assert_eq!(distance(&[1, 2, 3], &[4, 5, 6], 31), 1.0);
        assert_eq!(distance(&[], &[], 31), 1.0);
    }

    #[test]
    fn distance_grows_monotonically_as_sharing_falls() {
        let full: Vec<u64> = (0..1000).collect();
        let k = 31;
        // Walk from full sharing down to none: the distance must only grow.
        let mut prev = -1.0;
        for shared in (0..=10usize).rev().map(|i| i * 100) {
            let other: Vec<u64> = (0..shared as u64)
                .chain(10_000..10_000 + (1000 - shared) as u64)
                .collect();
            assert_eq!(shared_count(&full, &other), shared);
            let d = distance(&full, &other, k);
            assert!(d >= prev, "distance fell as sharing fell: {d} after {prev}");
            prev = d;
        }
        assert!(
            prev > 0.99,
            "no sharing at all must give distance 1, got {prev}"
        );
    }

    #[test]
    fn a_single_snp_is_barely_visible_in_a_sketch() {
        // Section 1: one SNP changes 21 of about 5 million k-mers, so a sketch
        // of a few thousand hashes notices it only rarely. This is the reason
        // Focus 0 can never be exact, and the test documents that.
        let k = 31;
        let a: Vec<u64> = (0..2500u64).collect();
        // One k-mer changed. The replacement is sorted into place so the sketch
        // stays sorted, which `shared_count` requires.
        let b: Vec<u64> = (1..2500u64).chain(std::iter::once(999_999)).collect();
        assert_eq!(shared_count(&a, &b), 2499);
        let d = distance(&a, &b, k);
        assert!(d > 0.0);
        assert!(
            d < 1e-4,
            "one SNP must move a sketch distance by almost nothing, got {d}"
        );
    }

    #[test]
    fn jaccard_handles_empty_input() {
        assert_eq!(jaccard(0, 0, 0), 0.0);
        assert_eq!(jaccard(2, 3, 3), 0.5);
    }
}
