//! Interval arithmetic on half-open `[start, end)` ranges.
//!
//! Coordinates are 0-based half-open inside the program (section 5). Coverage
//! (6.8) and the `comparable` denominator of the pairwise distance (6.9) are
//! both set operations on intervals, and both must be exact: an off-by-one
//! here changes every distance in the run.
//!
//! Every function that takes slices expects them **merged**: sorted, disjoint
//! and non-adjacent. [`merge`] produces that form.

use serde::{Deserialize, Serialize};

/// A half-open interval.
pub type Iv = (u64, u64);

/// Sort intervals and merge everything that overlaps or touches, so that
/// `[0,5)` and `[5,9)` become one interval `[0,9)`.
pub fn merge(mut v: Vec<Iv>) -> Vec<Iv> {
    if v.is_empty() {
        return v;
    }
    v.sort_unstable();
    let mut out: Vec<Iv> = Vec::with_capacity(v.len());
    let mut cur = v[0];
    for &(s, e) in v.iter().skip(1) {
        if s <= cur.1 {
            if e > cur.1 {
                cur.1 = e;
            }
        } else {
            out.push(cur);
            cur = (s, e);
        }
    }
    out.push(cur);
    out
}

/// Total length covered by a merged interval list.
pub fn total_len(merged: &[Iv]) -> u64 {
    merged.iter().map(|&(s, e)| e.saturating_sub(s)).sum()
}

/// Length of the union of two merged interval lists.
pub fn union_len(a: &[Iv], b: &[Iv]) -> u64 {
    let mut both = Vec::with_capacity(a.len() + b.len());
    both.extend_from_slice(a);
    both.extend_from_slice(b);
    total_len(&merge(both))
}

/// Length of the intersection of two merged interval lists.
pub fn intersection_len(a: &[Iv], b: &[Iv]) -> u64 {
    let (mut i, mut j) = (0, 0);
    let mut total = 0u64;
    while i < a.len() && j < b.len() {
        let s = a[i].0.max(b[j].0);
        let e = a[i].1.min(b[j].1);
        if s < e {
            total += e - s;
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    total
}

/// `from` minus `remove`, both merged. Result is merged.
pub fn subtract(from: &[Iv], remove: &[Iv]) -> Vec<Iv> {
    let mut out: Vec<Iv> = Vec::new();
    let mut j = 0usize;
    for &(mut s, e) in from {
        while j < remove.len() && remove[j].1 <= s {
            j += 1;
        }
        let mut k = j;
        let mut cur = s;
        while k < remove.len() && remove[k].0 < e {
            let (rs, re) = remove[k];
            if rs > cur {
                out.push((cur, rs.min(e)));
            }
            if re > cur {
                cur = re;
            }
            if cur >= e {
                break;
            }
            k += 1;
        }
        if cur < e {
            out.push((cur, e));
        }
        s = e; // silence unused-assignment lint; `s` is not reused
        let _ = s;
    }
    merge(out)
}

/// Is `pos` inside a merged interval list?
pub fn contains(merged: &[Iv], pos: u64) -> bool {
    // Binary search for the last interval starting at or before pos.
    let idx = merged.partition_point(|&(s, _)| s <= pos);
    if idx == 0 {
        return false;
    }
    pos < merged[idx - 1].1
}

/// Merge two already-merged lists into one merged list.
pub fn union(a: &[Iv], b: &[Iv]) -> Vec<Iv> {
    let mut both = Vec::with_capacity(a.len() + b.len());
    both.extend_from_slice(a);
    both.extend_from_slice(b);
    merge(both)
}

/// A merged interval list that can be serialised into `run.json` and events.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct IntervalSet {
    pub iv: Vec<Iv>,
}

impl IntervalSet {
    pub fn new(iv: Vec<Iv>) -> Self {
        Self { iv: merge(iv) }
    }
    pub fn len(&self) -> u64 {
        total_len(&self.iv)
    }
    pub fn is_empty(&self) -> bool {
        self.iv.is_empty()
    }
    pub fn contains(&self, pos: u64) -> bool {
        contains(&self.iv, pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_joins_touching_intervals() {
        assert_eq!(
            merge(vec![(5, 9), (0, 5), (20, 21)]),
            vec![(0, 9), (20, 21)]
        );
        assert_eq!(merge(vec![(0, 10), (2, 3)]), vec![(0, 10)]);
        assert!(merge(vec![]).is_empty());
    }

    #[test]
    fn lengths_are_exact() {
        let a = merge(vec![(0, 10), (20, 30)]);
        let b = merge(vec![(5, 25)]);
        assert_eq!(total_len(&a), 20);
        assert_eq!(union_len(&a, &b), 30);
        assert_eq!(intersection_len(&a, &b), 10);
    }

    #[test]
    fn subtract_handles_every_overlap_shape() {
        let from = merge(vec![(0, 100)]);
        assert_eq!(
            subtract(&from, &merge(vec![(20, 30), (50, 60)])),
            vec![(0, 20), (30, 50), (60, 100)]
        );
        assert_eq!(subtract(&from, &merge(vec![(0, 100)])), vec![]);
        assert_eq!(subtract(&from, &merge(vec![(100, 200)])), vec![(0, 100)]);
        assert_eq!(subtract(&from, &merge(vec![])), vec![(0, 100)]);
        let two = merge(vec![(0, 10), (20, 30)]);
        assert_eq!(
            subtract(&two, &merge(vec![(5, 25)])),
            vec![(0, 5), (25, 30)]
        );
    }

    #[test]
    fn contains_is_half_open() {
        let a = merge(vec![(10, 20), (30, 40)]);
        assert!(!contains(&a, 9));
        assert!(contains(&a, 10));
        assert!(contains(&a, 19));
        assert!(!contains(&a, 20));
        assert!(contains(&a, 30));
        assert!(!contains(&a, 40));
    }

    #[test]
    fn union_len_agrees_with_union_then_total() {
        let a = merge(vec![(0, 10), (20, 30)]);
        let b = merge(vec![(5, 25), (40, 50)]);
        assert_eq!(union_len(&a, &b), total_len(&union(&a, &b)));
    }
}
