//! The mutable shape of one genome during a simulation.
//!
//! A genome is not stored as a `Vec<u8>` while it is being evolved: that would
//! make every branch of the pedigree copy megabytes, and it would lose the link
//! back to reference coordinates that the truth files are written in. Instead it
//! is a *spine*: an ordered list of reference ranges and inserted segments, plus
//! a sorted list of substitutions in reference coordinates. Rendering to bases
//! happens once, at the leaves.
//!
//! Cloning a spine costs the size of its event list, not the size of the genome,
//! which is what makes a 50-leaf pedigree on a 15 Mb reference cheap.

/// One stretch of a genome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Piece {
    /// Reference positions `[start, end)`, still in reference coordinates.
    Ref { start: usize, end: usize },
    /// Sequence the simulator invented, with no reference coordinate.
    Ins { seq: Vec<u8> },
}

/// An ordered list of pieces: the genome without its substitutions.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct Spine {
    pieces: Vec<Piece>,
}

impl Spine {
    /// A spine holding the whole reference, `[0, len)`.
    pub(crate) fn whole(len: usize) -> Self {
        let pieces = if len == 0 {
            Vec::new()
        } else {
            vec![Piece::Ref { start: 0, end: len }]
        };
        Self { pieces }
    }

    /// Total number of bases this spine renders to, substitutions aside.
    pub(crate) fn total_len(&self) -> usize {
        self.pieces.iter().map(piece_len).sum()
    }

    /// Drop reference positions `[start, end)`.
    ///
    /// The caller guarantees - through [`crate::positions::FreeSpace`] - that the
    /// range is entirely present and contains no inserted segment, so this only
    /// ever trims or splits `Ref` pieces.
    pub(crate) fn delete(&mut self, start: usize, end: usize) {
        if start >= end {
            return;
        }
        let mut out: Vec<Piece> = Vec::with_capacity(self.pieces.len() + 1);
        for p in std::mem::take(&mut self.pieces) {
            match p {
                Piece::Ins { .. } => out.push(p),
                Piece::Ref { start: s, end: e } => {
                    if e <= start || s >= end {
                        out.push(Piece::Ref { start: s, end: e });
                    } else {
                        if s < start {
                            out.push(Piece::Ref {
                                start: s,
                                end: start,
                            });
                        }
                        if end < e {
                            out.push(Piece::Ref { start: end, end: e });
                        }
                    }
                }
            }
        }
        self.pieces = out;
    }

    /// Insert `seq` immediately before reference position `pos`.
    ///
    /// `pos` must be present in the spine; it stays present, the new sequence
    /// simply sits in front of it.
    pub(crate) fn insert_before(&mut self, pos: usize, seq: Vec<u8>) {
        let mut at = self.pieces.len();
        for (i, p) in self.pieces.iter().enumerate() {
            if let Piece::Ref { end, .. } = p {
                if pos < *end {
                    at = i;
                    break;
                }
            }
        }
        let span = match self.pieces.get(at) {
            Some(Piece::Ref { start, end }) => Some((*start, *end)),
            _ => None,
        };
        match span {
            Some((s, e)) if pos > s => {
                self.pieces[at] = Piece::Ref { start: s, end: pos };
                self.pieces.insert(at + 1, Piece::Ins { seq });
                self.pieces
                    .insert(at + 2, Piece::Ref { start: pos, end: e });
            }
            _ => self.pieces.insert(at, Piece::Ins { seq }),
        }
    }
}

fn piece_len(p: &Piece) -> usize {
    match p {
        Piece::Ref { start, end } => end - start,
        Piece::Ins { seq } => seq.len(),
    }
}

/// A substitution in reference coordinates: position and the base that replaces
/// the reference base there. Lists of these are always sorted and never contain
/// a position twice, which [`crate::positions::FreeSpace`] guarantees.
pub(crate) type Substitutions = [(usize, u8)];

/// Merge two sorted, disjoint substitution lists.
pub(crate) fn merge_subs(a: &Substitutions, b: &Substitutions) -> Vec<(usize, u8)> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        if a[i].0 <= b[j].0 {
            out.push(a[i]);
            i += 1;
        } else {
            out.push(b[j]);
            j += 1;
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Render a spine to bases.
///
/// `a` and `b` are the two sorted substitution lists that apply to it - the
/// cluster root's divergence and this branch line's own SNPs. A substitution
/// whose position is no longer present is dropped, which cannot happen given the
/// occupancy rules but must not corrupt the output if a caller ever relaxes
/// them.
pub(crate) fn render(
    spine: &Spine,
    ref_seq: &[u8],
    a: &Substitutions,
    b: &Substitutions,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(spine.total_len());
    let (mut ia, mut ib) = (0usize, 0usize);
    for piece in &spine.pieces {
        match piece {
            Piece::Ins { seq } => out.extend_from_slice(seq),
            Piece::Ref { start, end } => {
                while ia < a.len() && a[ia].0 < *start {
                    ia += 1;
                }
                while ib < b.len() && b[ib].0 < *start {
                    ib += 1;
                }
                for (offset, &base) in ref_seq[*start..*end].iter().enumerate() {
                    let p = *start + offset;
                    if ia < a.len() && a[ia].0 == p {
                        out.push(a[ia].1);
                        ia += 1;
                    } else if ib < b.len() && b[ib].0 == p {
                        out.push(b[ib].1);
                        ib += 1;
                    } else {
                        out.push(base);
                    }
                }
            }
        }
    }
    out
}

/// Render reference positions `[start, end)` with substitutions applied. Used to
/// record what a deletion actually removed.
pub(crate) fn render_ref_range(
    ref_seq: &[u8],
    start: usize,
    end: usize,
    a: &Substitutions,
    b: &Substitutions,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(end.saturating_sub(start));
    let mut ia = lower_bound(a, start);
    let mut ib = lower_bound(b, start);
    let from = start.min(ref_seq.len());
    let to = end.min(ref_seq.len());
    for (offset, &base) in ref_seq[from..to].iter().enumerate() {
        let p = from + offset;
        if ia < a.len() && a[ia].0 == p {
            out.push(a[ia].1);
            ia += 1;
        } else if ib < b.len() && b[ib].0 == p {
            out.push(b[ib].1);
            ib += 1;
        } else {
            out.push(base);
        }
    }
    out
}

fn lower_bound(list: &Substitutions, pos: usize) -> usize {
    list.partition_point(|&(p, _)| p < pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_seq(len: usize) -> Vec<u8> {
        (0..len).map(|i| b"ACGT"[i % 4]).collect()
    }

    #[test]
    fn whole_spine_renders_the_reference() {
        let r = ref_seq(100);
        let s = Spine::whole(100);
        assert_eq!(render(&s, &r, &[], &[]), r);
        assert_eq!(s.total_len(), 100);
    }

    #[test]
    fn empty_spine_renders_nothing() {
        let s = Spine::whole(0);
        assert_eq!(render(&s, &[], &[], &[]), Vec::<u8>::new());
        assert_eq!(s.total_len(), 0);
    }

    #[test]
    fn substitutions_are_applied_from_both_lists() {
        let r = ref_seq(40);
        let s = Spine::whole(40);
        let root = vec![(3, b'T'), (20, b'A')];
        let own = vec![(10, b'C'), (21, b'G')];
        let got = render(&s, &r, &root, &own);
        assert_eq!(got.len(), 40);
        for (i, &b) in got.iter().enumerate() {
            let want = match i {
                3 => b'T',
                20 => b'A',
                10 => b'C',
                21 => b'G',
                _ => r[i],
            };
            assert_eq!(b, want, "position {i}");
        }
    }

    #[test]
    fn merge_subs_keeps_order() {
        let a = vec![(1, b'A'), (5, b'C'), (9, b'G')];
        let b = vec![(0, b'T'), (5, b'A'), (20, b'C')];
        let m = merge_subs(&a, &b);
        assert_eq!(
            m.iter().map(|&(p, _)| p).collect::<Vec<_>>(),
            vec![0, 1, 5, 5, 9, 20]
        );
        assert!(a.windows(2).all(|w| w[0].0 < w[1].0));
    }

    #[test]
    fn deletion_splits_the_spine() {
        let r = ref_seq(100);
        let mut s = Spine::whole(100);
        s.delete(30, 40);
        assert_eq!(s.total_len(), 90);
        let got = render(&s, &r, &[], &[]);
        let mut want = r[..30].to_vec();
        want.extend_from_slice(&r[40..]);
        assert_eq!(got, want);
        // A second deletion inside the first piece.
        s.delete(10, 15);
        let got = render(&s, &r, &[], &[]);
        let mut want = r[..10].to_vec();
        want.extend_from_slice(&r[15..30]);
        want.extend_from_slice(&r[40..]);
        assert_eq!(got, want);
        assert_eq!(s.total_len(), 85);
    }

    #[test]
    fn deletion_at_both_ends() {
        let r = ref_seq(50);
        let mut s = Spine::whole(50);
        s.delete(0, 5);
        s.delete(45, 50);
        assert_eq!(render(&s, &r, &[], &[]), r[5..45].to_vec());
        s.delete(0, 50);
        assert_eq!(s.total_len(), 0);
        assert!(render(&s, &r, &[], &[]).is_empty());
    }

    #[test]
    fn insertions_sit_before_their_anchor() {
        let r = ref_seq(20);
        let mut s = Spine::whole(20);
        s.insert_before(5, b"NNN".to_vec());
        s.insert_before(0, b"X".to_vec());
        s.insert_before(20, b"YY".to_vec());
        let got = render(&s, &r, &[], &[]);
        let mut want = vec![b'X'];
        want.extend_from_slice(&r[0..5]);
        want.extend_from_slice(b"NNN");
        want.extend_from_slice(&r[5..20]);
        want.extend_from_slice(b"YY");
        assert_eq!(got, want);
        assert_eq!(s.total_len(), 26);
    }

    #[test]
    fn insertions_survive_a_neighbouring_deletion() {
        let r = ref_seq(30);
        let mut s = Spine::whole(30);
        s.insert_before(10, b"AAAA".to_vec());
        s.delete(10, 15);
        let got = render(&s, &r, &[], &[]);
        let mut want = r[..10].to_vec();
        want.extend_from_slice(b"AAAA");
        want.extend_from_slice(&r[15..]);
        assert_eq!(got, want);
    }

    #[test]
    fn render_drops_a_substitution_whose_position_is_gone() {
        let r = ref_seq(30);
        let mut s = Spine::whole(30);
        s.delete(10, 20);
        // Cannot happen through FreeSpace, but must not corrupt the output.
        let subs = vec![(5, b'Z'), (12, b'Y'), (25, b'X')];
        let got = render(&s, &r, &subs, &[]);
        let mut want = r[..10].to_vec();
        want.extend_from_slice(&r[20..]);
        want[5] = b'Z';
        want[15] = b'X';
        assert_eq!(got, want);
    }

    #[test]
    fn render_ref_range_applies_substitutions() {
        let r = ref_seq(50);
        let subs = vec![(10, b'Z'), (14, b'Q')];
        let mut got = render_ref_range(&r, 10, 15, &subs, &[]);
        assert_eq!(got, vec![b'Z', r[11], r[12], r[13], b'Q']);
        got = render_ref_range(&r, 0, 3, &subs, &[]);
        assert_eq!(got, r[..3]);
        assert!(render_ref_range(&r, 20, 20, &subs, &[]).is_empty());
    }
}
