//! Comparing a member to its pivot (6.6), and the piece cache (6.10).
//!
//! The four steps of 6.6 are implemented as separate functions over one contig
//! at a time, so that each can be tested on its own:
//!
//! - step A, [`match_blocks`]: the block shortcut. Blocks whose hash the pivot
//!   already knows are identical and are skipped without looking inside.
//! - step B, [`shared_points`]: give every anchor of the contig a pivot
//!   position, from the block mapping or from the pivot's unique-anchor map.
//! - step C, [`build_pieces`]: walk the shared points and turn the gaps between
//!   them into pieces, or into breakpoints when they are not collinear.
//! - step D, [`extend_contig_ends`]: ungapped x-drop extension past the last
//!   shared point at each contig end.
//!
//! Focus 1 is the output of A-D: a list of difference pieces, the breakpoints
//! and the unplaced bases, with no base-level work done. Focus 2 resolves each
//! difference piece with [`crate::variant::diff`].

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use serde::{Deserialize, Serialize};
use xxhash_rust::xxh3::xxh3_64;

use crate::align::{walked_range, xdrop_walk, AlignScores};
use crate::blocks::{block_ranges, hash_block};
use crate::counters::Counters;
use crate::fasta::GenomeData;
use crate::interval::{merge, subtract, total_len, Iv};
use crate::params::Parameters;
use crate::pivot::{BlockLoc, PivotIndex};
use crate::seq::{revcomp, scan_windows, Anchor, LEVELS};
use crate::variant::{flags, DiffParams, DiffResult, RawVariant, Variant, VariantKind};

/// A shared point (section 5): an anchor of a member genome whose position on
/// the pivot is known.
///
/// `gpos` and `ppos` are the starts of the same k-base anchor on the member
/// contig and on the pivot's forward strand. `rel` is true when the member runs
/// opposite to the pivot there.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SharedPoint {
    pub gpos: u64,
    pub ppos: u64,
    pub pcontig: u32,
    pub rel: bool,
    /// True when the position came from a matched block rather than from the
    /// anchor map. Only anchor-map points can be dropped as stray (6.6 B).
    pub from_block: bool,
}

/// A block of the member that was found, by hash, in the pivot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchedBlock {
    pub level: u8,
    /// Member contig interval.
    pub gs: u64,
    pub ge: u64,
    /// Pivot contig interval.
    pub pcontig: u32,
    pub ps: u64,
    pub pe: u64,
    pub rel: bool,
    /// The pivot holds this block hash more than once, so the chosen copy may
    /// be the wrong one (section 15).
    pub repeat: bool,
}

impl MatchedBlock {
    /// Length, which is the same on both sides.
    pub fn len(&self) -> u64 {
        self.ge - self.gs
    }

    /// True when the block spans no bases. A matched block always covers at
    /// least `k` bases, so this is the defensive partner of [`Self::len`]
    /// rather than a case the engine produces.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Map a member base position inside this block to the pivot (6.6 A).
    pub fn map_base(&self, x: u64) -> u64 {
        if self.rel {
            self.pe - 1 - (x - self.gs)
        } else {
            self.ps + (x - self.gs)
        }
    }

    /// Map the start of a member anchor inside this block to the start of the
    /// corresponding anchor on the pivot's forward strand.
    ///
    /// This differs from [`Self::map_base`] by `k - 1` in the reversed case:
    /// the member k-mer `[x, x+k)` occupies pivot bases `[pe-k, pe)` when
    /// reversed, so its forward-strand start is `pe - k`, not `pe - 1`.
    pub fn map_anchor(&self, x: u64, k: u64) -> u64 {
        if self.rel {
            self.pe - k - (x - self.gs)
        } else {
            self.ps + (x - self.gs)
        }
    }
}

/// A difference piece: a stretch where the member and the pivot differ, known
/// to differ, but not yet resolved (the Focus 1 record of 6.6).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiffPiece {
    pub gcontig: u32,
    /// Member interval `[gs, ge)`.
    pub gs: u64,
    pub ge: u64,
    pub pcontig: u32,
    /// Pivot interval `[ps, pe)`, always on the pivot's forward strand.
    pub ps: u64,
    pub pe: u64,
    /// True when the member runs opposite to the pivot across this piece.
    pub rel: bool,
    /// Hash of the member (alt) sequence, which with the pivot interval is the
    /// piece-cache key of 6.10.
    pub alt_hash: u64,
    /// Flags inherited from the piece's context: `REPEAT` when a bounding block
    /// matched a repetitive pivot hash.
    pub flags: u32,
}

impl DiffPiece {
    pub fn cache_key(&self) -> PieceKey {
        (self.alt_hash, self.pcontig, self.ps, self.pe)
    }

    /// The member-side sequence of this piece, oriented onto the pivot's
    /// forward strand.
    pub fn alt(&self, member_contig: &[u8]) -> Vec<u8> {
        let s = &member_contig[self.gs as usize..self.ge as usize];
        if self.rel {
            revcomp(s)
        } else {
            s.to_vec()
        }
    }
}

/// Piece-cache key: `(hash of alt, pivot contig, pivot start, pivot end)`,
/// exactly as 6.10 specifies. The pivot interval determines `ref` and the hash
/// determines `alt`, so together they determine the variants.
pub type PieceKey = (u64, u32, u64, u64);

/// A stretch of the pivot that one stretch of the member explains.
///
/// Both coordinate pairs are kept, and that is what makes the `dup` rule of 6.8
/// computable at all. Blocks overlap their neighbours by exactly `k` bases by
/// design (6.3), so two coverage intervals that overlap on the pivot are
/// normally just the seam between two consecutive blocks of the *same* part of
/// the member. A genuine duplication is when two places of the member that do
/// not overlap each other map onto the same place on the pivot. Without the
/// member coordinates those two cases are indistinguishable, and treating every
/// seam as a duplication silently deletes about 3% of the genome from the
/// comparable denominator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoverSeg {
    /// Pivot contig and interval.
    pub pcontig: u32,
    pub ps: u64,
    pub pe: u64,
    /// Member contig and interval that produced it.
    pub gcontig: u32,
    pub gs: u64,
    pub ge: u64,
}

impl CoverSeg {
    /// The pivot interval, in pivot-contig-local coordinates.
    pub fn pivot_iv(&self) -> Iv {
        (self.ps, self.pe)
    }

    /// True when the member intervals of two segments are disjoint, which is
    /// what makes their pivot overlap a real duplication rather than a seam.
    fn member_disjoint(&self, other: &Self) -> bool {
        self.gcontig != other.gcontig || self.ge <= other.gs || other.ge <= self.gs
    }
}

/// A place where the member's path through the pivot breaks: the two shared
/// points on either side are not collinear (6.6 C).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Breakpoint {
    pub gcontig: u32,
    pub gpos: u64,
    pub pcontig: u32,
    pub ppos: u64,
}

/// Resolved variants of one piece, shared through the cache.
#[derive(Debug, Clone)]
struct CachedPiece {
    variants: Vec<RawVariant>,
    n_positions: Vec<usize>,
    align_failed: bool,
}

/// The piece cache of 6.10.
///
/// A mutation shared by 500 members of a clade is then resolved once. The cache
/// is shared across the worker threads of one cluster; hits are counted because
/// the hit rate is one of the numbers the experiment of section 11 reports.
#[derive(Debug, Default)]
pub struct PieceCache {
    inner: Mutex<HashMap<PieceKey, Arc<CachedPiece>>>,
    enabled: bool,
    hits: AtomicU64,
}

impl PieceCache {
    pub fn new(enabled: bool) -> Self {
        Self {
            inner: Mutex::new(HashMap::new()),
            enabled,
            hits: AtomicU64::new(0),
        }
    }

    fn get(&self, key: &PieceKey) -> Option<Arc<CachedPiece>> {
        if !self.enabled {
            return None;
        }
        let guard = self.inner.lock().ok()?;
        let got = guard.get(key).cloned();
        drop(guard);
        if got.is_some() {
            self.hits.fetch_add(1, Ordering::Relaxed);
        }
        got
    }

    pub fn insert(&self, key: PieceKey, value: DiffResult) {
        if !self.enabled {
            return;
        }
        if let Ok(mut guard) = self.inner.lock() {
            guard.entry(key).or_insert_with(|| {
                Arc::new(CachedPiece {
                    variants: value.variants.clone(),
                    n_positions: value.n_positions.clone(),
                    align_failed: value.align_failed,
                })
            });
        }
    }

    pub fn hits(&self) -> u64 {
        self.hits.load(Ordering::Relaxed)
    }

    /// Look a piece up and return its resolved content as owned data.
    ///
    /// This exists so that callers outside this module can use the cache
    /// without the private `CachedPiece` type leaking into their signatures.
    pub fn cached(&self, key: &PieceKey) -> Option<(Vec<RawVariant>, Vec<usize>, bool)> {
        let hit = self.get(key)?;
        Some((
            hit.variants.clone(),
            hit.n_positions.clone(),
            hit.align_failed,
        ))
    }

    pub fn len(&self) -> usize {
        self.inner.lock().map(|g| g.len()).unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// The Focus 1 picture of one member: what differs, without base-level work.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Focus1 {
    pub pieces: Vec<DiffPiece>,
    pub breakpoints: Vec<Breakpoint>,
    /// Matched blocks, kept because coverage (6.8) needs them.
    pub matched: Vec<MatchedBlock>,
    /// Member intervals that could not be placed on the pivot, as
    /// `(contig, start, end)`.
    pub unplaced: Vec<(u32, u64, u64)>,
    /// Pivot intervals covered by a matched block or a compared piece, before
    /// the exclusions of 6.8. Each carries the member interval that produced it.
    pub covered: Vec<CoverSeg>,
    /// Pivot positions where nothing could be called because a base was N on
    /// either side, as `(pivot contig, start, end)`.
    pub not_called: Vec<(u32, u64, u64)>,
    pub novel_bp: u64,
    pub counters: Counters,
}

/// Everything a comparison of one member against its pivot produces.
#[derive(Debug, Clone)]
pub struct MemberComparison {
    pub genome_id: u32,
    pub focus1: Focus1,
    /// Resolved variants in pivot coordinates, sorted and de-duplicated.
    pub variants: Vec<Variant>,
    /// Pivot intervals this member does not cover (6.8), in concatenated pivot
    /// coordinates, merged.
    pub uncovered: Vec<Iv>,
    /// Pivot intervals covered by two or more places of the member.
    pub duplicated: Vec<Iv>,
    pub novel_bp: u64,
    pub counters: Counters,
}

impl MemberComparison {
    pub fn snps(&self) -> u64 {
        self.variants
            .iter()
            .filter(|v| v.kind == VariantKind::Snp)
            .count() as u64
    }

    pub fn indels(&self) -> u64 {
        self.variants.iter().filter(|v| v.kind.is_indel()).count() as u64
    }

    pub fn uncovered_bp(&self) -> u64 {
        total_len(&self.uncovered)
    }
}

// ---------------------------------------------------------------------------
// Step A: the block shortcut
// ---------------------------------------------------------------------------

/// A block whose hash the pivot knows, before the copy to use is chosen.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    level: usize,
    gs: u64,
    ge: u64,
    flipped: bool,
    hash: u64,
}

/// Find the blocks of one member contig that the pivot already knows (6.6 A).
///
/// In flat mode only leaf (level-3) blocks are used. In hierarchical mode the
/// search starts at level 1 and descends only into blocks that missed, which is
/// the "zoom" of 6.10. A coarse block that matches covers exactly the same
/// bases as its leaf blocks would, so the two modes must produce identical
/// variants (property 6 of section 10.2).
pub fn match_blocks(
    pivot: &PivotIndex,
    seq: &[u8],
    anchor_pos: &[Vec<u64>; LEVELS],
    params: &Parameters,
    counters: &mut Counters,
) -> Vec<MatchedBlock> {
    if !params.block_skip {
        // Experiment switch (section 11): no block hashing at all. Every stretch
        // between shared anchors is compared directly, so nothing is matched
        // here and all the work happens in the piece step.
        return Vec::new();
    }
    let k = params.k;
    // Blocks of this contig at every level, computed once and reused as the
    // descent narrows.
    let all: Vec<Vec<(u64, u64)>> = (0..LEVELS)
        .map(|l| block_ranges(&anchor_pos[l], seq.len() as u64, k as u64))
        .collect();

    let start_level = if params.hierarchical { 1 } else { LEVELS };
    let mut to_check: Vec<(u64, u64)> = all[start_level - 1].clone();
    let mut cands: Vec<Candidate> = Vec::new();

    for level in start_level..=LEVELS {
        if to_check.is_empty() {
            break;
        }
        let mut next: Vec<(u64, u64)> = Vec::new();
        let finer = if level < LEVELS {
            Some(&all[level])
        } else {
            None
        };
        for &(s, e) in &to_check {
            counters.blocks_looked_up += 1;
            let (hash, flipped) = hash_block(&seq[s as usize..e as usize]);
            if pivot.lookup_block(level, hash).is_some() {
                counters.blocks_matched[level - 1] += 1;
                counters.blocks_total += 1;
                cands.push(Candidate {
                    level,
                    gs: s,
                    ge: e,
                    flipped,
                    hash,
                });
            } else if let Some(finer) = finer {
                // Zoom in: the level+1 blocks of this block are exactly those
                // whose anchor start falls in [s, e-k].
                let from = finer.partition_point(|&(bs, _)| bs < s);
                let hi = e.saturating_sub(k as u64);
                for &(bs, be) in &finer[from..] {
                    if bs > hi {
                        break;
                    }
                    next.push((bs, be));
                }
            }
        }
        to_check = next;
    }

    // Choose which copy of a repeated block to use, left to right, so that the
    // "continues from the previous matched block" rule of 6.6 A can be applied.
    cands.sort_unstable_by_key(|c| (c.gs, c.ge, c.level));
    // Adjacent coarse blocks share their boundary anchor, so the leaf block
    // starting there is reached twice by the descent. Keep it once.
    cands.dedup_by_key(|c| (c.gs, c.ge, c.level));
    let mut out: Vec<MatchedBlock> = Vec::with_capacity(cands.len());
    for c in &cands {
        let hits = match pivot.lookup_block(c.level, c.hash) {
            Some(h) if !h.is_empty() => h,
            _ => continue,
        };
        let repeat = hits.len() > 1;
        let loc = if !repeat {
            hits[0]
        } else {
            pick_continuation(hits, c.flipped, out.last(), k as u64).unwrap_or(hits[0])
        };
        out.push(MatchedBlock {
            level: c.level as u8,
            gs: c.gs,
            ge: c.ge,
            pcontig: loc.contig,
            ps: loc.start,
            pe: loc.end,
            rel: c.flipped ^ loc.flipped,
            repeat,
        });
    }
    out
}

/// Among several copies of a repeated block, prefer the one that continues
/// directly from the previous matched block (6.6 A).
fn pick_continuation(
    hits: &[BlockLoc],
    flipped: bool,
    prev: Option<&MatchedBlock>,
    k: u64,
) -> Option<BlockLoc> {
    let p = prev?;
    hits.iter().copied().find(|h| {
        let rel = flipped ^ h.flipped;
        if rel != p.rel || h.contig != p.pcontig {
            return false;
        }
        // The two blocks overlap by exactly k bases on the pivot, in the
        // direction the member is running.
        if rel {
            h.end == p.ps + k
        } else {
            h.start + k == p.pe
        }
    })
}

// ---------------------------------------------------------------------------
// Step B: shared points
// ---------------------------------------------------------------------------

/// Give every anchor of one member contig a pivot position (6.6 B).
pub fn shared_points(
    pivot: &PivotIndex,
    anchors: &[Anchor],
    matched: &[MatchedBlock],
    params: &Parameters,
    counters: &mut Counters,
) -> Vec<SharedPoint> {
    let k = params.k as u64;

    // Boundary anchors of matched blocks, so the lookup below is O(1) rather
    // than a scan of every block for every anchor.
    let mut boundary: HashMap<u64, &MatchedBlock> = HashMap::with_capacity(matched.len() * 2);
    for m in matched {
        boundary.entry(m.gs).or_insert(m);
        if m.ge >= k {
            boundary.entry(m.ge - k).or_insert(m);
        }
    }

    let mut raw: Vec<SharedPoint> = Vec::with_capacity(anchors.len());
    for a in anchors {
        if let Some(m) = boundary.get(&a.pos) {
            raw.push(SharedPoint {
                gpos: a.pos,
                ppos: m.map_anchor(a.pos, k),
                pcontig: m.pcontig,
                rel: m.rel,
                from_block: true,
            });
            continue;
        }
        if let Some(loc) = pivot.anchor_map.get(&a.h) {
            raw.push(SharedPoint {
                gpos: a.pos,
                ppos: loc.pos,
                pcontig: loc.contig,
                // rel = strand_G XOR strand_P
                rel: a.plus != loc.plus,
                from_block: false,
            });
        }
    }
    raw.sort_unstable_by_key(|p| p.gpos);
    raw.dedup_by_key(|p| p.gpos);
    drop_stray_points(raw, params, counters)
}

/// Drop a shared point that came from the anchor map if it is not collinear
/// with its nearest shared neighbours on both sides while those two neighbours
/// are collinear with each other (6.6 B). This removes stray matches, which are
/// the main way a repeated k-mer sneaks past the unique-anchor filter.
fn drop_stray_points(
    points: Vec<SharedPoint>,
    params: &Parameters,
    counters: &mut Counters,
) -> Vec<SharedPoint> {
    if points.len() < 3 {
        return points;
    }
    let mut kept: Vec<SharedPoint> = Vec::with_capacity(points.len());
    for (i, p) in points.iter().enumerate() {
        if p.from_block {
            kept.push(*p);
            continue;
        }
        if let (Some(l), Some(r)) = (kept.last(), points.get(i + 1)) {
            if collinear(l, r, params) && !collinear(l, p, params) && !collinear(p, r, params) {
                counters.stray_points_dropped += 1;
                continue;
            }
        }
        kept.push(*p);
    }
    kept
}

/// Collinearity of two shared points (6.6 C): same pivot contig, same `rel`,
/// pivot order consistent with `rel`, and neither span longer than
/// `max_span_ratio` times the other plus `max_span_extra`.
pub fn collinear(x: &SharedPoint, y: &SharedPoint, params: &Parameters) -> bool {
    if x.pcontig != y.pcontig || x.rel != y.rel {
        return false;
    }
    let ordered = if x.rel {
        y.ppos <= x.ppos
    } else {
        y.ppos >= x.ppos
    };
    if !ordered {
        return false;
    }
    let gspan = y.gpos.saturating_sub(x.gpos) as f64;
    let pspan = x.ppos.abs_diff(y.ppos) as f64;
    let extra = params.max_span_extra as f64;
    gspan <= params.max_span_ratio * pspan + extra && pspan <= params.max_span_ratio * gspan + extra
}

// ---------------------------------------------------------------------------
// Step C: pieces
// ---------------------------------------------------------------------------

/// One gap between consecutive shared points, classified.
#[derive(Debug, Clone, PartialEq)]
enum Gap {
    /// Inside a single matched block: already known identical.
    InsideMatchedBlock,
    /// Compared and found equal. Still covers the pivot interval.
    IdenticalPiece { cover: CoverSeg },
    /// Compared and different: a Focus 1 record.
    Piece(DiffPiece),
    /// Not collinear: a breakpoint, and the member bases between are unplaced.
    Break {
        bp: Breakpoint,
        unplaced: (u32, u64, u64),
    },
}

/// What walking the shared points of one contig produced (6.6 step C).
#[derive(Debug, Default)]
pub struct ContigPieces {
    /// Compared and found different: the Focus 1 records.
    pub pieces: Vec<DiffPiece>,
    /// Shared points that were not collinear with their neighbours.
    pub breakpoints: Vec<Breakpoint>,
    /// Member intervals that could not be placed on the pivot, as
    /// `(contig, start, end)`.
    pub unplaced: Vec<(u32, u64, u64)>,
    /// Pivot intervals this contig covers, each carrying the member interval
    /// that produced it, which is what makes the `dup` rule of 6.8 work (D12).
    pub covered: Vec<CoverSeg>,
}

/// Walk the shared points of one contig and classify every gap (6.6 C).
///
/// Returns the difference pieces, the breakpoints, the unplaced member
/// intervals and the pivot intervals covered.
pub fn build_pieces(
    pivot: &PivotIndex,
    contig: u32,
    seq: &[u8],
    points: &[SharedPoint],
    matched: &[MatchedBlock],
    params: &Parameters,
    counters: &mut Counters,
) -> ContigPieces {
    let k = params.k as u64;
    let ctx = GapCtx {
        pivot,
        contig,
        seq,
        matched,
        k,
        params,
    };
    let mut diff_pieces = Vec::new();
    let mut breakpoints = Vec::new();
    let mut unplaced = Vec::new();
    // Matched blocks contribute coverage on their own, whether or not their
    // boundary anchors produced shared points.
    let mut covered: Vec<CoverSeg> = matched
        .iter()
        .map(|m| CoverSeg {
            pcontig: m.pcontig,
            ps: m.ps,
            pe: m.pe,
            gcontig: contig,
            gs: m.gs,
            ge: m.ge,
        })
        .collect();

    for w in points.windows(2) {
        let (x, y) = (w[0], w[1]);
        if y.gpos <= x.gpos {
            continue;
        }
        match classify_gap(&ctx, &x, &y) {
            Gap::InsideMatchedBlock => {}
            Gap::IdenticalPiece { cover } => covered.push(cover),
            Gap::Piece(p) => {
                counters.pieces_compared += 1;
                covered.push(CoverSeg {
                    pcontig: p.pcontig,
                    ps: p.ps,
                    pe: p.pe,
                    gcontig: contig,
                    gs: p.gs,
                    ge: p.ge,
                });
                diff_pieces.push(p);
            }
            Gap::Break { bp, unplaced: u } => {
                counters.breakpoints += 1;
                breakpoints.push(bp);
                unplaced.push(u);
            }
        }
    }
    ContigPieces {
        pieces: diff_pieces,
        breakpoints,
        unplaced,
        covered,
    }
}

/// Everything [`classify_gap`] needs that is the same for every gap in one
/// contig. Grouped so that the function's own signature lists only the two
/// shared points that actually vary from gap to gap.
#[derive(Clone, Copy)]
struct GapCtx<'a> {
    pivot: &'a PivotIndex,
    /// Member contig index.
    contig: u32,
    /// Member contig sequence.
    seq: &'a [u8],
    matched: &'a [MatchedBlock],
    /// Anchor length: `params.k` widened to `u64` once, per contig.
    k: u64,
    params: &'a Parameters,
}

fn classify_gap(ctx: &GapCtx<'_>, x: &SharedPoint, y: &SharedPoint) -> Gap {
    let GapCtx {
        pivot,
        contig,
        seq,
        matched,
        k,
        params,
    } = *ctx;
    // "If the only thing between them is one matched block, there is nothing to
    // do." Generalised to: if the whole gap lies inside one matched block, the
    // sequence is already known to be identical. This also covers anchors that
    // fall in the middle of a coarse matched block in hierarchical mode.
    if matched.iter().any(|m| m.gs <= x.gpos && y.gpos + k <= m.ge) {
        return Gap::InsideMatchedBlock;
    }
    if !collinear(x, y, params) {
        return Gap::Break {
            bp: Breakpoint {
                gcontig: contig,
                gpos: y.gpos,
                pcontig: y.pcontig,
                ppos: y.ppos,
            },
            // The member bases strictly between the two anchors are unplaced.
            unplaced: (contig, x.gpos + k, y.gpos.max(x.gpos + k)),
        };
    }

    let gs = x.gpos;
    let ge = y.gpos + k;
    let (ps, pe) = if x.rel {
        (y.ppos, x.ppos + k)
    } else {
        (x.ppos, y.ppos + k)
    };
    let plen = pivot
        .contig_lens
        .get(x.pcontig as usize)
        .copied()
        .unwrap_or(0);
    if ge > seq.len() as u64 || pe > plen || ps >= pe {
        return Gap::InsideMatchedBlock;
    }

    let member_slice = &seq[gs as usize..ge as usize];
    let alt: Vec<u8> = if x.rel {
        revcomp(member_slice)
    } else {
        member_slice.to_vec()
    };
    let ref_slice = &pivot.contig_seq(x.pcontig)[ps as usize..pe as usize];
    if ref_slice == alt.as_slice() {
        return Gap::IdenticalPiece {
            cover: CoverSeg {
                pcontig: x.pcontig,
                ps,
                pe,
                gcontig: contig,
                gs,
                ge,
            },
        };
    }

    // A repeat next door makes the mapping suspect (section 15): flag it rather
    // than trying to correct it.
    let mut piece_flags = 0u32;
    if matched
        .iter()
        .any(|m| m.repeat && ((m.gs <= gs && gs < m.ge) || (m.gs <= y.gpos && y.gpos < m.ge)))
    {
        piece_flags |= flags::REPEAT;
    }

    Gap::Piece(DiffPiece {
        gcontig: contig,
        gs,
        ge,
        pcontig: x.pcontig,
        ps,
        pe,
        rel: x.rel,
        alt_hash: xxh3_64(&alt),
        flags: piece_flags,
    })
}

// ---------------------------------------------------------------------------
// Step D: contig ends
// ---------------------------------------------------------------------------

/// What the two contig ends contributed.
#[derive(Debug, Default)]
pub struct EndResult {
    pub variants: Vec<Variant>,
    pub covered: Vec<CoverSeg>,
    pub unplaced: Vec<(u32, u64, u64)>,
    pub not_called: Vec<(u32, u64, u64)>,
}

/// Extend outward from the shared points nearest each contig end (6.6 D).
///
/// Between a contig end and its nearest shared point there is no second anchor,
/// so there is nothing to build a piece from. The extension is ungapped and
/// stops by the x-drop rule; mismatches inside the kept part are SNPs and
/// member bases beyond it are unplaced. A contig with no shared point at all is
/// entirely unplaced.
pub fn extend_contig_ends(
    pivot: &PivotIndex,
    contig: u32,
    seq: &[u8],
    points: &[SharedPoint],
    matched: &[MatchedBlock],
    params: &Parameters,
    counters: &mut Counters,
) -> EndResult {
    let mut out = EndResult::default();
    let k = params.k as u64;
    if points.is_empty() {
        if !seq.is_empty() {
            out.unplaced.push((contig, 0, seq.len() as u64));
        }
        return out;
    }
    let scores = AlignScores::default();

    // --- left end ---
    let first = points[0];
    let (m_at, p_at, rel) = left_frontier(&first, matched, k);
    if m_at > 0 {
        counters.end_extensions += 1;
        let pres = pivot.contig_seq(first.pcontig);
        // The member always walks backwards here; the pivot walks backwards when
        // the orientations agree and forwards when they do not.
        let p_dir: i64 = if rel { 1 } else { -1 };
        let p_start: i64 = if rel { p_at as i64 } else { p_at as i64 - 1 };
        let res = xdrop_walk(
            pres,
            seq,
            p_start,
            m_at as i64 - 1,
            p_dir,
            -1,
            rel,
            &scores,
            params.xdrop,
        );
        for &i in &res.mismatches {
            let ppos = p_start + i as i64 * p_dir;
            let gpos = m_at as i64 - 1 - i as i64;
            if ppos < 0 || gpos < 0 {
                continue;
            }
            let (rb, ab) = (pres[ppos as usize], seq[gpos as usize]);
            if rb == b'N' || ab == b'N' {
                out.not_called
                    .push((first.pcontig, ppos as u64, ppos as u64 + 1));
                continue;
            }
            out.variants.push(Variant {
                kind: VariantKind::Snp,
                pivot_contig: first.pcontig,
                pivot_start: ppos as u64,
                pivot_end: ppos as u64 + 1,
                ref_allele: (rb as char).to_string(),
                alt_allele: (ab as char).to_string(),
                genome_contig: contig,
                genome_start: gpos as u64,
                genome_end: gpos as u64 + 1,
                strand: if rel { '-' } else { '+' },
                flags: flags::END_EXT,
            });
        }
        let (cs, ce) = walked_range(p_start, p_dir, res.len);
        if ce > cs {
            out.covered.push(CoverSeg {
                pcontig: first.pcontig,
                ps: cs.max(0) as u64,
                pe: ce as u64,
                gcontig: contig,
                gs: m_at.saturating_sub(res.len as u64),
                ge: m_at,
            });
        }
        let front = m_at.saturating_sub(res.len as u64);
        if front > 0 {
            out.unplaced.push((contig, 0, front));
        }
    }

    // --- right end ---
    let last = *points.last().unwrap();
    let (m_at, p_at, rel) = right_frontier(&last, matched, k, seq.len() as u64);
    if m_at < seq.len() as u64 {
        counters.end_extensions += 1;
        let pres = pivot.contig_seq(last.pcontig);
        let p_dir: i64 = if rel { -1 } else { 1 };
        let p_start: i64 = if rel { p_at as i64 - 1 } else { p_at as i64 };
        let res = xdrop_walk(
            pres,
            seq,
            p_start,
            m_at as i64,
            p_dir,
            1,
            rel,
            &scores,
            params.xdrop,
        );
        for &i in &res.mismatches {
            let ppos = p_start + i as i64 * p_dir;
            let gpos = m_at as i64 + i as i64;
            if ppos < 0 || ppos >= pres.len() as i64 || gpos >= seq.len() as i64 {
                continue;
            }
            let (rb, ab) = (pres[ppos as usize], seq[gpos as usize]);
            if rb == b'N' || ab == b'N' {
                out.not_called
                    .push((last.pcontig, ppos as u64, ppos as u64 + 1));
                continue;
            }
            out.variants.push(Variant {
                kind: VariantKind::Snp,
                pivot_contig: last.pcontig,
                pivot_start: ppos as u64,
                pivot_end: ppos as u64 + 1,
                ref_allele: (rb as char).to_string(),
                alt_allele: (ab as char).to_string(),
                genome_contig: contig,
                genome_start: gpos as u64,
                genome_end: gpos as u64 + 1,
                strand: if rel { '-' } else { '+' },
                flags: flags::END_EXT,
            });
        }
        let (cs, ce) = walked_range(p_start, p_dir, res.len);
        if ce > cs && cs >= 0 {
            out.covered.push(CoverSeg {
                pcontig: last.pcontig,
                ps: cs as u64,
                pe: (ce as u64).min(pres.len() as u64),
                gcontig: contig,
                gs: m_at,
                ge: m_at + res.len as u64,
            });
        }
        let after = m_at + res.len as u64;
        if after < seq.len() as u64 {
            out.unplaced.push((contig, after, seq.len() as u64));
        }
    }
    out
}

/// Where the explained part of a contig begins, as `(member position, pivot
/// position, rel)`.
///
/// If the first shared point is the far end of a matched block, the block
/// already explains everything up to its own start, so the extension must begin
/// there rather than at the anchor.
fn left_frontier(first: &SharedPoint, matched: &[MatchedBlock], k: u64) -> (u64, u64, bool) {
    let covering = matched
        .iter()
        .filter(|m| m.gs <= first.gpos && first.gpos < m.ge)
        .min_by_key(|m| m.gs);
    match covering {
        Some(m) => (m.gs, m.pe, m.rel),
        None => (
            first.gpos,
            if first.rel {
                first.ppos + k
            } else {
                first.ppos
            },
            first.rel,
        ),
    }
}

/// The same, at the right-hand end of a contig. The member position returned is
/// the first base that is not yet explained.
fn right_frontier(
    last: &SharedPoint,
    matched: &[MatchedBlock],
    k: u64,
    contig_len: u64,
) -> (u64, u64, bool) {
    let anchor_end = (last.gpos + k).min(contig_len);
    let covering = matched
        .iter()
        .filter(|m| m.gs < anchor_end && anchor_end <= m.ge)
        .max_by_key(|m| m.ge);
    match covering {
        Some(m) => (m.ge, m.pe, m.rel),
        None => (
            anchor_end,
            if last.rel { last.ppos } else { last.ppos + k },
            last.rel,
        ),
    }
}

// ---------------------------------------------------------------------------
// The whole comparison of one member
// ---------------------------------------------------------------------------

/// Compare one member genome against its pivot index (6.6), and optionally
/// resolve every difference piece to base level (6.7).
///
/// Set `resolve` to false to stop at Focus 1: the difference pieces are found
/// but not aligned, which is what lets the staged run of 6.11 show a tree
/// before the exact work is done.
pub fn compare_member(
    pivot: &PivotIndex,
    member: &GenomeData,
    genome_id: u32,
    params: &Parameters,
    cache: &PieceCache,
    resolve: bool,
) -> MemberComparison {
    let mut counters = Counters::default();
    let mut focus1 = Focus1::default();
    let mut variants: Vec<Variant> = Vec::new();

    for (ci, contig) in member.contigs.iter().enumerate() {
        let cid = ci as u32;
        if contig.seq.is_empty() {
            continue;
        }
        let (anchors, anchor_pos) = scan_contig(&contig.seq, cid, params, &pivot.th);
        counters.kmers_scanned += contig.seq.len().saturating_sub(params.k) as u64 + 1;

        let matched = match_blocks(pivot, &contig.seq, &anchor_pos, params, &mut counters);
        let points = shared_points(pivot, &anchors, &matched, params, &mut counters);
        let ContigPieces {
            pieces: diff_pieces,
            breakpoints,
            unplaced,
            covered,
        } = build_pieces(
            pivot,
            cid,
            &contig.seq,
            &points,
            &matched,
            params,
            &mut counters,
        );
        let ends = extend_contig_ends(
            pivot,
            cid,
            &contig.seq,
            &points,
            &matched,
            params,
            &mut counters,
        );

        variants.extend(ends.variants);
        focus1.covered.extend(covered);
        focus1.covered.extend(ends.covered);
        focus1.not_called.extend(ends.not_called);
        focus1.unplaced.extend(unplaced);
        focus1.unplaced.extend(ends.unplaced);
        focus1.breakpoints.extend(breakpoints);
        focus1.matched.extend(matched);

        if resolve {
            for p in &diff_pieces {
                let alt = p.alt(&contig.seq);
                let res = resolve_piece(pivot, p, &alt, params, cache, &mut counters);
                for rv in &res.variants {
                    variants.push(piece_variant(p, rv));
                }
                for &npos in &res.n_positions {
                    focus1
                        .not_called
                        .push((p.pcontig, p.ps + npos as u64, p.ps + npos as u64 + 1));
                }
            }
        }
        focus1.pieces.extend(diff_pieces);
    }

    // Novel bp: member bases that are unplaced or inserted (6.8).
    let unplaced_bp: u64 = focus1
        .unplaced
        .iter()
        .map(|&(_, s, e)| e.saturating_sub(s))
        .sum();
    let inserted_bp: u64 = variants
        .iter()
        .filter(|v| v.kind == VariantKind::Insertion)
        .map(|v| v.alt_allele.len() as u64)
        .sum();
    focus1.novel_bp = unplaced_bp + inserted_bp;

    variants.sort_by(|a, b| {
        (a.pivot_contig, a.pivot_start, a.pivot_end, a.kind as u8)
            .cmp(&(b.pivot_contig, b.pivot_start, b.pivot_end, b.kind as u8))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });
    variants.dedup();

    focus1.counters = counters.clone();
    let novel_bp = focus1.novel_bp;
    let cov = coverage(pivot, &focus1, &variants);

    MemberComparison {
        genome_id,
        focus1,
        variants,
        uncovered: cov.uncovered,
        duplicated: cov.duplicated,
        novel_bp,
        counters,
    }
}

/// Scan one member contig, returning its anchors and the anchor start positions
/// at each of the three levels.
pub fn scan_contig(
    seq: &[u8],
    contig: u32,
    params: &Parameters,
    th: &crate::seq::Thresholds,
) -> (Vec<Anchor>, [Vec<u64>; LEVELS]) {
    let mut anchors: Vec<Anchor> = Vec::new();
    let mut anchor_pos: [Vec<u64>; LEVELS] = [Vec::new(), Vec::new(), Vec::new()];
    scan_windows(seq, params.k, |pos, h, plus| {
        if let Some(level) = th.level(h) {
            anchors.push(Anchor {
                h,
                contig,
                pos: pos as u64,
                plus,
                level,
            });
            for l in level as usize..=LEVELS {
                anchor_pos[l - 1].push(pos as u64);
            }
        }
    });
    (anchors, anchor_pos)
}

/// Resolve one difference piece into variants, through the cache (6.10).
pub fn resolve_piece(
    pivot: &PivotIndex,
    piece: &DiffPiece,
    alt: &[u8],
    params: &Parameters,
    cache: &PieceCache,
    counters: &mut Counters,
) -> DiffResult {
    let key = piece.cache_key();
    if let Some(hit) = cache.get(&key) {
        counters.cache_hits += 1;
        return DiffResult {
            variants: hit.variants.clone(),
            n_positions: hit.n_positions.clone(),
            align_failed: hit.align_failed,
            used_dp: false,
        };
    }
    let ref_slice = &pivot.contig_seq(piece.pcontig)[piece.ps as usize..piece.pe as usize];
    counters.bases_aligned += ref_slice.len().max(alt.len()) as u64;
    let result = crate::variant::diff(ref_slice, alt, &DiffParams::from(params));
    if result.used_dp {
        counters.alignments_run += 1;
    }
    cache.insert(key, result.clone());
    result
}

/// Translate a piece-local raw variant into pivot coordinates.
pub fn piece_variant(piece: &DiffPiece, rv: &RawVariant) -> Variant {
    let (gstart, gend, strand) = if piece.rel {
        // `alt` is the reverse complement of the member piece, so offset i in
        // alt sits at member position ge - 1 - i.
        (
            piece.ge.saturating_sub(rv.end as u64),
            piece.ge.saturating_sub(rv.start as u64),
            '-',
        )
    } else {
        (piece.gs + rv.start as u64, piece.gs + rv.end as u64, '+')
    };
    Variant {
        kind: rv.kind,
        pivot_contig: piece.pcontig,
        pivot_start: piece.ps + rv.start as u64,
        pivot_end: piece.ps + rv.end as u64,
        ref_allele: String::from_utf8_lossy(&rv.ref_allele).into_owned(),
        alt_allele: String::from_utf8_lossy(&rv.alt_allele).into_owned(),
        genome_contig: piece.gcontig,
        genome_start: gstart,
        genome_end: gend,
        strand,
        flags: rv.flags | piece.flags,
    }
}

// ---------------------------------------------------------------------------
// Coverage (6.8)
// ---------------------------------------------------------------------------

/// Start offset of each pivot contig in concatenated pivot coordinates.
pub fn contig_offsets(pivot: &PivotIndex) -> Vec<u64> {
    let mut out = Vec::with_capacity(pivot.contig_lens.len());
    let mut acc = 0u64;
    for &l in &pivot.contig_lens {
        out.push(acc);
        acc += l;
    }
    out
}

/// Runs of N in the pivot, per contig. A property of the pivot alone, so it is
/// computed once per cluster rather than once per member.
pub fn pivot_n_intervals(pivot: &PivotIndex) -> Vec<Vec<Iv>> {
    pivot
        .contigs
        .iter()
        .map(|c| {
            let mut out: Vec<Iv> = Vec::new();
            let mut start: Option<usize> = None;
            for (i, &b) in c.seq.iter().enumerate() {
                if b == b'N' {
                    start.get_or_insert(i);
                } else if let Some(s) = start.take() {
                    out.push((s as u64, i as u64));
                }
            }
            if let Some(s) = start {
                out.push((s as u64, c.seq.len() as u64));
            }
            out
        })
        .collect()
}

/// Pivot intervals covered by two or more places of the member, as
/// `(pivot contig, start, end)`.
///
/// Only overlaps between **disjoint member intervals** count. The `k`-base
/// overlap that consecutive blocks share by construction (6.3) is not a
/// duplication, and counting it as one would remove about 3% of every genome
/// from the comparable denominator of 6.9 — which is exactly what silently loses
/// real SNPs, because a position excluded from `comparable` is a position where
/// no difference can be counted.
pub fn duplicated_regions(segs: &[CoverSeg]) -> Vec<(u32, u64, u64)> {
    if segs.len() < 2 {
        return Vec::new();
    }
    let mut order: Vec<&CoverSeg> = segs.iter().collect();
    order.sort_by_key(|s| (s.pcontig, s.ps, s.pe));
    let mut out: Vec<(u32, u64, u64)> = Vec::new();
    // Segments still open, i.e. whose pivot end is beyond the current start.
    let mut active: Vec<&CoverSeg> = Vec::new();
    for s in order {
        active.retain(|a| a.pcontig == s.pcontig && a.pe > s.ps);
        for a in &active {
            if a.member_disjoint(s) {
                let end = a.pe.min(s.pe);
                if end > s.ps {
                    out.push((s.pcontig, s.ps, end));
                }
            }
        }
        active.push(s);
    }
    // Merge per contig; contigs are disjoint in flat coordinates anyway.
    let mut merged: Vec<(u32, u64, u64)> = Vec::new();
    out.sort_unstable();
    for (c, s, e) in out {
        match merged.last_mut() {
            Some((lc, ls, le)) if *lc == c && s <= *le => {
                if e > *le {
                    *le = e;
                }
                let _ = ls;
            }
            _ => merged.push((c, s, e)),
        }
    }
    merged
}

/// The result of applying the coverage rules of 6.8, in concatenated pivot
/// coordinates.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Coverage {
    pub uncovered: Vec<Iv>,
    pub duplicated: Vec<Iv>,
    pub covered: Vec<Iv>,
}

/// Apply 6.8: a pivot base is covered if it lies in a matched block or in a
/// compared piece, and is not inside a COMPLEX variant, an N, or an interval
/// covered by two or more places of the member.
///
/// `uncovered` is the complement of `covered` within the pivot, which is what
/// the `comparable` denominator of 6.9 needs.
pub fn coverage(pivot: &PivotIndex, focus1: &Focus1, variants: &[Variant]) -> Coverage {
    let offsets = contig_offsets(pivot);
    let total = pivot.total_len;
    let base_of = |c: u32| offsets.get(c as usize).copied().unwrap_or(0);

    // What a matched block or a compared piece explains.
    let explained: Vec<Iv> = focus1
        .covered
        .iter()
        .map(|s| (base_of(s.pcontig) + s.ps, base_of(s.pcontig) + s.pe))
        .collect();
    let duplicated: Vec<Iv> = duplicated_regions(&focus1.covered)
        .iter()
        .map(|&(c, s, e)| (base_of(c) + s, base_of(c) + e))
        .collect();
    let duplicated = merge(duplicated);

    // What must be taken back off it.
    let mut excluded: Vec<Iv> = Vec::new();
    for v in variants {
        if v.kind == VariantKind::Complex {
            excluded.push((
                base_of(v.pivot_contig) + v.pivot_start,
                base_of(v.pivot_contig) + v.pivot_end,
            ));
        }
    }
    for (ci, runs) in pivot_n_intervals(pivot).iter().enumerate() {
        let base = offsets.get(ci).copied().unwrap_or(0);
        for &(s, e) in runs {
            excluded.push((base + s, base + e));
        }
    }
    for &(c, s, e) in &focus1.not_called {
        excluded.push((base_of(c) + s, base_of(c) + e));
    }
    excluded.extend(duplicated.iter().copied());
    let excluded = merge(excluded);

    let covered = subtract(&merge(explained), &excluded);
    let uncovered = subtract(&[(0, total)], &covered);
    Coverage {
        uncovered,
        duplicated,
        covered,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matched_block_mapping_is_an_involution() {
        let m = MatchedBlock {
            level: 3,
            gs: 100,
            ge: 200,
            pcontig: 0,
            ps: 5000,
            pe: 5100,
            rel: false,
            repeat: false,
        };
        assert_eq!(m.map_base(100), 5000);
        assert_eq!(m.map_base(199), 5099);
        let r = MatchedBlock { rel: true, ..m };
        assert_eq!(r.map_base(100), 5099);
        assert_eq!(r.map_base(199), 5000);
    }

    #[test]
    fn anchor_mapping_accounts_for_the_kmer_width() {
        let k = 21;
        let fwd = MatchedBlock {
            level: 3,
            gs: 0,
            ge: 100,
            pcontig: 0,
            ps: 1000,
            pe: 1100,
            rel: false,
            repeat: false,
        };
        assert_eq!(fwd.map_anchor(0, k), 1000);
        assert_eq!(fwd.map_anchor(100 - k, k), 1100 - k);
        let rev = MatchedBlock { rel: true, ..fwd };
        assert_eq!(rev.map_anchor(0, k), 1100 - k);
        assert_eq!(rev.map_anchor(100 - k, k), 1000);
        // The anchor start is the base mapping minus (k - 1) when reversed.
        assert_eq!(rev.map_anchor(0, k), rev.map_base(0) - (k - 1));
    }

    #[test]
    fn collinearity_rejects_span_explosions() {
        let p = Parameters {
            max_span_ratio: 10.0,
            max_span_extra: 50_000,
            ..Default::default()
        };
        let a = SharedPoint {
            gpos: 0,
            ppos: 0,
            pcontig: 0,
            rel: false,
            from_block: true,
        };
        let near = SharedPoint {
            gpos: 1000,
            ppos: 1000,
            ..a
        };
        assert!(collinear(&a, &near, &p));
        // 1 kb on the member against 1 Mb on the pivot: outside the ratio, and
        // outside the 50 kb slack too.
        let far = SharedPoint {
            gpos: 1000,
            ppos: 1_000_000,
            ..a
        };
        assert!(!collinear(&a, &far, &p));
        assert!(!collinear(&a, &SharedPoint { pcontig: 1, ..near }, &p));
        assert!(!collinear(&a, &SharedPoint { rel: true, ..near }, &p));
    }

    #[test]
    fn collinearity_direction_follows_rel() {
        let p = Parameters::default();
        let a = SharedPoint {
            gpos: 100,
            ppos: 1000,
            pcontig: 0,
            rel: true,
            from_block: true,
        };
        // Reversed: the member advances while the pivot retreats.
        assert!(collinear(
            &a,
            &SharedPoint {
                gpos: 200,
                ppos: 900,
                ..a
            },
            &p
        ));
        assert!(!collinear(
            &a,
            &SharedPoint {
                gpos: 200,
                ppos: 1100,
                ..a
            },
            &p
        ));
    }

    /// A coverage segment from member `[gs, ge)` to pivot `[ps, pe)`.
    fn seg(gs: u64, ge: u64, ps: u64, pe: u64) -> CoverSeg {
        CoverSeg {
            pcontig: 0,
            ps,
            pe,
            gcontig: 0,
            gs,
            ge,
        }
    }

    #[test]
    fn a_designed_block_seam_is_not_a_duplication() {
        // Consecutive blocks overlap by k on BOTH sides (6.3). That is the
        // normal case and must not be removed from coverage: doing so cost
        // about 3% of every genome and silently lost real SNPs.
        let k = 31u64;
        let segs = vec![
            seg(0, 1000, 500, 1500),
            seg(1000 - k, 2000 - k, 1500 - k, 2500 - k),
        ];
        assert!(
            duplicated_regions(&segs).is_empty(),
            "the k-base seam between two consecutive blocks is not a duplication"
        );
    }

    #[test]
    fn two_disjoint_member_places_on_one_pivot_place_is_a_duplication() {
        // The same pivot stretch explained by two parts of the member that do
        // not overlap each other: a real `dup` (6.8).
        let segs = vec![seg(0, 100, 5000, 5100), seg(9000, 9100, 5050, 5150)];
        assert_eq!(duplicated_regions(&segs), vec![(0, 5050, 5100)]);
    }

    #[test]
    fn duplication_detection_ignores_other_contigs_and_disjoint_runs() {
        let mut a = seg(0, 100, 500, 600);
        let b = seg(9000, 9100, 550, 650);
        assert_eq!(duplicated_regions(&[a, b]), vec![(0, 550, 600)]);
        // Same pivot overlap but on different member contigs is still a
        // duplication, because the member places are disjoint.
        a.gcontig = 1;
        assert_eq!(duplicated_regions(&[a, b]), vec![(0, 550, 600)]);
        // Far apart on the pivot: nothing overlaps.
        assert!(duplicated_regions(&[seg(0, 100, 0, 100), seg(200, 300, 5000, 5100)]).is_empty());
        assert!(duplicated_regions(&[seg(0, 100, 0, 100)]).is_empty());
    }

    #[test]
    fn piece_cache_returns_what_was_inserted() {
        let c = PieceCache::new(true);
        let key = (1u64, 0u32, 10u64, 20u64);
        assert!(c.get(&key).is_none());
        c.insert(
            key,
            DiffResult {
                variants: vec![RawVariant {
                    kind: VariantKind::Snp,
                    start: 0,
                    end: 1,
                    ref_allele: b"A".to_vec(),
                    alt_allele: b"T".to_vec(),
                    flags: 0,
                }],
                ..Default::default()
            },
        );
        let hit = c.get(&key).expect("must hit");
        assert_eq!(hit.variants.len(), 1);
        assert_eq!(c.hits(), 1);
        assert_eq!(c.len(), 1);
    }

    #[test]
    fn a_disabled_cache_never_hits() {
        // The --no-piece-cache switch must remove the optimisation entirely,
        // not merely report it differently.
        let c = PieceCache::new(false);
        let key = (1u64, 0u32, 10u64, 20u64);
        c.insert(key, DiffResult::default());
        assert!(c.get(&key).is_none());
        assert_eq!(c.hits(), 0);
        assert!(c.is_empty());
    }

    #[test]
    fn piece_variant_flips_member_coords_with_orientation() {
        let piece = DiffPiece {
            gcontig: 0,
            gs: 100,
            ge: 200,
            pcontig: 3,
            ps: 5000,
            pe: 5100,
            rel: false,
            alt_hash: 0,
            flags: 0,
        };
        let rv = RawVariant {
            kind: VariantKind::Snp,
            start: 10,
            end: 11,
            ref_allele: b"A".to_vec(),
            alt_allele: b"T".to_vec(),
            flags: 0,
        };
        let v = piece_variant(&piece, &rv);
        assert_eq!(v.pivot_start, 5010);
        assert_eq!(v.pivot_end, 5011);
        assert_eq!(v.genome_start, 110);
        assert_eq!(v.strand, '+');

        let rev = DiffPiece {
            rel: true,
            ..piece.clone()
        };
        let vr = piece_variant(&rev, &rv);
        // alt offset 10..11 sits at member 200-11 .. 200-10.
        assert_eq!(vr.genome_start, 189);
        assert_eq!(vr.genome_end, 190);
        assert_eq!(vr.strand, '-');
        assert_eq!(vr.pivot_start, 5010, "pivot coordinates never flip");
    }

    #[test]
    fn alt_is_reverse_complemented_only_when_rel() {
        let piece = DiffPiece {
            gcontig: 0,
            gs: 0,
            ge: 4,
            pcontig: 0,
            ps: 0,
            pe: 4,
            rel: false,
            alt_hash: 0,
            flags: 0,
        };
        assert_eq!(piece.alt(b"ACGT"), b"ACGT".to_vec());
        let rev = DiffPiece { rel: true, ..piece };
        assert_eq!(rev.alt(b"ACGT"), b"ACGT".to_vec(), "revcomp(ACGT) == ACGT");
        let rev2 = DiffPiece {
            rel: true,
            ..DiffPiece {
                gs: 0,
                ge: 4,
                ps: 0,
                pe: 4,
                gcontig: 0,
                pcontig: 0,
                rel: false,
                alt_hash: 0,
                flags: 0,
            }
        };
        assert_eq!(rev2.alt(b"AAGG"), b"CCTT".to_vec());
    }
}
