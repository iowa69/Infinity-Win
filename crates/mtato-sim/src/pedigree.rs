//! Growing one cluster: its root, its pedigree and every leaf genome.
//!
//! The walk is a depth-first traversal from the cluster root. Each node carries
//! a [`NodeState`] - a spine plus the substitutions and events accumulated on the
//! path - and a branch turns its parent's state into its own by drawing its
//! events and applying them. Leaves are rendered to bases and handed to a sink as
//! they are produced, so the peak memory of a run is one genome per cluster, not
//! one per leaf.

use crate::positions::FreeSpace;
use crate::rng::{stream, Rng, TAG_EVENTS, TAG_REFERENCE, TAG_ROOT, TAG_TREE};
use crate::spine::{merge_subs, render, render_ref_range, Spine};
use crate::treegen::Pedigree;
use crate::{config::SimConfig, truth::Mutation, MutationKind};

/// What a cluster's simulation shares between its branches.
struct ClusterCtx<'a> {
    cfg: &'a SimConfig,
    ref_seq: &'a [u8],
    /// Substitutions that took the reference to this cluster's root. They are
    /// shared by every genome of the cluster, so they are kept once here instead
    /// of in each node's state.
    root_subs: Vec<(usize, u8)>,
    /// GC fraction of invented sequence.
    novel_gc: f64,
    space: FreeSpace,
    rng: Rng,
}

/// One node of the pedigree while the walk is running.
#[derive(Clone)]
struct NodeState {
    spine: Spine,
    /// Substitutions added on the branches from the cluster root down to here,
    /// sorted by position.
    subs: Vec<(usize, u8)>,
    /// Every event on that same path, in the order it will be written.
    mutations: Vec<Mutation>,
}

/// The result of simulating one cluster.
pub(crate) struct ClusterBuild {
    /// The pedigree, with branch lengths already equal to the SNP counts applied.
    pub pedigree: Pedigree,
    /// The cluster root's bases.
    pub root_seq: Vec<u8>,
    /// The substitutions that separate the root from the reference.
    pub root_mutations: Vec<Mutation>,
}

/// Simulate one cluster and hand every leaf genome to `on_leaf` as it is built.
///
/// The callback receives the leaf's index in `leaf_names`, its name, its rendered
/// sequence before fragmentation and its event log.
pub(crate) fn build_cluster<F>(
    cluster: u32,
    cfg: &SimConfig,
    ref_seq: &[u8],
    leaf_names: &[String],
    novel_gc: f64,
    on_leaf: &mut F,
) -> Result<ClusterBuild, String>
where
    F: FnMut(usize, &str, Vec<u8>, Vec<Mutation>) -> Result<(), String>,
{
    let mut space = FreeSpace::new(ref_seq.len());
    let mut rng_root = stream(cfg.seed, TAG_ROOT, u64::from(cluster));
    let (root_subs, root_mutations) = diverge(&mut rng_root, cfg, ref_seq, &mut space)?;
    let root_seq = render(&Spine::whole(ref_seq.len()), ref_seq, &root_subs, &[]);

    let mut rng_tree = stream(cfg.seed, TAG_TREE, u64::from(cluster));
    let pedigree = Pedigree::random(&mut rng_tree, leaf_names.len(), cfg.mean_snps_per_branch);

    let mut ctx = ClusterCtx {
        cfg,
        ref_seq,
        root_subs,
        novel_gc,
        space,
        rng: stream(cfg.seed, TAG_EVENTS, u64::from(cluster)),
    };

    let root_state = NodeState {
        spine: Spine::whole(ref_seq.len()),
        subs: Vec::new(),
        mutations: Vec::new(),
    };
    let mut stack: Vec<(usize, NodeState)> = vec![(pedigree.root, root_state)];
    while let Some((id, state)) = stack.pop() {
        if let Some(leaf) = pedigree.nodes[id].leaf {
            let seq = render(&state.spine, ctx.ref_seq, &ctx.root_subs, &state.subs);
            on_leaf(leaf, &leaf_names[leaf], seq, state.mutations)?;
            continue;
        }
        // Children are built left to right so that the generator is consumed in
        // a fixed order, then pushed reversed so the left one is popped first.
        let mut built = Vec::with_capacity(pedigree.nodes[id].children.len());
        for &child in &pedigree.nodes[id].children {
            let mut next = state.clone();
            let n_snps = pedigree.nodes[child].snps as usize;
            apply_branch(&mut ctx, n_snps, &mut next)?;
            built.push((child, next));
        }
        for item in built.into_iter().rev() {
            stack.push(item);
        }
    }

    Ok(ClusterBuild {
        pedigree,
        root_seq,
        root_mutations,
    })
}

/// The substitutions that separate a cluster root from the reference, and the
/// mutation log that records them.
type Divergence = (Vec<(usize, u8)>, Vec<Mutation>);

/// The substitutions that take the reference to a cluster root.
///
/// The count is `divergence` applied to the whole reference length, capped at
/// the number of positions that can carry a substitution at all - an N cannot
/// become a SNP, so counting one would inflate the requested divergence.
fn diverge(
    rng: &mut Rng,
    cfg: &SimConfig,
    ref_seq: &[u8],
    space: &mut FreeSpace,
) -> Result<Divergence, String> {
    let wanted = (cfg.divergence * ref_seq.len() as f64).round() as usize;
    let eligible = ref_seq
        .iter()
        .filter(|&&b| mtato_core::seq::code(b) < 4)
        .count();
    let count = wanted.min(eligible);
    // `block_events = false`: root substitutions must not stop a later branch
    // from deleting the segment they sit in, or segment loss would be impossible
    // at any realistic divergence. They still go into `snp_used`, so no branch
    // SNP can land on one.
    let positions = space.take_positions(rng, ref_seq, count, true, false)?;
    let mut subs = Vec::with_capacity(positions.len());
    let mut mutations = Vec::with_capacity(positions.len());
    for &p in &positions {
        let from = ref_seq[p];
        let to = alt_base(rng, from);
        subs.push((p, to));
        mutations.push(Mutation::snp(p, from, to));
    }
    Ok((subs, mutations))
}

/// Draw and apply everything that happens on one branch.
///
/// The counts come first, then the events, in a fixed order: losses, small
/// deletions, gains, small insertions, SNPs. Deletions go before insertions so
/// that a window is never rejected because of a segment this same branch just
/// added, and SNPs go last so that they can only land on sequence that survived.
fn apply_branch(ctx: &mut ClusterCtx, n_snps: usize, state: &mut NodeState) -> Result<(), String> {
    let cfg = ctx.cfg;
    let n_indels = ctx.rng.poisson(cfg.mean_indels_per_branch);
    let n_gains = ctx.rng.poisson(cfg.segment_events.mean_gains_per_branch);
    let n_losses = ctx
        .rng
        .poisson(ctx.cfg.segment_events.mean_losses_per_branch);

    let loss_lens = draw_lengths(
        &mut ctx.rng,
        n_losses,
        cfg.segment_events.min_len,
        cfg.segment_events.max_len,
    );
    let gain_lens = draw_lengths(
        &mut ctx.rng,
        n_gains,
        cfg.segment_events.min_len,
        cfg.segment_events.max_len,
    );
    let indels: Vec<(bool, usize)> = (0..n_indels)
        .map(|_| {
            let insertion = ctx.rng.coin();
            let len = ctx
                .rng
                .range_inclusive(cfg.indel_min_len, cfg.indel_max_len);
            (insertion, len)
        })
        .collect();

    for &len in &loss_lens {
        delete_event(ctx, state, len, MutationKind::Loss);
    }
    for &(insertion, len) in &indels {
        if !insertion {
            delete_event(ctx, state, len, MutationKind::Deletion);
        }
    }
    for &len in &gain_lens {
        insert_event(ctx, state, len, MutationKind::Gain);
    }
    for &(insertion, len) in &indels {
        if insertion {
            insert_event(ctx, state, len, MutationKind::Insertion);
        }
    }

    if n_snps > 0 {
        let positions = ctx
            .space
            .take_positions(&mut ctx.rng, ctx.ref_seq, n_snps, true, true)?;
        let mut fresh = Vec::with_capacity(positions.len());
        for &p in &positions {
            let from = ctx.ref_seq[p];
            let to = alt_base(&mut ctx.rng, from);
            fresh.push((p, to));
            state.mutations.push(Mutation::snp(p, from, to));
        }
        state.subs = merge_subs(&state.subs, &fresh);
    }

    state.mutations.sort_by_key(|m| (m.pos, kind_rank(m.kind)));
    Ok(())
}

/// One deletion or segment loss: find a window nothing else has claimed, remove
/// it from the spine and record exactly what went.
fn delete_event(ctx: &mut ClusterCtx, state: &mut NodeState, len: usize, kind: MutationKind) {
    // A window may be unobtainable in a small or heavily edited reference. The
    // event is then dropped rather than allowed to overlap another one, because
    // an overlap is what would make the truth pairwise counts wrong.
    let Some(start) = ctx.space.take_window(&mut ctx.rng, len) else {
        return;
    };
    let removed = render_ref_range(ctx.ref_seq, start, start + len, &ctx.root_subs, &state.subs);
    state.spine.delete(start, start + len);
    state.mutations.push(Mutation {
        kind,
        pos: start,
        ref_bases: removed,
        alt_bases: Vec::new(),
    });
}

/// One insertion or segment gain: invent a segment and anchor it in front of a
/// reference position nothing else has claimed.
fn insert_event(ctx: &mut ClusterCtx, state: &mut NodeState, len: usize, kind: MutationKind) {
    let Ok(positions) = ctx
        .space
        .take_positions(&mut ctx.rng, ctx.ref_seq, 1, false, true)
    else {
        return;
    };
    let Some(&pos) = positions.first() else {
        return;
    };
    let seq = random_bases(&mut ctx.rng, len, ctx.novel_gc);
    state.spine.insert_before(pos, seq.clone());
    state.mutations.push(Mutation {
        kind,
        pos,
        ref_bases: Vec::new(),
        alt_bases: seq,
    });
}

fn draw_lengths(rng: &mut Rng, count: usize, min: usize, max: usize) -> Vec<usize> {
    (0..count).map(|_| rng.range_inclusive(min, max)).collect()
}

fn kind_rank(kind: MutationKind) -> u8 {
    match kind {
        MutationKind::Snp => 0,
        MutationKind::Insertion => 1,
        MutationKind::Deletion => 2,
        MutationKind::Gain => 3,
        MutationKind::Loss => 4,
    }
}

/// A base that is not `from`, drawn uniformly from the other three.
fn alt_base(rng: &mut Rng, from: u8) -> u8 {
    const BASES: [u8; 4] = *b"ACGT";
    let others: [u8; 3] = {
        let mut out = [b'A'; 3];
        let mut i = 0;
        for b in BASES {
            if b != from {
                out[i] = b;
                i += 1;
            }
        }
        // `from` is a normalised ACGT base whenever this is called for a SNP;
        // if it is not, keep the three first bases rather than leave a hole.
        if i < 3 {
            out[2] = BASES[3];
        }
        out
    };
    others[rng.below_usize(3)]
}

/// `len` random bases at the requested GC fraction.
///
/// The GC decision is made on integer thresholds rather than a float comparison,
/// so the same seed gives the same bases on every platform.
fn random_bases(rng: &mut Rng, len: usize, gc: f64) -> Vec<u8> {
    let threshold = (gc.clamp(0.0, 1.0) * (u64::from(u32::MAX) + 1) as f64) as u64;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        let w = rng.next_u64();
        let gc_base = (w >> 32) < threshold;
        let first = (w >> 31) & 1 == 0;
        out.push(match (gc_base, first) {
            (true, true) => b'G',
            (true, false) => b'C',
            (false, true) => b'A',
            (false, false) => b'T',
        });
    }
    out
}

/// Draw the random reference sequence of `RefSpec::Random`.
pub(crate) fn random_reference(seed: u64, length: usize, gc: f64) -> Vec<u8> {
    let mut rng = stream(seed, TAG_REFERENCE, 0);
    random_bases(&mut rng, length, gc)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alt_base_never_returns_the_original() {
        let mut rng = Rng::from_seed(1);
        for &b in b"ACGT" {
            for _ in 0..200 {
                let a = alt_base(&mut rng, b);
                assert_ne!(a, b);
                assert!(matches!(a, b'A' | b'C' | b'G' | b'T'));
            }
        }
    }

    #[test]
    fn alt_base_is_uniform_over_the_three_others() {
        let mut rng = Rng::from_seed(2);
        let mut counts = [0usize; 256];
        for _ in 0..30_000 {
            counts[alt_base(&mut rng, b'A') as usize] += 1;
        }
        assert_eq!(counts[b'A' as usize], 0);
        for &b in b"CGT" {
            let c = counts[b as usize];
            assert!((9_000..11_000).contains(&c), "{b} got {c}");
        }
    }

    #[test]
    fn random_bases_hits_the_requested_gc() {
        let mut rng = Rng::from_seed(3);
        for &gc in &[0.0f64, 0.25, 0.5, 0.75, 1.0] {
            let seq = random_bases(&mut rng, 40_000, gc);
            let got =
                seq.iter().filter(|&&b| b == b'G' || b == b'C').count() as f64 / seq.len() as f64;
            assert!((got - gc).abs() < 0.02, "gc {gc}: got {got}");
            assert!(seq.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')));
        }
        assert!(random_bases(&mut rng, 0, 0.5).is_empty());
    }

    #[test]
    fn random_bases_is_deterministic_and_seed_sensitive() {
        let a = random_bases(&mut Rng::from_seed(4), 1000, 0.5);
        let b = random_bases(&mut Rng::from_seed(4), 1000, 0.5);
        let c = random_bases(&mut Rng::from_seed(5), 1000, 0.5);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn diverge_reaches_the_requested_count() {
        let cfg = SimConfig {
            divergence: 0.05,
            ..SimConfig::default()
        };
        let ref_seq = random_reference(1, 10_000, 0.5);
        let mut rng = Rng::from_seed(9);
        let mut space = FreeSpace::new(ref_seq.len());
        let (subs, muts) = diverge(&mut rng, &cfg, &ref_seq, &mut space).unwrap();
        assert_eq!(subs.len(), 500);
        assert_eq!(muts.len(), 500);
        assert!(muts.iter().all(|m| m.kind == MutationKind::Snp));
        assert!(subs.windows(2).all(|w| w[0].0 < w[1].0));
        for (m, &(p, to)) in muts.iter().zip(subs.iter()) {
            assert_eq!(m.pos, p);
            assert_eq!(m.ref_bases, vec![ref_seq[p]]);
            assert_eq!(m.alt_bases, vec![to]);
        }
    }

    #[test]
    fn diverge_skips_n_positions() {
        let cfg = SimConfig {
            divergence: 0.5,
            ..SimConfig::default()
        };
        let mut ref_seq = vec![b'A'; 1000];
        ref_seq.iter_mut().take(900).for_each(|b| *b = b'N');
        let mut rng = Rng::from_seed(10);
        let mut space = FreeSpace::new(ref_seq.len());
        let (subs, _) = diverge(&mut rng, &cfg, &ref_seq, &mut space).unwrap();
        // Only 100 positions can carry a substitution, so the request saturates.
        assert_eq!(subs.len(), 100);
        assert!(subs.iter().all(|&(p, _)| p >= 900));
    }

    #[test]
    fn diverge_zero_changes_nothing() {
        let cfg = SimConfig {
            divergence: 0.0,
            ..SimConfig::default()
        };
        let ref_seq = random_reference(2, 500, 0.5);
        let mut space = FreeSpace::new(ref_seq.len());
        let (subs, muts) = diverge(&mut Rng::from_seed(11), &cfg, &ref_seq, &mut space).unwrap();
        assert!(subs.is_empty());
        assert!(muts.is_empty());
    }
}
