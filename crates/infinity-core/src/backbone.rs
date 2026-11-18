//! Exact distances between clusters.
//!
//! The backbone of section 6.4 is built from sketch distances, so the branches
//! between clusters stay approximate even when every cluster is exact. Once the
//! clusters are done, this module compares the pivots to each other with the
//! same engine (blocks, anchors, base-level resolution and the SNP rules) and
//! rebuilds the backbone from SNP counts. There are few pivots, so this costs
//! little next to the clusters themselves.
//!
//! Pivots further apart than [`MAX_EXACT_SKETCH_DIST`] are not aligned: the
//! block shortcut and the anchors stop working at that divergence (section 15),
//! and an SNP count between different lineages is not what a user of this tool
//! needs. Those pairs keep a sketch estimate, converted to SNPs, and are listed
//! in [`ExactBackbone::estimated`].

use rayon::prelude::*;
use rayon::ThreadPool;

use crate::compare::{compare_member, PieceCache};
use crate::engine::RunResults;
use crate::filter;
use crate::params::Parameters;
use crate::pivot::PivotIndex;
use crate::profile::{concat_pivot, contig_offsets, pair_distance, Profile};
use crate::sketch;
use crate::tree::{midpoint_root, neighbor_joining, Tree};

/// Above this many clusters the exact backbone is skipped (it is quadratic).
pub const MAX_PIVOTS: usize = 200;

/// Sketch distance above which two pivots are not aligned (99% identity).
/// Further apart, clusters are different lineages: aligning them costs far
/// more than all the clusters together (measured: 105 s of a 116 s run on 137
/// *S. aureus*), and an exact count of thousands of SNPs between lineages is
/// not what this tool is for.
pub const MAX_EXACT_SKETCH_DIST: f64 = 0.01;

/// The backbone rebuilt from exact pivot-to-pivot comparisons.
#[derive(Debug, Clone, Default)]
pub struct ExactBackbone {
    /// Pivot genome ids, in matrix order.
    pub pivots: Vec<u32>,
    pub names: Vec<String>,
    /// Core SNPs between pivots (estimated for the pairs in `estimated`).
    pub snps: Vec<Vec<u64>>,
    /// Bases compared per pair (0 for estimated pairs).
    pub comparable: Vec<Vec<u64>>,
    /// Pairs too divergent to align, as (row, column) with row < column.
    pub estimated: Vec<(usize, usize)>,
    /// Neighbour joining on `snps`, midpoint rooted, branch lengths in SNPs.
    pub tree: Tree,
}

/// Compare every pair of pivots exactly. Returns `None` when there are fewer
/// than two pivots or more than [`MAX_PIVOTS`], or when cancelled.
pub(crate) fn refine(
    pool: &ThreadPool,
    results: &RunResults,
    pivot_ids: &[u32],
    params: &Parameters,
    cancelled: &(dyn Fn() -> bool + Sync),
) -> Result<Option<ExactBackbone>, String> {
    let p = pivot_ids.len();
    if !(2..=MAX_PIVOTS).contains(&p) {
        return Ok(None);
    }
    let names: Vec<String> = pivot_ids
        .iter()
        .map(|&id| results.names[id as usize].clone())
        .collect();
    let mut snps = vec![vec![0u64; p]; p];
    let mut comparable = vec![vec![0u64; p]; p];
    let mut estimated = Vec::new();

    for i in 0..p {
        if cancelled() {
            return Ok(None);
        }
        let rec_i = &results.records[pivot_ids[i] as usize];
        // Which partners are close enough to align.
        let mut close = Vec::new();
        for j in (i + 1)..p {
            let rec_j = &results.records[pivot_ids[j] as usize];
            let d = sketch::distance(&rec_i.sketch, &rec_j.sketch, params.k);
            if d <= MAX_EXACT_SKETCH_DIST {
                close.push(j);
            } else {
                let len = rec_i.stats.length.min(rec_j.stats.length) as f64;
                let est = (d * len).round() as u64;
                snps[i][j] = est;
                snps[j][i] = est;
                estimated.push((i, j));
            }
        }
        if close.is_empty() {
            continue;
        }
        let data_i = crate::fasta::read_genome(&rec_i.path)
            .map_err(|e| format!("reading pivot {}: {e}", rec_i.path.display()))?;
        let idx = PivotIndex::new(&rec_i.name, &data_i, params);
        drop(data_i);
        let offsets = contig_offsets(&idx);
        let concat = concat_pivot(&idx);
        let contig_names: Vec<String> = idx.contigs.iter().map(|c| c.name.clone()).collect();
        let cache = PieceCache::new(false);
        let row: Vec<(usize, u64, u64)> = pool.install(|| {
            close
                .par_iter()
                .filter_map(|&j| {
                    if cancelled() {
                        return None;
                    }
                    let rec_j = &results.records[pivot_ids[j] as usize];
                    let data_j = crate::fasta::read_genome(&rec_j.path).ok()?;
                    let cmp = compare_member(&idx, &data_j, pivot_ids[j], params, &cache, true);
                    let prof = Profile::from_comparison(&cmp, &offsets);
                    let fo = filter::apply(
                        &[Profile::pivot(pivot_ids[i]), prof],
                        &contig_names,
                        &idx.contig_lens,
                        &offsets,
                        idx.total_len,
                        &params.snp_filter,
                    );
                    let d = pair_distance(&fo.profiles[0], &fo.profiles[1], idx.total_len, &concat);
                    Some((j, d.snps, d.comparable))
                })
                .collect()
        });
        for (j, s, c) in row {
            snps[i][j] = s;
            snps[j][i] = s;
            comparable[i][j] = c;
            comparable[j][i] = c;
        }
    }
    if cancelled() {
        return Ok(None);
    }
    let dist: Vec<Vec<f64>> = snps
        .iter()
        .map(|r| r.iter().map(|&v| v as f64).collect())
        .collect();
    let mut t = neighbor_joining(&names, &dist);
    t.clamp_negative_lengths();
    Ok(Some(ExactBackbone {
        pivots: pivot_ids.to_vec(),
        names,
        snps,
        comparable,
        estimated,
        tree: midpoint_root(t),
    }))
}

/// Scale every branch of a tree by `factor` (sketch distances to SNPs).
pub fn scaled(t: &Tree, factor: f64) -> Tree {
    match t {
        Tree::Leaf { name, len } => Tree::Leaf {
            name: name.clone(),
            len: len * factor,
        },
        Tree::Internal { len, children } => Tree::Internal {
            len: len * factor,
            children: children.iter().map(|c| scaled(c, factor)).collect(),
        },
    }
}
