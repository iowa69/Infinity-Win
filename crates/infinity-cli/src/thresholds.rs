//! `infinity thresholds`: help choose the cluster distance for a data set.
//!
//! The cluster distance decides which genomes are compared base by base. Only
//! genomes in the same cluster get exact SNP distances to each other, so the
//! threshold must be wide enough to put every pair you care about in one
//! cluster, and narrow enough that a cluster does not mix lineages (which makes
//! the comparison slower and shrinks the core genome every member shares).
//!
//! This command sketches the genomes (seconds, the same Focus 0 stage a run
//! starts with), then clusters them exactly as a run would at a range of
//! thresholds and reports, for each one, how many clusters it makes and how
//! many genome pairs sit close to the cut. A threshold with few pairs near it
//! falls in a natural gap of the data: small changes to it, or to the sketch
//! estimates, do not move genomes between clusters.

use std::path::PathBuf;

use clap::Parser;
use rayon::prelude::*;

use infinity_core::cluster::cluster_genomes;
use infinity_core::engine::{run_collecting, RunConfig};
use infinity_core::params::{Parameters, Resolution};
use infinity_core::sketch;

#[derive(Parser, Debug)]
#[command(about = "Suggest a cluster distance for a set of genomes")]
pub struct ThresholdArgs {
    /// Genome files (`.fa`, `.fasta`, `.fna`, optionally `.gz`) or folders of them.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Worker threads. 0 uses every core.
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// The largest SNP distance you need exactly (for example your outbreak
    /// threshold, or 10x it to see the context). Used to warn when a cluster
    /// distance is too tight for it.
    #[arg(long, default_value_t = 50)]
    max_snps: u64,
}

/// Thresholds tried, as sketch distances (1 - identity).
const CANDIDATES: [f64; 9] = [
    0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.015, 0.02,
];

/// Below this the sketch distance is too noisy to cluster on reliably (with
/// the default sketch, about one hash per 2 kb, the estimate's spread at
/// 0.0005 is about a quarter of its value).
const NOISE_FLOOR: f64 = 0.001;

/// Above this the block shortcut and the anchors stop working well (PROJECT.md
/// section 15): comparisons get slower and less exact.
const DIVERGENCE_CEILING: f64 = 0.01;

pub fn cmd_thresholds(args: ThresholdArgs) -> i32 {
    let params = Parameters {
        threads: args.threads,
        resolution: Resolution::Sketch,
        refine_backbone: false,
        ..Default::default()
    };
    eprintln!("sketching the genomes ...");
    let (results, _) = match run_collecting(RunConfig {
        inputs: args.inputs.clone(),
        params: params.clone(),
        outdir: None,
    }) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("infinity thresholds: {e}");
            return 1;
        }
    };
    let recs = &results.records;
    let n = recs.len();
    if n < 2 {
        eprintln!("infinity thresholds: need at least two genomes, got {n}");
        return 1;
    }
    let mut lens: Vec<u64> = recs.iter().map(|r| r.stats.length).collect();
    lens.sort_unstable();
    let median_len = lens[n / 2] as f64;

    // All pairwise sketch distances.
    let k = results.params.k;
    let dists: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            let recs = &recs;
            ((i + 1)..n).map(move |j| sketch::distance(&recs[i].sketch, &recs[j].sketch, k))
        })
        .collect();
    let mut nearest = vec![f64::INFINITY; n];
    {
        let mut idx = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = dists[idx];
                idx += 1;
                nearest[i] = nearest[i].min(d);
                nearest[j] = nearest[j].min(d);
            }
        }
    }

    println!(
        "{n} genomes, median length {:.2} Mb. A sketch distance d means about \
         d x {:.0} differing bases between two genomes of that size.",
        median_len / 1e6,
        median_len
    );
    println!();
    println!(
        "{:>9}  {:>9}  {:>11}  {:>8}  {:>10}  {:>7}  {:>11}  note",
        "distance", "identity", "~SNPs apart", "clusters", "singletons", "largest", "near the cut"
    );
    let mut rows = Vec::new();
    for &t in &CANDIDATES {
        let mut p = results.params.clone();
        p.cluster_dist = t;
        let c = cluster_genomes(recs, &p);
        let singletons = c.clusters.iter().filter(|x| x.size() == 1).count();
        let largest = c.clusters.iter().map(|x| x.size()).max().unwrap_or(0);
        // Pairs whose distance is within a factor 1.5 of the threshold: these
        // are the ones a slightly different threshold would split or join.
        let near = dists
            .iter()
            .filter(|&&d| d >= t / 1.5 && d <= t * 1.5)
            .count();
        let snps = t * median_len;
        let mut note = Vec::new();
        if t < NOISE_FLOOR {
            note.push("below the sketch's precision");
        }
        if t > DIVERGENCE_CEILING {
            note.push("mixes lineages: slower, smaller core");
        }
        if snps < args.max_snps as f64 * 3.0 {
            note.push("too tight for --max-snps");
        }
        println!(
            "{:>9}  {:>8.2}%  {:>11.0}  {:>8}  {:>10}  {:>7}  {:>11}  {}",
            t,
            (1.0 - t) * 100.0,
            snps,
            c.clusters.len(),
            singletons,
            largest,
            near,
            note.join("; ")
        );
        rows.push((t, near, c.clusters.len(), snps));
    }

    // Recommendation: among thresholds that are precise enough, wide enough for
    // --max-snps and not past the divergence ceiling, the one with the fewest
    // pairs near the cut; ties go to the one closest to the default (0.005), so
    // a data set without any gap keeps the tested default.
    let usable: Vec<&(f64, usize, usize, f64)> = rows
        .iter()
        .filter(|r| {
            r.0 >= NOISE_FLOOR && r.0 <= DIVERGENCE_CEILING && r.3 >= args.max_snps as f64 * 3.0
        })
        .collect();
    println!();
    match usable.iter().min_by(|a, b| {
        let off = |t: f64| (t / Parameters::default().cluster_dist).ln().abs();
        a.1.cmp(&b.1).then(off(a.0).total_cmp(&off(b.0)))
    }) {
        Some(best) => {
            println!(
                "Suggested cluster distance: {} ({:.2}% identity, {} clusters).",
                best.0,
                (1.0 - best.0) * 100.0,
                best.2
            );
            println!(
                "It sits in the quietest gap of your data: {} genome pair(s) lie near it, \
                 so the clusters do not depend on small changes.",
                best.1
            );
            println!("Use it with:  infinity run ... --cluster-dist {}", best.0);
        }
        None => println!(
            "No candidate fits: --max-snps {} needs a distance above {:.4}, which is past \
             the point where clusters stay exact. Lower --max-snps or accept slower, \
             less exact clusters.",
            args.max_snps,
            args.max_snps as f64 * 3.0 / median_len
        ),
    }

    let mut nn: Vec<f64> = nearest.into_iter().filter(|d| d.is_finite()).collect();
    nn.sort_by(|a, b| a.total_cmp(b));
    let q = |p: f64| nn[((nn.len() - 1) as f64 * p).round() as usize];
    println!();
    println!(
        "Nearest neighbour of each genome (sketch distance): median {:.5} (~{:.0} SNPs), \
         90th percentile {:.5} (~{:.0} SNPs).",
        q(0.5),
        q(0.5) * median_len,
        q(0.9),
        q(0.9) * median_len
    );
    println!(
        "Sketch distances are estimates: below about {NOISE_FLOOR} they cannot separate genomes \
         reliably. The exact SNP distances come from the run itself."
    );
    0
}
