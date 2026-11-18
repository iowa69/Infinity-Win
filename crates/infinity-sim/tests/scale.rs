//! Scale checks, ignored by default.
//!
//! These write real datasets at the sizes section 10.3 and 10.4 ask for, so they
//! produce hundreds of megabytes of FASTA. Run them deliberately:
//!
//! ```text
//! cargo test -p infinity-sim --test scale -- --ignored --nocapture --test-threads=2
//! ```
//!
//! What they pin down is not just "it finishes": every one of them reads the
//! truth back off the disk and requires it to equal what `simulate` returned, so
//! a dataset that is too large for the writers to keep consistent does not slip
//! through.

use std::time::Instant;

use infinity_sim::{simulate, RefSpec, SegmentEvents, SimConfig, Truth};

fn run(tag: &str, cfg: SimConfig) {
    let dir = std::env::temp_dir().join(format!("infinity-scale-{tag}"));
    let _ = std::fs::remove_dir_all(&dir);
    let t = Instant::now();
    let out = simulate(&cfg, &dir).unwrap();
    let sim = t.elapsed();

    let t = Instant::now();
    let back = Truth::read_from(&dir).unwrap();
    let read = t.elapsed();
    assert_eq!(
        back, out.truth,
        "{tag}: the truth on disk is not the truth returned"
    );

    let bytes: u64 = out
        .genomes
        .iter()
        .map(|p| std::fs::metadata(p).unwrap().len())
        .sum();
    println!(
        "{tag}: {} genomes, {:.1} MB of FASTA, simulate {sim:?}, read_from {read:?}",
        out.genomes.len(),
        bytes as f64 / 1e6
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// T1 of section 10.3: a 2 Mb random reference, one cluster, 50 leaves, three
/// SNPs per branch.
#[test]
#[ignore = "writes about 100 MB of FASTA"]
fn t1_shape() {
    run(
        "T1",
        SimConfig {
            seed: 1,
            reference: RefSpec::Random {
                length: 2_000_000,
                gc: 0.5,
            },
            clusters: 1,
            leaves_per_cluster: 50,
            mean_snps_per_branch: 3.0,
            ..SimConfig::default()
        },
    );
}

/// T3 of section 10.3: T2 with every genome cut into 100 contigs, shuffled and
/// flipped. Scaled to 20 leaves so the check stays quick; the shape that matters
/// is the contig count.
#[test]
#[ignore = "writes about 40 MB of FASTA"]
fn t3_shape() {
    run(
        "T3",
        SimConfig {
            seed: 3,
            reference: RefSpec::Random {
                length: 2_000_000,
                gc: 0.5,
            },
            clusters: 1,
            leaves_per_cluster: 20,
            mean_snps_per_branch: 3.0,
            mean_indels_per_branch: 0.5,
            contigs: 100,
            random_orientation: true,
            ..SimConfig::default()
        },
    );
}

/// T4 of section 10.3: T2 plus segment gain and loss.
#[test]
#[ignore = "writes about 40 MB of FASTA"]
fn t4_shape() {
    run(
        "T4",
        SimConfig {
            seed: 4,
            reference: RefSpec::Random {
                length: 2_000_000,
                gc: 0.5,
            },
            clusters: 1,
            leaves_per_cluster: 20,
            mean_snps_per_branch: 3.0,
            mean_indels_per_branch: 0.5,
            segment_events: SegmentEvents {
                mean_gains_per_branch: 0.3,
                mean_losses_per_branch: 0.3,
                min_len: 500,
                max_len: 5000,
            },
            ..SimConfig::default()
        },
    );
}

/// T5 of section 10.3: five clusters of twenty, 2% apart.
#[test]
#[ignore = "writes about 200 MB of FASTA"]
fn t5_shape() {
    run(
        "T5",
        SimConfig {
            seed: 5,
            reference: RefSpec::Random {
                length: 2_000_000,
                gc: 0.5,
            },
            clusters: 5,
            leaves_per_cluster: 20,
            mean_snps_per_branch: 3.0,
            mean_indels_per_branch: 0.5,
            ..SimConfig::default()
        },
    );
}

/// A small reference at 30% divergence with heavy segment traffic. This is the
/// configuration that drives the occupancy sampler off rejection sampling and
/// onto exact enumeration, and the one most likely to run out of room, so it is
/// the interesting edge of the design rather than a large run.
#[test]
#[ignore = "exercises the sampler fallback; writes about 2 MB"]
fn dense_events_on_a_small_reference() {
    run(
        "stress",
        SimConfig {
            seed: 7,
            reference: RefSpec::Random {
                length: 60_000,
                gc: 0.5,
            },
            clusters: 2,
            leaves_per_cluster: 15,
            divergence: 0.30,
            mean_snps_per_branch: 40.0,
            mean_indels_per_branch: 3.0,
            segment_events: SegmentEvents {
                mean_gains_per_branch: 1.0,
                mean_losses_per_branch: 1.0,
                min_len: 500,
                max_len: 2000,
            },
            contigs: 12,
            random_orientation: true,
            ..SimConfig::default()
        },
    );
}

/// Section 10.4's dataset shape - 5 Mb genomes in 20 clusters - at a tenth of the
/// genome count. The full 1,000 genomes take about 20 s and 5 GB of FASTA on a
/// laptop; raise `leaves_per_cluster` to 50 to reproduce that, and check first
/// that the output directory has the room.
#[test]
#[ignore = "writes about 500 MB of FASTA"]
fn performance_dataset_shape() {
    run(
        "perf",
        SimConfig {
            seed: 6,
            reference: RefSpec::Random {
                length: 5_000_000,
                gc: 0.5,
            },
            clusters: 20,
            leaves_per_cluster: 5,
            mean_snps_per_branch: 3.0,
            ..SimConfig::default()
        },
    );
}
