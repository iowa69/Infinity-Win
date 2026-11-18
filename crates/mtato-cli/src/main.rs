//! `mtato`, the command-line consumer of the engine (section 9).
//!
//! The engine is a library and knows nothing about this binary (section 0 rule
//! 2); everything shown here comes off the same event stream the GUI reads.

mod bench;

use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use clap::{Parser, Subcommand, ValueEnum};

use mtato_core::engine::{run, RunConfig};
use mtato_core::events::{event_bus, Command, Envelope, Event};
use mtato_core::params::Parameters;
use mtato_core::seq::LEVELS;

#[derive(Parser, Debug)]
#[command(
    name = "mtato",
    version = mtato_core::VERSION,
    about = "Progressive, cluster-first, position-aware comparison of nearly identical genomes",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    command: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Compare a set of assembled genomes.
    Run(RunArgs),
    /// Compare an output folder with a truth folder and report the measurements.
    Bench(bench::BenchArgs),
    /// Simulate genomes with known truth.
    Simulate(SimulateArgs),
}

/// How much the run narrates while it works.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Verbosity {
    /// Only errors and the final summary.
    Quiet,
    /// Stage progress, one line per stage (the default).
    Normal,
    /// Every event, as it happens.
    Verbose,
}

#[derive(Parser, Debug)]
#[command(about = "Compare a set of assembled genomes")]
struct RunArgs {
    /// Genome files (`.fa`, `.fasta`, `.fna`, optionally `.gz`) or folders of them.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Output folder. Created if it does not exist.
    #[arg(short, long)]
    output: PathBuf,

    /// Worker threads. 0 uses every core.
    #[arg(long, default_value_t = 0)]
    threads: usize,

    /// Anchor k-mer length. Must be odd and at most 31.
    #[arg(long, default_value_t = 31)]
    k: usize,

    /// Anchor densities at levels 1, 2 and 3, for example `16384,2048,256`.
    #[arg(long, value_parser = parse_densities, default_value = "16384,2048,256")]
    densities: [u64; LEVELS],

    /// Level whose anchors form the sketch (1-3).
    #[arg(long, default_value_t = 2)]
    sketch_level: u8,

    /// Sketch distance threshold for joining a cluster.
    #[arg(long, default_value_t = 0.005)]
    cluster_dist: f64,

    /// Collinearity limit: maximum span ratio between member and pivot.
    #[arg(long, default_value_t = 10.0)]
    max_span_ratio: f64,

    /// Collinearity limit: slack added to the ratio test, in bases.
    #[arg(long, default_value_t = 50000)]
    max_span_extra: u64,

    /// Longest trimmed piece aligned with dynamic programming.
    #[arg(long, default_value_t = 5000)]
    dp_max_len: usize,

    /// Below this identity an aligned piece becomes one complex variant.
    #[arg(long, default_value_t = 0.70)]
    min_identity: f64,

    /// x-drop stop rule for ungapped extension at contig ends.
    #[arg(long, default_value_t = 15)]
    xdrop: i32,

    // --- experiment switches (section 11) ---
    /// Use leaf blocks only; no nested levels.
    #[arg(long)]
    flat: bool,

    /// Do not hash blocks: compare every stretch between shared anchors.
    #[arg(long)]
    no_block_skip: bool,

    /// Do not reuse a resolved difference piece.
    #[arg(long)]
    no_piece_cache: bool,

    /// Do not handle identical genomes once.
    #[arg(long)]
    no_profile_grouping: bool,

    /// Compare every genome to one pivot instead of clustering first.
    #[arg(long)]
    single_cluster: bool,

    /// Sequences of non-pivot genomes are kept while they fit this many bytes.
    /// 0 means a quarter of physical RAM.
    #[arg(long, default_value_t = 0)]
    memory_budget: u64,

    #[arg(long, value_enum, default_value_t = Verbosity::Normal)]
    verbosity: Verbosity,

    /// Stop after the blurry tree, without refining any cluster.
    #[arg(long)]
    focus0_only: bool,
}

#[derive(Parser, Debug)]
#[command(about = "Simulate genomes with known truth")]
struct SimulateArgs {
    /// Where to write the simulated genomes and the truth files.
    #[arg(short, long)]
    output: PathBuf,

    /// Random seed. The same seed gives byte-identical output.
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Reference FASTA to mutate. Omit to generate a random reference.
    #[arg(long)]
    reference: Option<PathBuf>,

    /// Length of the generated reference, when none is supplied.
    #[arg(long, default_value_t = 2_000_000)]
    length: usize,

    /// GC fraction of the generated reference.
    #[arg(long, default_value_t = 0.5)]
    gc: f64,

    /// Number of clusters.
    #[arg(short = 'm', long, default_value_t = 1)]
    clusters: usize,

    /// Divergence of each cluster root from the reference, as a fraction.
    #[arg(long, default_value_t = 0.02)]
    divergence: f64,

    /// Leaves per cluster.
    #[arg(short = 'n', long, default_value_t = 10)]
    leaves: usize,

    /// Mean SNPs per branch, drawn from a Poisson distribution.
    #[arg(long, default_value_t = 3.0)]
    mean_snps: f64,

    /// Mean small indels per branch. 0 switches indels off.
    #[arg(long, default_value_t = 0.0)]
    mean_indels: f64,

    /// Mean segment gain/loss events per branch. 0 switches them off.
    #[arg(long, default_value_t = 0.0)]
    mean_segments: f64,

    /// Cut each genome into this many contigs, shuffled and flipped. 1 keeps it
    /// unfragmented.
    #[arg(short = 'c', long, default_value_t = 1)]
    contigs: usize,
}

fn parse_densities(s: &str) -> Result<[u64; LEVELS], String> {
    let parts: Vec<&str> = s.split(',').map(str::trim).collect();
    if parts.len() != LEVELS {
        return Err(format!(
            "expected {LEVELS} comma-separated densities, got {}",
            parts.len()
        ));
    }
    let mut out = [0u64; LEVELS];
    for (i, p) in parts.iter().enumerate() {
        out[i] = p
            .parse::<u64>()
            .map_err(|e| format!("density {i} is not a number: {e}"))?;
    }
    Ok(out)
}

fn main() {
    let cli = Cli::parse();
    let code = match cli.command {
        Cmd::Run(args) => cmd_run(args),
        Cmd::Bench(args) => bench::cmd_bench(args),
        Cmd::Simulate(args) => cmd_simulate(args),
    };
    std::process::exit(code);
}

fn cmd_run(args: RunArgs) -> i32 {
    let mut params = Parameters {
        k: args.k,
        densities: args.densities,
        sketch_level: args.sketch_level,
        cluster_dist: args.cluster_dist,
        max_span_ratio: args.max_span_ratio,
        max_span_extra: args.max_span_extra,
        dp_max_len: args.dp_max_len,
        min_identity: args.min_identity,
        xdrop: args.xdrop,
        threads: args.threads,
        block_skip: !args.no_block_skip,
        hierarchical: !args.flat,
        piece_cache: !args.no_piece_cache,
        profile_grouping: !args.no_profile_grouping,
        single_cluster: args.single_cluster,
        memory_budget_bytes: args.memory_budget,
    };
    if let Err(e) = params.validate() {
        eprintln!("mtato: invalid parameters: {e}");
        return 2;
    }
    if args.focus0_only {
        // Stopping at the blurry tree means never refining, which is expressed
        // as a cluster threshold of nothing: every genome is its own cluster.
        // The engine still reports that honestly rather than pretending to be
        // exact.
        params.cluster_dist = f64::MIN_POSITIVE;
    }

    let cfg = RunConfig {
        inputs: args.inputs.clone(),
        params,
        outdir: Some(args.output.clone()),
    };

    let (sink, rx) = event_bus();
    sink.attach_log(&args.output).unwrap_or_else(|e| {
        eprintln!(
            "mtato: warning: cannot write events.jsonl in {}: {e}",
            args.output.display()
        );
    });

    let started = Instant::now();
    let handle = std::thread::spawn(move || {
        let (_tx, cmd_rx) = crossbeam_channel::unbounded::<Command>();
        run(cfg, sink, cmd_rx)
    });

    let verbosity = args.verbosity;
    let stdout = BufWriter::new(std::io::stdout().lock());
    let mut out = stdout;
    let mut last_stage = String::new();
    let mut summary_line: Option<String> = None;
    let mut exit_code = 0;

    while let Ok(env) = rx.recv() {
        if let Some(line) = render(&env, verbosity, &mut last_stage) {
            let _ = writeln!(out, "{line}");
            let _ = out.flush();
        }
        if let Event::Finished(s) = &env.event {
            summary_line = Some(format_summary(s, started.elapsed().as_secs_f64()));
            if s.cancelled {
                exit_code = 130;
            }
        }
    }

    match handle.join() {
        Ok(Ok(_)) => {}
        Ok(Err(e)) => {
            let _ = writeln!(out, "mtato: error: {e}");
            let _ = out.flush();
            return 1;
        }
        Err(_) => {
            let _ = writeln!(out, "mtato: error: the engine thread panicked");
            let _ = out.flush();
            return 1;
        }
    }
    if let Some(s) = summary_line {
        let _ = writeln!(out, "{s}");
    }
    let _ = writeln!(out, "output written to {}", args.output.display());
    let _ = out.flush();
    exit_code
}

/// Turn one event into a line for the terminal, or `None` when this verbosity
/// level shows nothing for it.
fn render(env: &Envelope, verbosity: Verbosity, last_stage: &mut String) -> Option<String> {
    let t = env.t_ms / 1000.0;
    match &env.event {
        Event::Progress {
            stage, done, total, ..
        } => {
            if verbosity == Verbosity::Quiet {
                return None;
            }
            // One line per stage change, plus a completion line, so a long scan
            // does not scroll the terminal away.
            if *stage == *last_stage && !(*total > 0 && done == total) {
                return None;
            }
            *last_stage = stage.clone();
            Some(if *total > 0 {
                format!("[{t:8.2}s] {stage:<9} {done}/{total}")
            } else {
                format!("[{t:8.2}s] {stage}")
            })
        }
        Event::Warning { text } => Some(format!("[{t:8.2}s] warning: {text}")),
        Event::Clustered {
            threshold,
            clusters,
        } => {
            if verbosity == Verbosity::Quiet {
                return None;
            }
            let largest = clusters.iter().map(|c| c.members.len()).max().unwrap_or(0);
            Some(format!(
                "[{t:8.2}s] clustered {} genomes into {} clusters at d<={threshold} (largest {largest})",
                clusters.iter().map(|c| c.members.len()).sum::<usize>(),
                clusters.len()
            ))
        }
        Event::BackboneTree { pivots, .. } => {
            if verbosity == Verbosity::Quiet {
                return None;
            }
            Some(format!(
                "[{t:8.2}s] backbone tree ready over {} clusters  <-- first blurry tree",
                pivots.len()
            ))
        }
        Event::ClusterFocus(info) => {
            if verbosity != Verbosity::Verbose {
                return None;
            }
            let marker = if info.approximate {
                " (approximate)"
            } else {
                ""
            };
            Some(format!(
                "[{t:8.2}s] cluster {} focus {}{}: {} members",
                info.cluster,
                info.focus,
                marker,
                info.members.len()
            ))
        }
        Event::GenomeScanned(g) => {
            if verbosity != Verbosity::Verbose {
                return None;
            }
            Some(format!(
                "[{t:8.2}s] scanned {}: {} bp, {} contigs, N50 {}, sketch {}",
                g.name, g.length, g.contigs, g.n50, g.sketch_size
            ))
        }
        Event::MapUpdate(m) => {
            if verbosity != Verbosity::Verbose {
                return None;
            }
            Some(format!(
                "[{t:8.2}s] map {}: {} pieces, {} variants",
                m.genome_name,
                m.pieces.len(),
                m.variants.len()
            ))
        }
        Event::Counters(c) => {
            if verbosity != Verbosity::Verbose {
                return None;
            }
            Some(format!(
                "[{t:8.2}s] counters: blocks {} ({} matched), pieces {}, cache hits {}, alignments {}",
                c.blocks_looked_up,
                c.total_matched(),
                c.pieces_compared,
                c.cache_hits,
                c.alignments_run
            ))
        }
        Event::Finished(_) => None,
    }
}

fn format_summary(s: &mtato_core::summary::RunSummary, secs: f64) -> String {
    let first = s
        .timings
        .first_tree_ms
        .map(|m| format!("{:.2}s", m as f64 / 1000.0))
        .unwrap_or_else(|| "n/a".to_string());
    let exact = s
        .timings
        .exact_tree_ms
        .map(|m| format!("{:.2}s", m as f64 / 1000.0))
        .unwrap_or_else(|| "not reached".to_string());
    let mem = s
        .peak_memory_bytes
        .map(|b| format!("{:.0} MB", b as f64 / 1e6))
        .unwrap_or_else(|| "unknown".to_string());
    let identity = s
        .counters
        .block_identity()
        .map(|f| format!("{:.2}%", f * 100.0))
        .unwrap_or_else(|| "n/a".to_string());
    format!(
        "genomes {} | clusters {} | focus 0/1/2 = {}/{}/{} | first tree {first} | exact {exact} | \
         wall {secs:.2}s | peak RSS {mem} | identical blocks {identity} | pieces {} | cache hits {} | alignments {}",
        s.genomes,
        s.clusters,
        s.focus_counts[0],
        s.focus_counts[1],
        s.focus_counts[2],
        s.counters.pieces_compared,
        s.counters.cache_hits,
        s.counters.alignments_run,
    )
}

fn cmd_simulate(args: SimulateArgs) -> i32 {
    let cfg = mtato_sim::SimConfig {
        seed: args.seed,
        reference: match args.reference {
            Some(p) => mtato_sim::RefSpec::File(p),
            None => mtato_sim::RefSpec::Random {
                length: args.length,
                gc: args.gc,
            },
        },
        clusters: args.clusters,
        divergence: args.divergence,
        leaves_per_cluster: args.leaves,
        mean_snps_per_branch: args.mean_snps,
        mean_indels_per_branch: args.mean_indels,
        segment_events: mtato_sim::SegmentEvents {
            mean_gains_per_branch: args.mean_segments,
            mean_losses_per_branch: args.mean_segments,
            ..Default::default()
        },
        contigs: args.contigs,
        ..Default::default()
    };
    if let Err(e) = cfg.validate() {
        eprintln!("mtato: invalid simulation: {e}");
        return 2;
    }
    let started = Instant::now();
    match mtato_sim::simulate(&cfg, &args.output) {
        Ok(out) => {
            println!(
                "wrote {} genomes and truth to {} in {:.2}s",
                out.genomes.len(),
                args.output.display(),
                started.elapsed().as_secs_f64()
            );
            0
        }
        Err(e) => {
            eprintln!("mtato: simulation failed: {e}");
            1
        }
    }
}
