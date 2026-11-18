//! `mtato bench`: compare an output folder with a truth folder (section 11).
//!
//! Everything the experiment needs is measured here rather than in the report
//! script, so that the numbers come from the program that produced them and not
//! from a second implementation that could disagree with it.

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};

use clap::Parser;
use serde::Serialize;

use mtato_core::summary::RunSummary;
use mtato_core::tree::{collapse_short, parse_newick, rf_distance, to_newick, Tree};

#[derive(Parser, Debug)]
#[command(about = "Compare an output folder with a truth folder and report the measurements")]
pub struct BenchArgs {
    /// An M-Tato output folder, containing run.json and the trees.
    #[arg(long)]
    pub output: PathBuf,

    /// A simulator truth folder, containing tree_*.nwk, pairwise_snps.tsv and
    /// clusters.tsv.
    #[arg(long)]
    pub truth: PathBuf,

    /// Where to write the measurements as JSON. Defaults to
    /// `<output>/bench.json`.
    #[arg(short, long)]
    pub json: Option<PathBuf>,

    /// Branches shorter than this are collapsed before comparing trees, so that
    /// a branch of less than half a SNP cannot count as a topological
    /// difference (section 11).
    #[arg(long, default_value_t = 0.5)]
    pub collapse_below: f64,

    /// Do not print the table, only write the JSON.
    #[arg(long)]
    pub quiet: bool,
}

/// Accuracy of the pairwise distances against truth.
#[derive(Debug, Clone, Default, Serialize)]
pub struct DistanceAccuracy {
    pub pairs_in_truth: usize,
    pub pairs_compared: usize,
    pub pairs_missing_from_output: usize,
    pub pairs_exact: usize,
    pub max_abs_error: u64,
    pub sum_abs_error: u64,
    pub sum_true_snps: u64,
    pub sum_reported_snps: u64,
    pub worst_pair: Option<(String, String, u64, u64)>,
}

impl DistanceAccuracy {
    /// Fraction of pairs whose reported distance equals the truth exactly.
    pub fn fraction_exact(&self) -> Option<f64> {
        if self.pairs_compared == 0 {
            return None;
        }
        Some(self.pairs_exact as f64 / self.pairs_compared as f64)
    }

    /// Sum of reported distances over sum of true distances. Below 1 means the
    /// tool is missing differences; section 10.3 T3 asks for at least 99%.
    pub fn sum_ratio(&self) -> Option<f64> {
        if self.sum_true_snps == 0 {
            return None;
        }
        Some(self.sum_reported_snps as f64 / self.sum_true_snps as f64)
    }
}

/// Accuracy of the clustering against truth.
#[derive(Debug, Clone, Default, Serialize)]
pub struct ClusterAccuracy {
    pub genomes_in_truth: usize,
    pub genomes_compared: usize,
    pub true_clusters: usize,
    pub reported_clusters: usize,
    /// True clusters whose members landed in exactly one reported cluster with
    /// nothing extra in it.
    pub recovered_exactly: usize,
    /// Pairwise co-clustering F1: does the tool put the same pairs together?
    pub pairwise_f1: Option<f64>,
}

/// One tree comparison.
#[derive(Debug, Clone, Serialize)]
pub struct TreeComparison {
    /// Which emitted tree this is, e.g. `tree.nwk` or `cluster_3.nwk`.
    pub source: String,
    /// Which truth tree it was compared against.
    pub truth: String,
    pub leaves: usize,
    pub rf: usize,
    pub splits_emitted: usize,
    pub splits_true: usize,
    pub rf_normalised: Option<f64>,
    /// Emitted minus truth, in minutes:seconds since the run started, taken from
    /// `events.jsonl`. This is the x axis of the section 11 figure.
    pub t_ms: Option<f64>,
    /// Whether the distances behind this tree were approximate.
    pub approximate: bool,
    pub focus: u8,
    /// What the branch lengths mean: `snps` for a cluster tree, or
    /// `subs_per_site` for the full tree (6.9). Recorded because the collapse
    /// threshold is only meaningful in SNP units.
    #[serde(default)]
    pub units: String,
}

/// Everything `mtato bench` measures.
#[derive(Debug, Clone, Serialize)]
pub struct BenchReport {
    pub output: PathBuf,
    pub truth: PathBuf,
    pub run: Option<RunSummary>,
    pub distance: DistanceAccuracy,
    pub clusters: ClusterAccuracy,
    pub trees: Vec<TreeComparison>,
    /// The accuracy-over-time series: every tree the run emitted, in order.
    pub timeline: Vec<TreeComparison>,
    pub problems: Vec<String>,
}

pub fn cmd_bench(args: BenchArgs) -> i32 {
    match bench(&args) {
        Ok(report) => {
            let path = args
                .json
                .clone()
                .unwrap_or_else(|| args.output.join("bench.json"));
            if let Some(parent) = path.parent() {
                if let Err(e) = std::fs::create_dir_all(parent) {
                    eprintln!("mtato bench: cannot create {}: {e}", parent.display());
                    return 1;
                }
            }
            match std::fs::File::create(&path) {
                Ok(mut f) => {
                    if let Err(e) = f.write_all(
                        serde_json::to_string_pretty(&report)
                            .unwrap_or_else(|_| "{}".to_string())
                            .as_bytes(),
                    ) {
                        eprintln!("mtato bench: cannot write {}: {e}", path.display());
                        return 1;
                    }
                }
                Err(e) => {
                    eprintln!("mtato bench: cannot write {}: {e}", path.display());
                    return 1;
                }
            }
            if !args.quiet {
                print_table(&report, &path);
            }
            // A benchmark that silently found nothing to compare would look like
            // a pass, so that case is a failure.
            if report.problems.iter().any(|p| p.starts_with("no ")) {
                return 1;
            }
            0
        }
        Err(e) => {
            eprintln!("mtato bench: {e}");
            1
        }
    }
}

pub fn bench(args: &BenchArgs) -> Result<BenchReport, String> {
    if !args.output.is_dir() {
        return Err(format!(
            "output folder not found: {}",
            args.output.display()
        ));
    }
    if !args.truth.is_dir() {
        return Err(format!("truth folder not found: {}", args.truth.display()));
    }
    let mut problems = Vec::new();

    let run = read_run_json(&args.output);
    if run.is_none() {
        problems.push("run.json missing or unreadable: timings and counters not measured".into());
    }

    let true_pairs = read_pairwise(&args.truth.join("pairwise_snps.tsv"), &mut problems);
    let reported_pairs = read_reported_pairs(&args.output, &mut problems);
    let distance = compare_distances(&true_pairs, &reported_pairs, &mut problems);

    let (true_clusters, reported_clusters) =
        read_clusters(&args.truth, &args.output, &mut problems);
    let clusters = compare_clusters(&true_clusters, &reported_clusters);

    let timeline = event_timeline(&args.output);
    let trees = compare_trees(
        &args.output,
        &args.truth,
        args.collapse_below,
        &timeline,
        &mut problems,
    );

    Ok(BenchReport {
        output: args.output.clone(),
        truth: args.truth.clone(),
        run,
        distance,
        clusters,
        trees,
        timeline,
        problems,
    })
}

// ---------------------------------------------------------------------------
// Reading the two sides
// ---------------------------------------------------------------------------

fn read_run_json(dir: &Path) -> Option<RunSummary> {
    let txt = std::fs::read_to_string(dir.join("run.json")).ok()?;
    serde_json::from_str(&txt).ok()
}

/// `truth/pairwise_snps.tsv`: genome_a, genome_b, snps.
fn read_pairwise(path: &Path, problems: &mut Vec<String>) -> HashMap<(String, String), u64> {
    let mut out = HashMap::new();
    let Ok(txt) = std::fs::read_to_string(path) else {
        problems.push(format!("no pairwise truth at {}", path.display()));
        return out;
    };
    for (i, line) in txt.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if i == 0 && (f[0] == "genome_a" || f[0].starts_with('#')) {
            continue;
        }
        if f.len() < 3 {
            continue;
        }
        let Ok(s) = f[2].trim().parse::<f64>() else {
            continue;
        };
        out.insert(pair_key(f[0], f[1]), s.round() as u64);
    }
    if out.is_empty() {
        problems.push(format!(
            "no pairwise distances parsed from {}",
            path.display()
        ));
    }
    out
}

/// An unordered pair key, so `a,b` and `b,a` are the same entry.
fn pair_key(a: &str, b: &str) -> (String, String) {
    if a <= b {
        (a.trim().to_string(), b.trim().to_string())
    } else {
        (b.trim().to_string(), a.trim().to_string())
    }
}

/// Every `cluster_<id>.dist.tsv` in an output folder, keyed by genome name.
fn read_reported_pairs(dir: &Path, problems: &mut Vec<String>) -> HashMap<(String, String), u64> {
    let mut out = HashMap::new();
    let Ok(entries) = std::fs::read_dir(dir) else {
        problems.push(format!("cannot read output folder {}", dir.display()));
        return out;
    };
    let mut files: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .map(|n| {
                    let n = n.to_string_lossy();
                    n.starts_with("cluster_") && n.ends_with(".dist.tsv")
                })
                .unwrap_or(false)
        })
        .collect();
    files.sort();
    if files.is_empty() {
        problems.push("no cluster_*.dist.tsv files in the output folder".into());
    }
    for f in files {
        let Ok(txt) = std::fs::read_to_string(&f) else {
            continue;
        };
        for (i, line) in txt.lines().enumerate() {
            let cols: Vec<&str> = line.split('\t').collect();
            if i == 0 || cols.len() < 5 {
                continue;
            }
            // genome_a name_a genome_b name_b snps ...
            let Ok(snps) = cols[4].trim().parse::<u64>() else {
                continue;
            };
            out.insert(pair_key(cols[1], cols[3]), snps);
        }
    }
    out
}

fn read_clusters(
    truth_dir: &Path,
    out_dir: &Path,
    problems: &mut Vec<String>,
) -> (HashMap<String, String>, HashMap<String, String>) {
    let mut t = HashMap::new();
    if let Ok(txt) = std::fs::read_to_string(truth_dir.join("clusters.tsv")) {
        for (i, line) in txt.lines().enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if i == 0 || f.len() < 2 {
                continue;
            }
            t.insert(f[0].trim().to_string(), f[1].trim().to_string());
        }
    } else {
        problems.push("no truth clusters.tsv".into());
    }

    let mut r = HashMap::new();
    if let Ok(txt) = std::fs::read_to_string(out_dir.join("genomes.tsv")) {
        for (i, line) in txt.lines().enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if i == 0 || f.len() < 8 {
                continue;
            }
            // id name path length contigs n50 gc cluster ...
            r.insert(f[1].trim().to_string(), f[7].trim().to_string());
        }
    } else {
        problems.push("no output genomes.tsv".into());
    }
    (t, r)
}

// ---------------------------------------------------------------------------
// Comparing
// ---------------------------------------------------------------------------

fn compare_distances(
    truth: &HashMap<(String, String), u64>,
    reported: &HashMap<(String, String), u64>,
    problems: &mut Vec<String>,
) -> DistanceAccuracy {
    let mut acc = DistanceAccuracy {
        pairs_in_truth: truth.len(),
        ..Default::default()
    };
    if truth.is_empty() {
        return acc;
    }
    let mut keys: Vec<&(String, String)> = truth.keys().collect();
    keys.sort();
    for k in keys {
        let tv = truth[k];
        acc.sum_true_snps += tv;
        match reported.get(k) {
            Some(rv) => {
                acc.pairs_compared += 1;
                acc.sum_reported_snps += rv;
                let err = rv.abs_diff(tv);
                acc.sum_abs_error += err;
                if err == 0 {
                    acc.pairs_exact += 1;
                }
                if err > acc.max_abs_error {
                    acc.max_abs_error = err;
                    acc.worst_pair = Some((k.0.clone(), k.1.clone(), tv, *rv));
                }
            }
            None => acc.pairs_missing_from_output += 1,
        }
    }
    if acc.pairs_missing_from_output > 0 {
        problems.push(format!(
            "{} of {} true pairs were not reported by the run (they may sit in different clusters)",
            acc.pairs_missing_from_output, acc.pairs_in_truth
        ));
    }
    acc
}

fn compare_clusters(
    truth: &HashMap<String, String>,
    reported: &HashMap<String, String>,
) -> ClusterAccuracy {
    let mut acc = ClusterAccuracy {
        genomes_in_truth: truth.len(),
        true_clusters: unique_values(truth).len(),
        reported_clusters: unique_values(reported).len(),
        ..Default::default()
    };
    if truth.is_empty() {
        return acc;
    }
    // Group by true cluster and check each one landed whole and alone.
    let mut by_true: HashMap<&str, Vec<&str>> = HashMap::new();
    for (g, c) in truth {
        by_true.entry(c.as_str()).or_default().push(g.as_str());
    }
    for members in by_true.values() {
        let reported_ids: HashSet<&str> = members
            .iter()
            .filter_map(|g| reported.get(*g).map(|s| s.as_str()))
            .collect();
        if reported_ids.len() != 1 {
            continue;
        }
        let rid = *reported_ids.iter().next().unwrap();
        let extra = reported
            .iter()
            .filter(|(g, c)| c.as_str() == rid && !truth.contains_key(*g))
            .count();
        let members_not_here = members
            .iter()
            .filter(|g| reported.get(**g).map(|s| s.as_str()) != Some(rid))
            .count();
        if extra == 0 && members_not_here == 0 {
            acc.recovered_exactly += 1;
        }
    }

    // Pairwise co-clustering F1 over the genomes present on both sides.
    let shared: Vec<&String> = truth
        .keys()
        .filter(|g| reported.contains_key(g.as_str()))
        .collect();
    acc.genomes_compared = shared.len();
    if shared.len() >= 2 {
        let (mut tp, mut fp, mut fn_) = (0u64, 0u64, 0u64);
        let mut sorted = shared.clone();
        sorted.sort();
        for i in 0..sorted.len() {
            for j in (i + 1)..sorted.len() {
                let (a, b) = (sorted[i], sorted[j]);
                let same_true = truth[a] == truth[b];
                let same_rep = reported[a] == reported[b];
                match (same_true, same_rep) {
                    (true, true) => tp += 1,
                    (false, true) => fp += 1,
                    (true, false) => fn_ += 1,
                    (false, false) => {}
                }
            }
        }
        let precision = if tp + fp == 0 {
            0.0
        } else {
            tp as f64 / (tp + fp) as f64
        };
        let recall = if tp + fn_ == 0 {
            0.0
        } else {
            tp as f64 / (tp + fn_) as f64
        };
        acc.pairwise_f1 = if precision + recall == 0.0 {
            Some(0.0)
        } else {
            Some(2.0 * precision * recall / (precision + recall))
        };
    }
    acc
}

fn unique_values(m: &HashMap<String, String>) -> HashSet<&str> {
    m.values().map(|s| s.as_str()).collect()
}

// ---------------------------------------------------------------------------
// Trees
// ---------------------------------------------------------------------------

/// Every tree the run emitted, in time order, read from `events.jsonl`.
fn event_timeline(output: &Path) -> Vec<TreeComparison> {
    let path = output.join("events.jsonl");
    let Ok(events) = mtato_core::events::read_event_log(&path) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for env in events {
        match &env.event {
            mtato_core::events::Event::BackboneTree { newick, .. } => {
                out.push(TreeComparison {
                    source: "backbone (live)".into(),
                    truth: String::new(),
                    leaves: leaf_count_of_newick(newick),
                    rf: 0,
                    splits_emitted: 0,
                    splits_true: 0,
                    rf_normalised: None,
                    t_ms: Some(env.t_ms),
                    approximate: true,
                    focus: 0,
                    units: "subs_per_site".into(),
                });
            }
            mtato_core::events::Event::ClusterFocus(info) => {
                out.push(TreeComparison {
                    source: format!("cluster {} (live)", info.cluster),
                    truth: String::new(),
                    leaves: leaf_count_of_newick(&info.newick),
                    rf: 0,
                    splits_emitted: 0,
                    splits_true: 0,
                    rf_normalised: None,
                    t_ms: Some(env.t_ms),
                    approximate: info.approximate,
                    focus: info.focus,
                    units: if info.focus == 0 {
                        "subs_per_site"
                    } else if info.focus == 1 {
                        "diff_pieces"
                    } else {
                        "snps"
                    }
                    .into(),
                });
            }
            _ => {}
        }
    }
    out
}

fn leaf_count_of_newick(s: &str) -> usize {
    parse_newick(s).map(|t| t.leaf_count()).unwrap_or(0)
}

/// Compare each emitted tree with the truth trees.
///
/// Truth trees are per cluster, emitted trees are per reported cluster, and the
/// numbering need not agree, so trees are matched by their leaf sets rather than
/// by name.
fn compare_trees(
    output: &Path,
    truth: &Path,
    collapse_below: f64,
    timeline: &[TreeComparison],
    problems: &mut Vec<String>,
) -> Vec<TreeComparison> {
    let mut out = Vec::new();

    let truth_trees = read_trees(truth, |n| n.starts_with("tree"));
    if truth_trees.is_empty() {
        problems.push("no truth tree_*.nwk files to compare against".into());
    }
    let mut emitted = read_trees(output, |n| n.starts_with("cluster_") && n.ends_with(".nwk"));
    if let Some(full) = read_one_tree(&output.join("tree.nwk")) {
        emitted.push(("tree.nwk".to_string(), full));
    }
    if emitted.is_empty() {
        problems.push("no emitted trees in the output folder".into());
    }

    // A timeline entry carries the time a tree was emitted; the final trees are
    // matched to the last emission of the same cluster so the figure has an x.
    let time_of = |source: &str| -> Option<f64> {
        let key = source.trim_end_matches(".nwk");
        timeline
            .iter()
            .rev()
            .find(|t| t.source.starts_with(key) || key.starts_with("cluster"))
            .and_then(|t| t.t_ms)
    };

    for (tname, ttree) in &truth_trees {
        let tleaves: HashSet<String> = ttree.leaf_names().into_iter().collect();
        if tleaves.is_empty() {
            continue;
        }
        // The emitted tree that covers the same genomes. Per-cluster trees win
        // ties over the full `tree.nwk`, because a cluster tree is in SNP counts
        // while the full tree is in substitutions per site (6.9) — and the
        // collapse threshold below is expressed in SNPs. Collapsing a
        // subs/site tree at 0.5 "SNPs" would flatten it to a star and report a
        // meaningless Robinson-Foulds distance.
        let best = emitted.iter().max_by_key(|(ename, et)| {
            let leaves: HashSet<String> = et.leaf_names().into_iter().collect();
            (
                leaves.intersection(&tleaves).count(),
                ename.starts_with("cluster_"),
            )
        });
        let Some((ename, etree)) = best else {
            continue;
        };
        let eleaves: HashSet<String> = etree.leaf_names().into_iter().collect();
        if eleaves != tleaves {
            problems.push(format!(
                "truth {tname} has {} leaves but the closest emitted tree {ename} has {}; \
                 comparing on the shared {} only",
                tleaves.len(),
                eleaves.len(),
                tleaves.intersection(&eleaves).count()
            ));
        }
        let tc = collapse_short(ttree.clone(), collapse_below);
        // Only a tree in SNP counts can be collapsed at a SNP threshold.
        let snp_units = ename.starts_with("cluster_");
        let ec = collapse_short(etree.clone(), if snp_units { collapse_below } else { 0.0 });
        // Restricting both to the shared leaf set keeps the comparison fair when
        // the clustering did not agree.
        let shared: Vec<String> = tleaves.intersection(&eleaves).cloned().collect();
        let (tc, ec) = (prune(&tc, &shared), prune(&ec, &shared));
        let (rf, na, nb) = rf_distance(&tc, &ec);
        out.push(TreeComparison {
            source: ename.clone(),
            truth: tname.clone(),
            leaves: shared.len(),
            rf,
            splits_emitted: na,
            splits_true: nb,
            rf_normalised: if na + nb == 0 {
                None
            } else {
                Some(rf as f64 / (na + nb) as f64)
            },
            t_ms: time_of(ename),
            approximate: false,
            focus: 2,
            units: if snp_units { "snps" } else { "subs_per_site" }.to_string(),
        });
    }
    out
}

fn read_trees(dir: &Path, accept: impl Fn(&str) -> bool) -> Vec<(String, Tree)> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut paths: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.extension().map(|e| e == "nwk").unwrap_or(false)
                && p.file_name()
                    .map(|n| accept(&n.to_string_lossy()))
                    .unwrap_or(false)
        })
        .collect();
    paths.sort();
    paths
        .iter()
        .filter_map(|p| read_one_tree(p).map(|t| (file_stem(p), t)))
        .collect()
}

fn read_one_tree(path: &Path) -> Option<Tree> {
    let txt = std::fs::read_to_string(path).ok()?;
    parse_newick(txt.trim()).ok()
}

fn file_stem(p: &Path) -> String {
    p.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.to_string_lossy().into_owned())
}

/// Keep only the named leaves, dropping internal nodes that end up with one
/// child. Used to compare two trees over their shared leaf set.
fn prune(t: &Tree, keep: &[String]) -> Tree {
    let set: HashSet<&str> = keep.iter().map(|s| s.as_str()).collect();
    prune_rec(t, &set).unwrap_or_default()
}

fn prune_rec(t: &Tree, keep: &HashSet<&str>) -> Option<Tree> {
    match t {
        Tree::Leaf { name, len } => {
            if keep.contains(name.as_str()) {
                Some(Tree::Leaf {
                    name: name.clone(),
                    len: *len,
                })
            } else {
                None
            }
        }
        Tree::Internal { len, children } => {
            let kept: Vec<Tree> = children.iter().filter_map(|c| prune_rec(c, keep)).collect();
            match kept.len() {
                0 => None,
                1 => Some(kept.into_iter().next().unwrap()),
                _ => Some(Tree::Internal {
                    len: *len,
                    children: kept,
                }),
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

fn print_table(r: &BenchReport, json_path: &Path) {
    let mut o = std::io::stdout();
    let _ = writeln!(o, "M-Tato benchmark");
    let _ = writeln!(o, "  output: {}", r.output.display());
    let _ = writeln!(o, "  truth:  {}", r.truth.display());
    if let Some(s) = &r.run {
        let t = &s.timings;
        let fmt = |v: Option<u64>| match v {
            Some(ms) => format!("{:.2}s", ms as f64 / 1000.0),
            None => "not reached".to_string(),
        };
        let _ = writeln!(
            o,
            "  version {} | genomes {} | clusters {} | first tree {} | exact {} | peak RSS {}",
            s.version,
            s.genomes,
            s.clusters,
            fmt(t.first_tree_ms),
            fmt(t.exact_tree_ms),
            s.peak_memory_bytes
                .map(|b| format!("{:.0} MB", b as f64 / 1e6))
                .unwrap_or_else(|| "unknown".to_string())
        );
        let _ = writeln!(
            o,
            "  counters: blocks looked up {} | matched {} | pieces {} | cache hits {} | alignments {} | bases aligned {}",
            s.counters.blocks_looked_up,
            s.counters.total_matched(),
            s.counters.pieces_compared,
            s.counters.cache_hits,
            s.counters.alignments_run,
            s.counters.bases_aligned
        );
    }
    let d = &r.distance;
    let _ = writeln!(o, "  distances:");
    let _ = writeln!(
        o,
        "    pairs in truth {} | compared {} | missing {}",
        d.pairs_in_truth, d.pairs_compared, d.pairs_missing_from_output
    );
    let _ = writeln!(
        o,
        "    exact {} | max |error| {} | mean |error| {}",
        d.pairs_exact,
        d.max_abs_error,
        if d.pairs_compared == 0 {
            "n/a".to_string()
        } else {
            format!("{:.3}", d.sum_abs_error as f64 / d.pairs_compared as f64)
        }
    );
    let _ = writeln!(
        o,
        "    fraction exact {} | sum reported / sum true {}",
        pct(d.fraction_exact()),
        pct(d.sum_ratio())
    );
    if let Some((a, b, tv, rv)) = &d.worst_pair {
        let _ = writeln!(o, "    worst pair: {a} vs {b}: true {tv}, reported {rv}");
    }
    let c = &r.clusters;
    let _ = writeln!(o, "  clusters:");
    let _ = writeln!(
        o,
        "    true {} | reported {} | recovered exactly {} | pairwise F1 {}",
        c.true_clusters,
        c.reported_clusters,
        c.recovered_exactly,
        pct(c.pairwise_f1)
    );
    if !r.trees.is_empty() {
        let _ = writeln!(
            o,
            "  trees (branches shorter than the collapse threshold removed):"
        );
        for t in &r.trees {
            let _ = writeln!(
                o,
                "    {:<18} vs {:<18} leaves {:>4}  RF {:>4}  normalised {}",
                t.source,
                t.truth,
                t.leaves,
                t.rf,
                pct(t.rf_normalised)
            );
        }
    }
    if !r.problems.is_empty() {
        let _ = writeln!(o, "  problems:");
        for p in &r.problems {
            let _ = writeln!(o, "    - {p}");
        }
    }
    let _ = writeln!(o, "  written to {}", json_path.display());
    let _ = o.flush();
}

fn pct(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{:.4}", x),
        None => "n/a".to_string(),
    }
}

#[allow(dead_code)]
fn newick_of(t: &Tree) -> String {
    to_newick(t)
}
