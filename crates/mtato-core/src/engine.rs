//! The order of work (6.11) and the engine handle of section 7.
//!
//! The engine is a library that runs in background threads and communicates
//! only through messages, so the CLI and the GUI are interchangeable consumers
//! (section 0 rule 2).
//!
//! Determinism (section 0 rule 3) is preserved by construction: scanning and
//! refinement use rayon for speed but collect into index-addressed slots and
//! are merged in genome-id order, so the thread count changes nothing in the
//! output.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Instant;

use crossbeam_channel::{unbounded, Receiver, Sender};
use rayon::prelude::*;
use rayon::ThreadPool;

use crate::cluster::{backbone_tree, cluster_genomes, unique_names, Cluster, Clustering};
use crate::compare::{
    compare_member, coverage, piece_variant, Focus1, MemberComparison, PieceCache,
};
use crate::counters::Counters;
use crate::events::{
    event_bus, stage, ClusterFocusInfo, ClusterInfo, Command, Envelope, Event, EventSink,
    GenomeScannedInfo, MapUpdateInfo, MemberNumbers, PieceView,
};
use crate::fasta::GenomeData;
use crate::genome::GenomeRecord;
use crate::interval::Iv;
use crate::params::Parameters;
use crate::pivot::PivotIndex;
use crate::profile::{
    cluster_tree, concat_pivot, contig_offsets, distance_matrix, expand_groups, group_identical,
    star_tree, DistanceMatrix, Profile, TreeBasis,
};
use crate::summary::{RunSummary, StageTimings};
use crate::tree::{graft, to_newick, Tree};
use crate::variant::{DiffParams, DiffResult};

/// Everything the engine needs to start.
#[derive(Debug, Clone)]
pub struct RunConfig {
    pub inputs: Vec<PathBuf>,
    pub params: Parameters,
    /// Where the output files go. `None` runs the engine without writing
    /// anything, which is what the tests and a GUI live preview use.
    pub outdir: Option<PathBuf>,
}

/// A running engine. Consumers read `events` and may send `commands`.
#[derive(Debug)]
pub struct EngineHandle {
    pub events: Receiver<Envelope>,
    pub commands: Sender<Command>,
    pub worker: std::thread::JoinHandle<Result<RunResults, String>>,
}

/// Everything a finished run produced, in the shape the writers need.
#[derive(Debug, Default)]
pub struct RunResults {
    pub records: Vec<GenomeRecord>,
    /// Display name per genome id, de-duplicated for Newick output.
    pub names: Vec<String>,
    pub clusters: Vec<Cluster>,
    /// Sketch distance from each genome to its cluster pivot.
    pub dist_to_pivot: Vec<f64>,
    pub backbone: Tree,
    pub full_tree: Tree,
    /// Cluster id -> tree in substitutions per site, used inside the full tree.
    pub cluster_trees: HashMap<u32, Tree>,
    /// Cluster id -> tree in SNP counts, written as `cluster_<id>.nwk`.
    pub cluster_snp_trees: HashMap<u32, Tree>,
    pub cluster_focus: HashMap<u32, u8>,
    pub matrices: HashMap<u32, DistanceMatrix>,
    /// Genome id -> exact comparison, for members that reached Focus 2.
    pub comparisons: HashMap<u32, MemberComparison>,
    /// Genome id -> its contig names, so that writing the variant files never
    /// has to open a FASTA file again.
    pub contig_names: HashMap<u32, Vec<String>>,
    /// Cluster id -> the contig names of its pivot.
    pub pivot_contig_names: HashMap<u32, Vec<String>>,
    pub counters: Counters,
    pub timings: StageTimings,
    pub notes: Vec<String>,
    /// Cluster id -> total length of its pivot, the denominator of 6.9.
    pub pivot_len: HashMap<u32, u64>,
    pub params: Parameters,
    pub cancelled: bool,
}

/// Start the engine on a background thread.
pub fn spawn(cfg: RunConfig) -> EngineHandle {
    let (sink, events) = event_bus();
    let (commands, cmd_rx) = unbounded::<Command>();
    let worker = std::thread::Builder::new()
        .name("mtato-engine".into())
        .spawn(move || run(cfg, sink, cmd_rx))
        .expect("failed to spawn the engine thread");
    EngineHandle {
        events,
        commands,
        worker,
    }
}

/// Run the engine on the calling thread, collecting every event.
///
/// Convenient for tests and for the CLI's `--no-stream` path.
pub fn run_collecting(cfg: RunConfig) -> Result<(RunResults, Vec<Envelope>), String> {
    let (sink, rx) = event_bus();
    let results = run(cfg, sink, unbounded::<Command>().1)?;
    let mut events = Vec::new();
    while let Ok(e) = rx.try_recv() {
        events.push(e);
    }
    Ok((results, events))
}

/// Shared state between the engine thread and the command receiver.
#[derive(Debug)]
struct Controls {
    cancel: Arc<AtomicBool>,
    cluster_dist_bits: Arc<AtomicU64>,
}

impl Controls {
    fn new(params: &Parameters) -> Self {
        Self {
            cancel: Arc::new(AtomicBool::new(false)),
            cluster_dist_bits: Arc::new(AtomicU64::new(params.cluster_dist.to_bits())),
        }
    }

    fn cluster_dist(&self) -> f64 {
        f64::from_bits(self.cluster_dist_bits.load(Ordering::Relaxed))
    }

    fn cancelled(&self) -> bool {
        self.cancel.load(Ordering::Relaxed)
    }

    /// Drain pending commands. Returns the cluster ids to move to the front and
    /// whether the threshold changed, which forces a re-cluster.
    fn poll(&self, rx: &Receiver<Command>, sink: &EventSink) -> (Vec<u32>, bool) {
        let mut prioritised = Vec::new();
        let mut recluster = false;
        while let Ok(cmd) = rx.try_recv() {
            match cmd {
                Command::Cancel => {
                    self.cancel.store(true, Ordering::Relaxed);
                    sink.warn("cancellation requested; stopping after the current genome");
                }
                Command::Prioritize { cluster } => prioritised.push(cluster),
                Command::SetClusterDist { value } => {
                    if value > 0.0 && value < 1.0 {
                        self.cluster_dist_bits
                            .store(value.to_bits(), Ordering::Relaxed);
                        recluster = true;
                        sink.warn(format!("cluster_dist changed to {value}; re-clustering"));
                    } else {
                        sink.warn(format!("ignoring cluster_dist {value}: must be in (0, 1)"));
                    }
                }
            }
        }
        (prioritised, recluster)
    }
}

/// Run the whole pipeline. Both the CLI and the GUI call this.
pub fn run(
    cfg: RunConfig,
    sink: EventSink,
    cmd_rx: Receiver<Command>,
) -> Result<RunResults, String> {
    let started = Instant::now();
    let controls = Controls::new(&cfg.params);
    let mut params = cfg.params.clone();
    params.validate()?;
    if params.memory_budget_bytes == 0 {
        params.memory_budget_bytes = Parameters::default_memory_budget();
    }

    // A per-run pool, not the global one: a second run in the same process (the
    // GUI, or a test suite) must still honour its own thread count.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.resolved_threads())
        .thread_name(|i| format!("mtato-worker-{i}"))
        .build()
        .map_err(|e| format!("building the thread pool: {e}"))?;

    let mut results = RunResults {
        params: params.clone(),
        ..Default::default()
    };

    // ---- discover ----
    sink.progress(stage::DISCOVER, 0, 0);
    let files = crate::fasta::collect_inputs(&cfg.inputs).map_err(|e| e.to_string())?;
    if files.is_empty() {
        return Err(format!(
            "no genome files found in the given inputs: {:?}",
            cfg.inputs
        ));
    }
    sink.progress(stage::DISCOVER, files.len() as u64, files.len() as u64);

    // ---- scan (6.11.1) ----
    let scanned = scan_all(&pool, &files, &params, &sink, &controls)?;
    results.records = scanned.records;
    results.contig_names = scanned.contig_names;
    let mut notes = scanned.notes;
    let mut scan_ms_total = scanned.elapsed_ms;

    // The small-genome rule of 6.1 can only be judged once the lengths are
    // known, so it is applied after a first scan and the scan is repeated with
    // the denser sampling. For bacterial assemblies, which is the design
    // target, this second pass never runs.
    let median = median_length(&results.records);
    let mut effective = params.clone();
    if let Some(msg) = effective.apply_small_genome_rule(median) {
        sink.warn(msg.clone());
        notes.push(msg);
        let rescan = scan_all(&pool, &files, &effective, &sink, &controls)?;
        results.records = rescan.records;
        results.contig_names = rescan.contig_names;
        scan_ms_total = scan_ms_total.saturating_add(rescan.elapsed_ms);
        params = effective.clone();
        results.params = params.clone();
    }
    results.notes = notes;
    results.timings.scan_ms = Some(scan_ms_total);
    // `genomes_scanned` counts the genomes that really contributed, so it is
    // set once, here: after unreadable files have been dropped and after the
    // small-genome rescan of 6.1 (D8) has replaced the records. Setting it per
    // pass would double-count a rescan; setting it from `files.len()` would
    // count a file that could not be read. The per-cluster counters merged in
    // during refinement leave this field at zero, so the value survives.
    results.counters.genomes_scanned = results.records.len() as u64;
    if controls.cancelled() {
        return Ok(finish(results, &sink, started));
    }

    // ---- cluster (6.11.2) ----
    let t = Instant::now();
    let (names, name_warnings) = unique_names(&results.records);
    results.names = names;
    for w in name_warnings {
        sink.warn(w.clone());
        results.notes.push(w);
    }
    let mut clustering = cluster_genomes(&results.records, &params);
    emit_clusters(&sink, &clustering, &results.names);
    results.timings.cluster_ms = Some(t.elapsed().as_millis() as u64);

    // ---- backbone tree and the Focus 0 picture (6.11.2) ----
    let t = Instant::now();
    results.backbone = backbone_tree(
        &clustering.clusters,
        &results.records,
        &results.names,
        &params,
    );
    let pivot_names: Vec<String> = clustering
        .clusters
        .iter()
        .map(|c| results.names[c.pivot as usize].clone())
        .collect();
    sink.emit(Event::BackboneTree {
        newick: to_newick(&results.backbone),
        pivots: pivot_names,
    });
    results.timings.backbone_ms = Some(t.elapsed().as_millis() as u64);
    // "A usable picture early is worth more than waiting for the exact one"
    // (section 1): this is the moment the first blurry tree exists.
    results.timings.first_tree_ms = Some(started.elapsed().as_millis() as u64);
    assign_clusters(&mut results, &clustering);
    emit_focus0(&sink, &clustering, &results);

    // ---- refine clusters (6.11.3) ----
    let t = Instant::now();
    refine_clusters(
        &pool,
        &mut results,
        &mut clustering,
        &params,
        &sink,
        &controls,
        &cmd_rx,
    )?;
    results.timings.refine_ms = Some(t.elapsed().as_millis() as u64);

    results.clusters = clustering.clusters.clone();
    results.dist_to_pivot = clustering.dist_to_pivot.clone();

    // ---- the full tree (6.9) ----
    results.full_tree = build_full_tree(&results);
    if !results.cancelled && results.cluster_focus.values().all(|&f| f >= 2 || f == 0) {
        results.timings.exact_tree_ms = Some(started.elapsed().as_millis() as u64);
    }

    // ---- output files (section 9) ----
    if let Some(dir) = &cfg.outdir {
        let t = Instant::now();
        sink.progress(stage::OUTPUT, 0, 1);
        // The tables are written first and `run.json` last, so the timings it
        // carries are complete. Writing everything in one call used to
        // serialise `run.json` before `output_ms` was known and before
        // `total_ms` was ever set, leaving both `null` in a file that section
        // 9 says holds "timings per stage" and that section 11 reads
        // time-to-first-tree and time-to-exact from.
        crate::output::write_tables(dir, &results)
            .map_err(|e| format!("writing output to {}: {e}", dir.display()))?;
        results.timings.output_ms = Some(t.elapsed().as_millis() as u64);
        // `total_ms` in `run.json` necessarily excludes the writing of
        // `run.json` itself: that step has not happened yet, and the file
        // cannot report its own duration. It is the whole run up to this point,
        // which is the honest figure, and it must not be "fixed" later by
        // rewriting the file. The `Finished` event's `total_ms`, set in
        // `finish()`, does include the `run.json` write.
        results.timings.total_ms = Some(started.elapsed().as_millis() as u64);
        let summary = summary_of(&results);
        crate::output::write_run_json(dir, &results, &summary)
            .map_err(|e| format!("writing output to {}: {e}", dir.display()))?;
        sink.progress(stage::OUTPUT, 1, 1);
    }

    Ok(finish(results, &sink, started))
}

fn finish(mut results: RunResults, sink: &EventSink, started: Instant) -> RunResults {
    results.timings.total_ms = Some(started.elapsed().as_millis() as u64);
    let summary = summary_of(&results);
    sink.emit(Event::Counters(results.counters.clone()));
    sink.emit(Event::Finished(Box::new(summary)));
    sink.flush_log();
    results
}

/// Build the `Finished` payload from the results.
pub fn summary_of(results: &RunResults) -> RunSummary {
    let mut s = RunSummary::new(crate::VERSION, results.params.clone());
    s.notes = results.notes.clone();
    s.genomes = results.records.len() as u32;
    s.clusters = results.clusters.len() as u32;
    s.total_length = results.records.iter().map(|r| r.stats.length).sum();
    for r in &results.records {
        s.focus_counts[r.focus.min(2) as usize] += 1;
    }
    s.timings = results.timings.clone();
    s.counters = results.counters.clone();
    s.peak_memory_bytes = crate::summary::peak_memory_bytes();
    s.cancelled = results.cancelled;
    s
}

fn median_length(records: &[GenomeRecord]) -> u64 {
    let mut lens: Vec<u64> = records.iter().map(|r| r.stats.length).collect();
    if lens.is_empty() {
        return 0;
    }
    lens.sort_unstable();
    let mid = lens.len() / 2;
    if lens.len() % 2 == 1 {
        return lens[mid];
    }
    // An even count has no single middle element, so average the two. Taking
    // the upper one alone made the 500 kb threshold of 6.1 depend on whether
    // the genome count happened to be even: for {499_998, 500_001} the true
    // median is 499_999.5, which is under the threshold and must fire the
    // small-genome rule, but the upper element is 500_001 and does not.
    // Integer division rounds down, which cannot make the rule fire when the
    // true median is at or above the threshold.
    (lens[mid - 1] + lens[mid]) / 2
}

// ---------------------------------------------------------------------------
// Scan
// ---------------------------------------------------------------------------

/// What one scan pass produced (6.11.1).
struct ScanOutcome {
    /// One record per genome that could be read. Ids are renumbered densely,
    /// so an id and a position in this list always agree even when a file
    /// failed to read.
    records: Vec<GenomeRecord>,
    /// Contig names per genome id, captured here so that output writing never
    /// has to reopen the inputs (D9).
    contig_names: HashMap<u32, Vec<String>>,
    /// Read failures and anything else worth telling the user about.
    notes: Vec<String>,
    /// Wall-clock milliseconds this pass took.
    elapsed_ms: u64,
}

fn scan_all(
    pool: &ThreadPool,
    files: &[PathBuf],
    params: &Parameters,
    sink: &EventSink,
    controls: &Controls,
) -> Result<ScanOutcome, String> {
    let t = Instant::now();
    let total = files.len() as u64;
    sink.progress(stage::SCAN, 0, total);
    let done = Arc::new(AtomicU64::new(0));

    // Collected into an index-addressed vector, so the merge order below does
    // not depend on which worker finished first.
    let scanned: Vec<Result<(GenomeRecord, Vec<String>), String>> = pool.install(|| {
        files
            .par_iter()
            .enumerate()
            .map(|(i, path)| {
                if controls.cancelled() {
                    return Err("cancelled".to_string());
                }
                let (data, scan) = crate::genome::scan_file(path, params)
                    .map_err(|e| format!("reading {}: {e}", path.display()))?;
                let names: Vec<String> = data.contigs.iter().map(|c| c.name.clone()).collect();
                let name = crate::fasta::genome_name(path);
                let rec = GenomeRecord::new(
                    i as u32,
                    name.clone(),
                    path.clone(),
                    scan.stats,
                    scan.sketch,
                );
                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                if n % 16 == 0 || n == total {
                    sink.progress(stage::SCAN, n, total);
                }
                sink.emit(Event::GenomeScanned(GenomeScannedInfo {
                    id: rec.id,
                    name,
                    length: rec.stats.length,
                    contigs: rec.stats.contigs,
                    n50: rec.stats.n50,
                    gc: rec.stats.gc,
                    sketch_size: rec.sketch.len(),
                }));
                Ok((rec, names))
            })
            .collect()
    });

    let mut pairs: Vec<(GenomeRecord, Vec<String>)> = Vec::with_capacity(files.len());
    let mut notes = Vec::new();
    for r in scanned {
        match r {
            Ok(pair) => pairs.push(pair),
            Err(e) if e == "cancelled" => {}
            Err(e) => {
                sink.warn(e.clone());
                notes.push(e);
            }
        }
    }
    if pairs.is_empty() {
        return Err("no genome could be scanned".to_string());
    }
    // Ids were assigned by input index; renumber densely so that an id and a
    // position in `records` always agree, even if some file failed to read.
    let mut records = Vec::with_capacity(pairs.len());
    let mut contig_names: HashMap<u32, Vec<String>> = HashMap::with_capacity(pairs.len());
    for (i, (mut rec, names)) in pairs.into_iter().enumerate() {
        rec.id = i as u32;
        contig_names.insert(rec.id, names);
        records.push(rec);
    }
    sink.progress(stage::SCAN, total, total);
    Ok(ScanOutcome {
        records,
        contig_names,
        notes,
        elapsed_ms: t.elapsed().as_millis() as u64,
    })
}

// ---------------------------------------------------------------------------
// Clustering events
// ---------------------------------------------------------------------------

fn emit_clusters(sink: &EventSink, c: &Clustering, names: &[String]) {
    let infos: Vec<ClusterInfo> = c
        .clusters
        .iter()
        .map(|cl| ClusterInfo {
            id: cl.id,
            pivot: cl.pivot,
            pivot_name: names[cl.pivot as usize].clone(),
            members: cl.members.clone(),
            member_names: cl
                .members
                .iter()
                .map(|&m| names[m as usize].clone())
                .collect(),
            focus: 0,
            cluster_dist: c.threshold,
        })
        .collect();
    sink.emit(Event::Clustered {
        threshold: c.threshold,
        clusters: infos,
    });
}

fn assign_clusters(results: &mut RunResults, c: &Clustering) {
    for r in results.records.iter_mut() {
        r.cluster = u32::MAX;
        r.is_pivot = false;
        r.focus = 0;
    }
    for cl in &c.clusters {
        for &m in &cl.members {
            if let Some(r) = results.records.get_mut(m as usize) {
                r.cluster = cl.id;
                r.is_pivot = m == cl.pivot;
            }
        }
    }
}

/// The blurry cluster view of 6.4: until a cluster is refined, its members hang
/// from the pivot as a star, each at its sketch distance.
fn emit_focus0(sink: &EventSink, c: &Clustering, results: &RunResults) {
    for cl in &c.clusters {
        let pivot_name = &results.names[cl.pivot as usize];
        let members: Vec<(String, f64)> = cl
            .members
            .iter()
            .filter(|&&m| m != cl.pivot)
            .map(|&m| {
                (
                    results.names[m as usize].clone(),
                    c.dist_to_pivot[m as usize],
                )
            })
            .collect();
        let tree = star_tree(pivot_name, &members);
        sink.emit(Event::ClusterFocus(ClusterFocusInfo {
            cluster: cl.id,
            focus: 0,
            approximate: true,
            newick: to_newick(&tree),
            members: cl
                .members
                .iter()
                .map(|&m| MemberNumbers {
                    id: m,
                    name: results.names[m as usize].clone(),
                    diff_pieces: 0,
                    snps: 0,
                    indels: 0,
                    uncovered_bp: 0,
                    novel_bp: 0,
                })
                .collect(),
        }));
    }
}

// ---------------------------------------------------------------------------
// Refinement
// ---------------------------------------------------------------------------

/// What one cluster's refinement produced.
#[derive(Debug, Clone, Default)]
pub struct ClusterResult {
    pub focus: u8,
    /// Neighbour joining on SNP counts, for `cluster_<id>.nwk`.
    pub snp_tree: Tree,
    /// Neighbour joining on substitutions per site, for the full tree.
    pub subs_tree: Tree,
    pub matrix: DistanceMatrix,
    pub comparisons: Vec<(u32, MemberComparison)>,
    pub counters: Counters,
    pub pivot_len: u64,
    /// Contig names of the cluster's pivot, for the variant files.
    pub pivot_contig_names: Vec<String>,
    /// Members dropped at Focus 2 because their FASTA could not be read. Each
    /// is also emitted as a warning; carrying it here lets `apply_cluster_result`
    /// put it in `run.json`'s notes, so a run that silently lost a member is
    /// impossible.
    pub notes: Vec<String>,
}

/// Move the clusters the user selected to the front of the work order.
///
/// 6.11.3: "a cluster the user selects jumps to the front". `order` holds
/// indices into `clusters` while `wanted` holds cluster ids, so the lookup
/// translates one into the other. Walking the requests in reverse and inserting
/// each at position 0 leaves them at the front in the order they were made.
///
/// The element removed from `order` is the element put back at the front:
/// nothing is dropped and nothing is duplicated. Inserting the *position* that
/// was found instead of the cluster index that was removed silently drops the
/// selected cluster from this pass and refines an unrelated one first, which
/// cannot panic (a position is always a valid index here) and so fails quietly.
fn prioritize(order: &mut Vec<usize>, clusters: &[Cluster], wanted: &[u32]) {
    for want in wanted.iter().rev() {
        if let Some(pos) = order.iter().position(|&i| clusters[i].id == *want) {
            let ci = order.remove(pos);
            order.insert(0, ci);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn refine_clusters(
    pool: &ThreadPool,
    results: &mut RunResults,
    clustering: &mut Clustering,
    params: &Parameters,
    sink: &EventSink,
    controls: &Controls,
    cmd_rx: &Receiver<Command>,
) -> Result<(), String> {
    // Results keyed by cluster membership, so that a re-cluster which does not
    // change a cluster keeps what was already computed (section 7).
    let mut done: HashMap<Vec<u32>, ClusterResult> = HashMap::new();
    let mut emitted_focus: HashSet<u32> = HashSet::new();
    let n_clusters_total = clustering.clusters.len();

    loop {
        if controls.cancelled() {
            results.cancelled = true;
            break;
        }
        // Work order: largest cluster first, except that a cluster the user
        // selected jumps to the front (6.11.3).
        let (prioritised, recluster) = controls.poll(cmd_rx, sink);
        if recluster {
            let mut p = params.clone();
            p.cluster_dist = controls.cluster_dist();
            *clustering = cluster_genomes(&results.records, &p);
            emit_clusters(sink, clustering, &results.names);
            assign_clusters(results, clustering);
            emit_focus0(sink, clustering, results);
            emitted_focus.clear();
            results.params.cluster_dist = p.cluster_dist;
        }

        let mut order: Vec<usize> = (0..clustering.clusters.len()).collect();
        order.sort_by(|&a, &b| {
            clustering.clusters[b]
                .size()
                .cmp(&clustering.clusters[a].size())
                .then(a.cmp(&b))
        });
        prioritize(&mut order, &clustering.clusters, &prioritised);

        let mut worked = false;
        for ci in order {
            if controls.cancelled() {
                results.cancelled = true;
                break;
            }
            let cluster = clustering.clusters[ci].clone();
            let signature = cluster.members.clone();
            if done.contains_key(&signature) {
                continue;
            }
            worked = true;
            let res = match done.remove(&cluster.members.clone()) {
                Some(prev) => prev,
                None => refine_one_cluster(pool, results, &cluster, params, sink, controls)?,
            };
            apply_cluster_result(results, &cluster, &res);
            if emitted_focus.insert(cluster.id) {
                sink.progress(
                    stage::REFINE,
                    emitted_focus.len() as u64,
                    n_clusters_total as u64,
                );
            }
            done.insert(signature, res);
            // Re-poll between clusters so a Prioritize takes effect quickly.
            break;
        }
        if !worked {
            break;
        }
    }
    Ok(())
}

fn refine_one_cluster(
    pool: &ThreadPool,
    results: &RunResults,
    cluster: &Cluster,
    params: &Parameters,
    sink: &EventSink,
    controls: &Controls,
) -> Result<ClusterResult, String> {
    let pivot_id = cluster.pivot;
    let pivot_rec = &results.records[pivot_id as usize];
    let members: Vec<u32> = cluster.non_pivot().collect();

    if members.is_empty() {
        // "Clusters with a single genome need no index" (6.5). It still needs a
        // tree. `build_full_tree` grafts the cluster tree over the backbone
        // leaf that names the pivot, and `Tree::default()` is the *empty* tree,
        // so leaving it at the default deleted this genome from `tree.nwk` —
        // and under `--focus0-only` (D11) every cluster is a singleton, which
        // produced a `tree.nwk` with no leaves at all.
        let leaf = Tree::Leaf {
            name: results.names[pivot_id as usize].clone(),
            len: 0.0,
        };
        return Ok(ClusterResult {
            focus: 2,
            pivot_len: pivot_rec.stats.length,
            pivot_contig_names: results
                .contig_names
                .get(&pivot_id)
                .cloned()
                .unwrap_or_default(),
            snp_tree: leaf.clone(),
            subs_tree: leaf,
            ..Default::default()
        });
    }

    let pivot_data: GenomeData = crate::fasta::read_genome(&pivot_rec.path)
        .map_err(|e| format!("reading pivot {}: {e}", pivot_rec.path.display()))?;
    let pivot = PivotIndex::new(&pivot_rec.name, &pivot_data, params);
    if pivot.repeated_anchor_hashes > 0 {
        sink.warn(format!(
            "cluster {}: {} anchor hashes occur more than once in pivot '{}'; \
             variants near them are flagged 'repeat' rather than corrected",
            cluster.id, pivot.repeated_anchor_hashes, pivot_rec.name
        ));
    }
    let offsets = contig_offsets(&pivot);
    let pivot_concat = concat_pivot(&pivot);
    let cache = Arc::new(PieceCache::new(params.piece_cache));
    let pivot_names: Vec<String> = pivot.contigs.iter().map(|c| c.name.clone()).collect();

    // ---- Focus 1: find the differences, do not resolve them (6.11.3) ----
    let focus1: Vec<Option<MemberComparison>> = pool.install(|| {
        members
            .par_iter()
            .map(|&mid| {
                if controls.cancelled() {
                    return None;
                }
                let rec = &results.records[mid as usize];
                match crate::fasta::read_genome(&rec.path) {
                    Ok(data) => Some(compare_member(&pivot, &data, mid, params, &cache, false)),
                    Err(e) => {
                        sink.warn(format!("reading {}: {e}", rec.path.display()));
                        None
                    }
                }
            })
            .collect()
    });

    let mut counters = Counters::default();
    let mut f1: Vec<MemberComparison> = Vec::new();
    for c in focus1.into_iter().flatten() {
        counters.merge(&c.counters);
        f1.push(c);
    }
    f1.sort_by_key(|c| c.genome_id);

    let profiles1: Vec<Profile> = std::iter::once(Profile::pivot(pivot_id))
        .chain(
            f1.iter()
                .map(|c| Profile::from_focus1(c.genome_id, &c.focus1.pieces, c.uncovered.clone())),
        )
        .collect();
    let refs1: Vec<&Profile> = profiles1.iter().collect();
    let m1 = distance_matrix(&refs1, pivot.total_len, &pivot_concat);
    let names1: Vec<String> = profiles1
        .iter()
        .map(|p| results.names[p.genome_id as usize].clone())
        .collect();
    let tree1 = cluster_tree(&names1, &m1, TreeBasis::Focus1Estimate);
    sink.emit(Event::ClusterFocus(ClusterFocusInfo {
        cluster: cluster.id,
        focus: 1,
        // Focus 1 distances are counts of differing difference pieces, not SNP
        // counts, so they are approximate and must be shown as such.
        approximate: true,
        newick: to_newick(&tree1),
        members: profiles1
            .iter()
            .map(|p| MemberNumbers {
                id: p.genome_id,
                name: results.names[p.genome_id as usize].clone(),
                diff_pieces: p.pieces.len() as u64,
                snps: 0,
                indels: 0,
                uncovered_bp: crate::interval::total_len(&p.uncovered),
                novel_bp: 0,
            })
            .collect(),
    }));
    if controls.cancelled() {
        return Ok(ClusterResult {
            focus: 1,
            counters,
            pivot_len: pivot.total_len,
            matrix: m1,
            snp_tree: tree1.clone(),
            subs_tree: tree1,
            pivot_contig_names: pivot_names.clone(),
            ..Default::default()
        });
    }

    // ---- Focus 2: resolve every difference piece (6.11.3) ----
    let resolved: Vec<Result<Option<MemberComparison>, (u32, String)>> = pool.install(|| {
        f1.par_iter()
            .map(|c1| {
                if controls.cancelled() {
                    return Ok(None);
                }
                let rec = &results.records[c1.genome_id as usize];
                match resolve_focus1(&pivot, rec, &c1.focus1, params, &cache) {
                    Ok(cmp) => Ok(Some(cmp)),
                    // Carry the genome id out with the error so that the
                    // warnings and the notes can be emitted in genome order.
                    Err(e) => Err((c1.genome_id, e)),
                }
            })
            .collect()
    });

    let mut comparisons: Vec<(u32, MemberComparison)> = Vec::new();
    let mut dropped: Vec<(u32, String)> = Vec::new();
    for r in resolved {
        match r {
            Ok(Some(c)) => {
                counters.merge(&c.counters);
                comparisons.push((c.genome_id, c));
            }
            // Cancelled: skip quietly, the whole run is being wound down.
            Ok(None) => {}
            Err(e) => dropped.push(e),
        }
    }
    // Warn in genome order rather than in worker-completion order, so that
    // neither `events.jsonl` nor `run.json`'s notes depend on the thread count
    // (section 0 rule 3). The Focus 1 pass warns from inside its workers, which
    // is fine there only because a failed read removes the genome from `f1`
    // entirely; here the genome is in the tree, so the warning has to be
    // attributable.
    dropped.sort_by_key(|(id, _)| *id);
    let mut notes: Vec<String> = Vec::with_capacity(dropped.len());
    for (id, e) in dropped {
        let msg = format!(
            "cluster {}: dropping member {} from Focus 2: {e}; \
             it is reported as having no variants, which is not a measurement",
            cluster.id, results.names[id as usize]
        );
        sink.warn(msg.clone());
        notes.push(msg);
    }
    comparisons.sort_by_key(|x| x.0);
    // `cache_hits` is deliberately NOT taken from `cache.hits()`. The atomic
    // counts every successful lookup by every worker, while the merged
    // per-worker counters already record exactly one hit per piece resolved
    // from the cache. Overwriting the sum with the atomic doubled the hit count
    // and made it depend on how the workers raced.

    let mut profiles: Vec<Profile> = vec![Profile::pivot(pivot_id)];
    for (_, cmp) in &comparisons {
        profiles.push(Profile::from_comparison(cmp, &offsets));
    }
    let refs: Vec<&Profile> = profiles.iter().collect();
    let matrix = distance_matrix(&refs, pivot.total_len, &pivot_concat);
    let names: Vec<String> = profiles
        .iter()
        .map(|p| results.names[p.genome_id as usize].clone())
        .collect();

    let (snp_tree, subs_tree) = if params.profile_grouping {
        let (reps, groups) = group_identical(&refs);
        let rep_names: Vec<String> = reps.iter().map(|&i| names[i].clone()).collect();
        let sub = sub_matrix(&matrix, &reps);
        (
            expand_groups(
                cluster_tree(&rep_names, &sub, TreeBasis::Snps),
                &groups,
                &names,
            ),
            expand_groups(
                cluster_tree(&rep_names, &sub, TreeBasis::SubsPerSite),
                &groups,
                &names,
            ),
        )
    } else {
        (
            cluster_tree(&names, &matrix, TreeBasis::Snps),
            cluster_tree(&names, &matrix, TreeBasis::SubsPerSite),
        )
    };

    sink.emit(Event::ClusterFocus(ClusterFocusInfo {
        cluster: cluster.id,
        focus: 2,
        approximate: false,
        newick: to_newick(&snp_tree),
        members: std::iter::once(MemberNumbers {
            id: pivot_id,
            name: results.names[pivot_id as usize].clone(),
            diff_pieces: 0,
            snps: 0,
            indels: 0,
            uncovered_bp: 0,
            novel_bp: 0,
        })
        .chain(comparisons.iter().map(|(mid, cmp)| MemberNumbers {
            id: *mid,
            name: results.names[*mid as usize].clone(),
            diff_pieces: cmp.focus1.pieces.len() as u64,
            snps: cmp.snps(),
            indels: cmp.indels(),
            uncovered_bp: cmp.uncovered_bp(),
            novel_bp: cmp.novel_bp,
        }))
        .collect(),
    }));
    for (mid, cmp) in &comparisons {
        sink.emit(Event::MapUpdate(Box::new(MapUpdateInfo {
            cluster: cluster.id,
            genome: *mid,
            genome_name: results.names[*mid as usize].clone(),
            pivot_name: pivot_rec.name.clone(),
            pieces: cmp
                .focus1
                .pieces
                .iter()
                .map(|p| piece_view(p, &pivot_names))
                .collect(),
            variants: cmp.variants.clone(),
            uncovered: cmp.uncovered.clone(),
            pivot_contigs: pivot_names.clone(),
            pivot_contig_lens: pivot.contig_lens.clone(),
        })));
    }
    sink.emit(Event::Counters(counters.clone()));

    Ok(ClusterResult {
        focus: if controls.cancelled() { 1 } else { 2 },
        snp_tree,
        subs_tree,
        matrix,
        comparisons,
        counters,
        pivot_len: pivot.total_len,
        pivot_contig_names: pivot_names,
        notes,
    })
}

fn piece_view(p: &crate::compare::DiffPiece, pivot_names: &[String]) -> PieceView {
    PieceView {
        pivot_contig: p.pcontig,
        pivot_contig_name: pivot_names
            .get(p.pcontig as usize)
            .cloned()
            .unwrap_or_default(),
        start: p.ps,
        end: p.pe,
        resolved: true,
    }
}

/// Read a member's FASTA for the Focus 2 pass, naming the file in any error.
///
/// A member that cannot be read has to be reported, never silently dropped:
/// dropping it would leave a genome that looks identical to its pivot in
/// `genomes.tsv`, in `cluster_<id>.dist.tsv` and in the tree, which is a wrong
/// answer presented as a right one (section 0 rule 8).
fn read_member(rec: &GenomeRecord) -> Result<GenomeData, String> {
    crate::fasta::read_genome(&rec.path).map_err(|e| format!("reading {}: {e}", rec.path.display()))
}

/// Turn a Focus 1 result into an exact one without repeating steps A-D.
///
/// A member whose difference pieces are all already in the cache is not read
/// from disk again at all, which is the memory rule of 6.10.
///
/// The error is a read failure that names the file; the caller warns about it
/// and drops the member from the cluster rather than passing it off as
/// identical to the pivot.
fn resolve_focus1(
    pivot: &PivotIndex,
    rec: &GenomeRecord,
    f1: &Focus1,
    params: &Parameters,
    cache: &PieceCache,
) -> Result<MemberComparison, String> {
    // Fresh counters: the Focus 1 pass already accounted for steps A-D and the
    // caller merged those. Counting them again here would double every number
    // the experiment of section 11 reports.
    let mut counters = Counters::default();
    let dp = DiffParams::from(params);

    // `contains`, not `cached`: this pre-pass only decides whether the member's
    // FASTA has to be read, and probing must not record a hit on top of the one
    // the resolution loop records for the same piece. Counting both drove
    // `cache_hits` to 2N for N pieces in exactly the case the cache exists for
    // — a mutation shared by a whole clade, where `.any()` never short-circuits
    // and so walks every piece before the loop walks every piece again.
    let needs_data = f1.pieces.iter().any(|p| !cache.contains(&p.cache_key()));
    // Loaded eagerly when some piece is going to miss, and lazily below as a
    // fallback, so that no piece is ever skipped for want of sequence data.
    let mut data: Option<GenomeData> = if needs_data {
        Some(read_member(rec)?)
    } else {
        None
    };

    // Seed from the Focus 1 variants. The SNPs that the ungapped end extension
    // of 6.6 step D found are settled at Focus 1 and are not re-derived here,
    // so starting from an empty list would throw every one of them away —
    // 6.6 D says "mismatches inside the kept part are SNPs".
    let mut variants = f1.variants.clone();
    let mut not_called: Vec<(u32, u64, u64)> = Vec::new();
    // Cache accounting that cannot move with the thread count. Counting a hit
    // per successful lookup is racy: two workers that reach a mutation shared
    // by their two members at the same moment both miss, both resolve it, and
    // only one wins the insert, so the total depends on the schedule — measured
    // as 30 or 31 for the same 6-genome input run twice at 8 threads, and 24
    // against 31 across thread counts. What *is* schedule-independent is that
    // each distinct piece key is won exactly once, so count resolutions and
    // wins and report the difference. At one thread that is exactly the number
    // of lookups that hit; at any thread count it is the number of piece
    // resolutions the cache made redundant, which is the benefit section 11
    // reports. With the cache off, `insert` wins every time and this is 0.
    let mut resolutions = 0u64;
    let mut first_computations = 0u64;
    for p in &f1.pieces {
        let key = p.cache_key();
        resolutions += 1;
        let result = match cache.cached(&key) {
            Some((vs, ns, failed)) => DiffResult {
                variants: vs,
                n_positions: ns,
                align_failed: failed,
                used_dp: false,
            },
            None => {
                if data.is_none() {
                    data = Some(read_member(rec)?);
                }
                let d = data.as_ref().expect("just loaded");
                let alt = match d.contigs.get(p.gcontig as usize) {
                    Some(c) => p.alt(&c.seq),
                    // Cannot happen: the piece came from this same genome's
                    // Focus 1 pass, so the contig index is in range. This is a
                    // last resort, not a way to swallow a read error.
                    None => continue,
                };
                let ref_slice = &pivot.contig_seq(p.pcontig)[p.ps as usize..p.pe as usize];
                let r = crate::variant::diff(ref_slice, &alt, &dp);
                // Charge the alignment to whoever populated the slot, not to
                // every worker that missed it. These two counters therefore
                // mean "distinct pieces aligned", deliberately, so that they
                // are identical at any thread count (section 0 rule 3).
                if cache.insert(key, r.clone()) {
                    first_computations += 1;
                    counters.bases_aligned += ref_slice.len().max(alt.len()) as u64;
                    if r.used_dp {
                        counters.alignments_run += 1;
                    }
                }
                r
            }
        };
        for rv in &result.variants {
            variants.push(piece_variant(p, rv));
        }
        for &npos in &result.n_positions {
            not_called.push((p.pcontig, p.ps + npos as u64, p.ps + npos as u64 + 1));
        }
    }
    counters.cache_hits = resolutions.saturating_sub(first_computations);
    // The union of the end-extension SNPs and the resolved pieces is sorted and
    // de-duplicated, so a position found by both routes is emitted once.
    variants.sort_by(|a, b| {
        (a.pivot_contig, a.pivot_start, a.pivot_end, a.kind as u8)
            .cmp(&(b.pivot_contig, b.pivot_start, b.pivot_end, b.kind as u8))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });
    variants.dedup();

    let mut focus1 = f1.clone();
    focus1.not_called.extend(not_called);
    // `novel_bp` is recomputed here rather than inherited from the Focus 1
    // clone: at Focus 1 no insertion had been resolved yet, so its
    // inserted-bp term was zero by construction, while 6.8 counts member bases
    // that are unplaced OR inserted. Inheriting it made `genomes.tsv` omit
    // every inserted base.
    let unplaced_bp: u64 = focus1
        .unplaced
        .iter()
        .map(|&(_, s, e)| e.saturating_sub(s))
        .sum();
    let inserted_bp: u64 = variants
        .iter()
        .filter(|v| v.kind == crate::variant::VariantKind::Insertion)
        .map(|v| v.alt_allele.len() as u64)
        .sum();
    let novel_bp = unplaced_bp + inserted_bp;
    focus1.novel_bp = novel_bp;
    let cov = coverage(pivot, &focus1, &variants);

    Ok(MemberComparison {
        genome_id: rec.id,
        novel_bp,
        focus1,
        uncovered: cov.uncovered,
        duplicated: cov.duplicated,
        variants,
        counters,
    })
}

fn sub_matrix(m: &DistanceMatrix, reps: &[usize]) -> DistanceMatrix {
    let pick_u = |v: &[Vec<u64>]| -> Vec<Vec<u64>> {
        reps.iter()
            .map(|&i| reps.iter().map(|&j| v[i][j]).collect())
            .collect()
    };
    let pick_f = |v: &[Vec<f64>]| -> Vec<Vec<f64>> {
        reps.iter()
            .map(|&i| reps.iter().map(|&j| v[i][j]).collect())
            .collect()
    };
    DistanceMatrix {
        ids: reps.iter().map(|&i| m.ids[i]).collect(),
        snps: pick_u(&m.snps),
        indels: pick_u(&m.indels),
        comparable: pick_u(&m.comparable),
        dist: pick_f(&m.dist),
        estimates: pick_u(&m.estimates),
        basis: m.basis,
    }
}

fn apply_cluster_result(results: &mut RunResults, cluster: &Cluster, res: &ClusterResult) {
    results.counters.merge(&res.counters);
    // Members dropped at Focus 2 because their FASTA could not be read have to
    // show up in `run.json`, not only in the transient warning stream.
    results.notes.extend(res.notes.iter().cloned());
    results
        .cluster_trees
        .insert(cluster.id, res.subs_tree.clone());
    results
        .cluster_snp_trees
        .insert(cluster.id, res.snp_tree.clone());
    results.cluster_focus.insert(cluster.id, res.focus);
    results.matrices.insert(cluster.id, res.matrix.clone());
    results.pivot_len.insert(cluster.id, res.pivot_len);
    results
        .pivot_contig_names
        .insert(cluster.id, res.pivot_contig_names.clone());
    for (mid, cmp) in &res.comparisons {
        if let Some(r) = results.records.get_mut(*mid as usize) {
            r.focus = res.focus;
            r.snps_vs_pivot = cmp.snps();
            r.indels_vs_pivot = cmp.indels();
            r.uncovered_bp = cmp.uncovered_bp();
            r.novel_bp = cmp.novel_bp;
        }
        results.comparisons.insert(*mid, cmp.clone());
    }
    if let Some(r) = results.records.get_mut(cluster.pivot as usize) {
        r.focus = res.focus;
    }
}

/// The full tree: the backbone with every pivot leaf replaced by its cluster
/// tree, whose branch lengths are substitutions per site so that the units match
/// the backbone (6.9).
pub fn build_full_tree(results: &RunResults) -> Tree {
    let mut full = results.backbone.clone();
    let mut ids: Vec<u32> = results.clusters.iter().map(|c| c.id).collect();
    ids.sort_unstable();
    for id in ids {
        let cluster = match results.clusters.iter().find(|c| c.id == id) {
            Some(c) => c,
            None => continue,
        };
        let pivot_name = &results.names[cluster.pivot as usize];
        if let Some(sub) = results.cluster_trees.get(&id) {
            // Never graft an empty subtree. `Tree::default()` has no leaves, and
            // grafting it would replace the backbone leaf with nothing and
            // delete that genome from `tree.nwk`. A singleton cluster now gets
            // a real one-leaf tree from `refine_one_cluster`, so this is the
            // belt-and-braces half of the same fix.
            if sub.leaf_count() > 0 {
                graft(&mut full, pivot_name, sub.clone());
            }
        }
    }
    full
}

/// Flattened uncovered intervals of one comparison, for callers that need a
/// single list rather than per-contig ones.
pub fn uncovered_flat(cmp: &MemberComparison) -> Vec<Iv> {
    cmp.uncovered.clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_length_is_the_middle_value() {
        let recs: Vec<GenomeRecord> = [5u64, 1, 9, 3, 7]
            .iter()
            .enumerate()
            .map(|(i, &l)| {
                let mut r = GenomeRecord::new(
                    i as u32,
                    format!("g{i}"),
                    PathBuf::from("x"),
                    Default::default(),
                    vec![],
                );
                r.stats.length = l;
                r
            })
            .collect();
        assert_eq!(median_length(&recs), 5);
        assert_eq!(median_length(&[]), 0);
    }

    #[test]
    fn median_length_averages_the_two_middle_values_of_an_even_count() {
        // 6.1's small-genome rule turns on "under 500 kb", so a two-genome run
        // straddling the threshold has to be judged on the real median.
        let recs: Vec<GenomeRecord> = [499_998u64, 500_001]
            .iter()
            .enumerate()
            .map(|(i, &l)| {
                let mut r = GenomeRecord::new(
                    i as u32,
                    format!("g{i}"),
                    PathBuf::from("x"),
                    Default::default(),
                    vec![],
                );
                r.stats.length = l;
                r
            })
            .collect();
        assert_eq!(median_length(&recs), 499_999);
        assert!(
            Parameters::default()
                .apply_small_genome_rule(median_length(&recs))
                .is_some(),
            "the true median is under 500 kb, so the rule must fire"
        );
        // At or above the threshold the rule must stay quiet: rounding the
        // average down may not invent a small-genome run.
        assert_eq!(median_length(&recs[..1]), 499_998);
        let big: Vec<GenomeRecord> = [500_000u64, 500_001]
            .iter()
            .enumerate()
            .map(|(i, &l)| {
                let mut r = GenomeRecord::new(
                    i as u32,
                    format!("g{i}"),
                    PathBuf::from("x"),
                    Default::default(),
                    vec![],
                );
                r.stats.length = l;
                r
            })
            .collect();
        assert_eq!(median_length(&big), 500_000);
        assert!(Parameters::default()
            .apply_small_genome_rule(500_000)
            .is_none());
    }

    #[test]
    fn prioritize_puts_the_selected_cluster_first_and_loses_nothing() {
        // Cluster ids are deliberately NOT equal to their positions in the
        // vector. The previous implementation inserted the position it had
        // found instead of the cluster index it had just removed, which passes
        // unnoticed whenever ids happen to equal indices — so this test uses
        // ids 7, 3, 9 where that cannot hold.
        let clusters = vec![
            Cluster {
                id: 7,
                pivot: 0,
                members: vec![0],
            },
            Cluster {
                id: 3,
                pivot: 1,
                members: vec![1],
            },
            Cluster {
                id: 9,
                pivot: 2,
                members: vec![2],
            },
        ];
        // Work order as `refine_clusters` builds it: indices into `clusters`.
        let mut order = vec![0usize, 1, 2];

        prioritize(&mut order, &clusters, &[9]);
        assert_eq!(clusters[order[0]].id, 9, "the selected cluster jumps first");
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2], "no cluster lost or duplicated");

        // A fresh batch of requests lands at the front in the order it was
        // made, and nothing is lost.
        prioritize(&mut order, &clusters, &[9, 3]);
        assert_eq!(clusters[order[0]].id, 9);
        assert_eq!(clusters[order[1]].id, 3);
        assert_eq!(clusters[order[2]].id, 7, "the unselected cluster stays");
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2]);

        // An id that is not in this clustering changes nothing.
        let before = order.clone();
        prioritize(&mut order, &clusters, &[42]);
        assert_eq!(order, before);
        prioritize(&mut order, &clusters, &[]);
        assert_eq!(order, before);
    }

    #[test]
    fn sub_matrix_selects_the_representative_rows_and_columns() {
        let m = DistanceMatrix {
            ids: vec![0, 1, 2, 3],
            snps: vec![
                vec![0, 1, 2, 3],
                vec![1, 0, 4, 5],
                vec![2, 4, 0, 6],
                vec![3, 5, 6, 0],
            ],
            indels: vec![vec![0; 4]; 4],
            comparable: vec![vec![100; 4]; 4],
            dist: vec![vec![0.0; 4]; 4],
            estimates: vec![vec![0; 4]; 4],
            basis: TreeBasis::Snps,
        };
        let s = sub_matrix(&m, &[0, 2]);
        assert_eq!(s.ids, vec![0, 2]);
        assert_eq!(s.snps, vec![vec![0, 2], vec![2, 0]]);
    }

    #[test]
    fn summary_counts_focus_levels() {
        let r = RunResults {
            params: Parameters::default(),
            records: (0..2u32)
                .map(|i| {
                    let mut rec = GenomeRecord::new(
                        i,
                        format!("g{i}"),
                        PathBuf::from("a.fa"),
                        Default::default(),
                        vec![],
                    );
                    rec.focus = if i == 0 { 2 } else { 0 };
                    rec.stats.length = 1000;
                    rec
                })
                .collect(),
            ..Default::default()
        };
        let s = summary_of(&r);
        assert_eq!(s.focus_counts, [1, 0, 1]);
        assert_eq!(s.genomes, 2);
        assert_eq!(s.total_length, 2000);
        assert!(!s.cancelled);
    }

    #[test]
    fn controls_round_trip_a_cluster_distance() {
        let c = Controls::new(&Parameters::default());
        assert!((c.cluster_dist() - 0.005).abs() < 1e-12);
        c.cluster_dist_bits
            .store(0.01f64.to_bits(), Ordering::Relaxed);
        assert!((c.cluster_dist() - 0.01).abs() < 1e-12);
        assert!(!c.cancelled());
        c.cancel.store(true, Ordering::Relaxed);
        assert!(c.cancelled());
    }

    #[test]
    fn commands_are_drained_and_translated() {
        let c = Controls::new(&Parameters::default());
        let (tx, rx) = unbounded::<Command>();
        let sink = EventSink::null();
        tx.send(Command::Prioritize { cluster: 7 }).unwrap();
        tx.send(Command::SetClusterDist { value: 0.01 }).unwrap();
        let (prio, recluster) = c.poll(&rx, &sink);
        assert_eq!(prio, vec![7]);
        assert!(recluster);
        assert!((c.cluster_dist() - 0.01).abs() < 1e-12);
        // An out-of-range threshold must be refused, not applied.
        tx.send(Command::SetClusterDist { value: 5.0 }).unwrap();
        let (_, recluster2) = c.poll(&rx, &sink);
        assert!(!recluster2);
        assert!((c.cluster_dist() - 0.01).abs() < 1e-12);
        tx.send(Command::Cancel).unwrap();
        c.poll(&rx, &sink);
        assert!(c.cancelled());
    }

    #[test]
    fn an_empty_input_list_is_an_error_not_an_empty_run() {
        // A run over nothing must fail loudly: an empty output folder would look
        // like a successful run with no differences.
        let dir = std::env::temp_dir().join(format!("mtato-empty-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let cfg = RunConfig {
            inputs: vec![dir.clone()],
            params: Parameters::default(),
            outdir: None,
        };
        let (sink, _rx) = event_bus();
        let (_tx, cmd_rx) = unbounded::<Command>();
        let err = run(cfg, sink, cmd_rx).unwrap_err();
        assert!(err.contains("no genome files"), "{err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_missing_input_path_is_reported_with_the_path() {
        let cfg = RunConfig {
            inputs: vec![PathBuf::from("/nonexistent/mtato-test-dir")],
            params: Parameters::default(),
            outdir: None,
        };
        let (sink, _rx) = event_bus();
        let (_tx, cmd_rx) = unbounded::<Command>();
        let err = run(cfg, sink, cmd_rx).unwrap_err();
        assert!(err.contains("/nonexistent/mtato-test-dir"), "{err}");
    }

    #[test]
    fn invalid_parameters_are_rejected_before_any_work() {
        let cfg = RunConfig {
            inputs: vec![PathBuf::from(".")],
            params: Parameters {
                k: 30,
                ..Default::default()
            },
            outdir: None,
        };
        let (sink, _rx) = event_bus();
        let (_tx, cmd_rx) = unbounded::<Command>();
        assert!(run(cfg, sink, cmd_rx).is_err());
    }

    /// A pivot, a member that differs from it in the middle, and a Focus 1
    /// comparison of the two, for the `resolve_focus1` tests below.
    fn focus1_with_pieces() -> (PivotIndex, GenomeData, Focus1, Parameters) {
        let params = Parameters {
            k: 21,
            densities: [64, 16, 4],
            ..Default::default()
        };
        // Deterministic pseudo-random and N-free. The k-mers have to be unique
        // enough for the pivot's anchor map to keep them — it drops repeated
        // hashes (6.5) — otherwise no shared point, and so no piece, would ever
        // arise and the tests below would prove nothing.
        let mut x = 0x2545_F491_4F6C_DD1Du64;
        let mut base: Vec<u8> = Vec::with_capacity(4000);
        for _ in 0..4000 {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            base.push(b"ACGT"[((x >> 33) % 4) as usize]);
        }
        let pivot_data = GenomeData {
            name: "pivot".into(),
            contigs: vec![crate::fasta::Contig {
                name: "c0".into(),
                seq: base.clone(),
            }],
        };
        let pivot = PivotIndex::new("pivot", &pivot_data, &params);
        // The member differs over a stretch in the middle, so that stretch becomes
        // a difference piece rather than a matched block.
        let mut member = base;
        for b in &mut member[1800..1860] {
            *b = b"ACGT"[(crate::seq::code(*b) + 1) as usize % 4];
        }
        let member_data = GenomeData {
            name: "member".into(),
            contigs: vec![crate::fasta::Contig {
                name: "c0".into(),
                seq: member,
            }],
        };
        let cache = PieceCache::new(true);
        let cmp = compare_member(&pivot, &member_data, 1, &params, &cache, false);
        (pivot, member_data, cmp.focus1, params)
    }

    #[test]
    fn a_member_that_cannot_be_read_at_focus_2_is_an_error_naming_the_file() {
        // The Focus 2 pass used to swallow the read error with `.ok()` and then
        // skip every piece, which produced a member with zero variants and no
        // warning at all — indistinguishable from a genome identical to its
        // pivot in `genomes.tsv`, in `cluster_<id>.dist.tsv` and in the tree.
        let (pivot, _member, f1, params) = focus1_with_pieces();
        assert!(
            !f1.pieces.is_empty(),
            "the test needs at least one difference piece"
        );
        let missing = GenomeRecord::new(
            1,
            "member".into(),
            PathBuf::from("/nonexistent/mtato-focus2-member.fasta"),
            Default::default(),
            vec![],
        );
        let cache = PieceCache::new(true);
        let err = resolve_focus1(&pivot, &missing, &f1, &params, &cache).unwrap_err();
        assert!(err.contains("reading"), "{err}");
        assert!(
            err.contains("/nonexistent/mtato-focus2-member.fasta"),
            "the message must name the file: {err}"
        );
    }

    #[test]
    fn a_fully_cached_member_is_not_read_from_disk_at_focus_2() {
        // 6.10: "a member whose difference pieces are all in the piece cache
        // does not need to be re-read for Focus 2". Warm the cache, then point
        // the record at a file that does not exist: the resolution must still
        // succeed, which is only possible if no read was attempted.
        let (pivot, member_data, f1, params) = focus1_with_pieces();
        let cache = PieceCache::new(true);
        let dp = DiffParams::from(&params);
        let seq = &member_data.contigs[0].seq;
        for p in &f1.pieces {
            let alt = p.alt(seq);
            let ref_slice = &pivot.contig_seq(p.pcontig)[p.ps as usize..p.pe as usize];
            cache.insert(p.cache_key(), crate::variant::diff(ref_slice, &alt, &dp));
        }
        assert_eq!(cache.len(), f1.pieces.len(), "every piece is warm");

        let gone = GenomeRecord::new(
            1,
            "member".into(),
            PathBuf::from("/nonexistent/mtato-focus2-cached.fasta"),
            Default::default(),
            vec![],
        );
        let cmp = resolve_focus1(&pivot, &gone, &f1, &params, &cache)
            .expect("a fully cached member needs no read");
        assert_eq!(
            cmp.counters.cache_hits as usize,
            f1.pieces.len(),
            "every piece was served from the cache"
        );
        assert_eq!(cmp.counters.bases_aligned, 0, "nothing was aligned");
        assert_eq!(cmp.counters.alignments_run, 0);
        assert!(
            !cmp.variants.is_empty(),
            "the cached variants still have to be reported"
        );
    }

    #[test]
    fn cache_hits_are_derived_so_they_do_not_depend_on_the_schedule() {
        // Half of one member's pieces are already warm. The hit count must be
        // exactly the number of resolutions that did not have to be the first
        // computation of their key, and the alignment work must be charged only
        // for the ones that were — which is what makes both numbers identical
        // at any thread count.
        let (pivot, member_data, f1, params) = focus1_with_pieces();
        let dir = std::env::temp_dir().join(format!("mtato-focus2-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("member.fasta");
        {
            use std::io::Write;
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, ">c0").unwrap();
            let seq = &member_data.contigs[0].seq;
            for chunk in seq.chunks(60) {
                writeln!(f, "{}", std::str::from_utf8(chunk).unwrap()).unwrap();
            }
        }
        let rec = GenomeRecord::new(1, "member".into(), path, Default::default(), vec![]);

        let cache = PieceCache::new(true);
        let dp = DiffParams::from(&params);
        let seq = &member_data.contigs[0].seq;
        let half = f1.pieces.len() / 2;
        for p in &f1.pieces[..half] {
            let alt = p.alt(seq);
            let r = &pivot.contig_seq(p.pcontig)[p.ps as usize..p.pe as usize];
            cache.insert(p.cache_key(), crate::variant::diff(r, &alt, &dp));
        }
        // What a strictly sequential resolution would hit: a piece hits when its
        // key was already warm or was computed earlier in this same walk.
        let mut seen: HashSet<_> = f1.pieces[..half].iter().map(|p| p.cache_key()).collect();
        let mut expected = 0u64;
        for p in &f1.pieces {
            if !seen.insert(p.cache_key()) {
                expected += 1;
            }
        }

        let cmp = resolve_focus1(&pivot, &rec, &f1, &params, &cache).unwrap();
        assert_eq!(
            cmp.counters.cache_hits, expected,
            "hits are the redundancies"
        );
        assert!(
            cmp.counters.bases_aligned > 0,
            "the cold pieces were aligned"
        );

        // Everything is warm now, so a second pass hits on every piece, aligns
        // nothing, and must give the same answer.
        let again = resolve_focus1(&pivot, &rec, &f1, &params, &cache).unwrap();
        assert_eq!(again.counters.cache_hits, f1.pieces.len() as u64);
        assert_eq!(again.counters.bases_aligned, 0);
        assert_eq!(again.counters.alignments_run, 0);
        assert_eq!(again.variants, cmp.variants, "the answer does not change");
        std::fs::remove_dir_all(&dir).ok();
    }
}
