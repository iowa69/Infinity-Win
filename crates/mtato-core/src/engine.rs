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
    let (records, contig_names, mut notes, scan_ms) =
        scan_all(&pool, &files, &params, &sink, &controls)?;
    results.records = records;
    results.contig_names = contig_names;
    let mut scan_ms_total = scan_ms;

    // The small-genome rule of 6.1 can only be judged once the lengths are
    // known, so it is applied after a first scan and the scan is repeated with
    // the denser sampling. For bacterial assemblies, which is the design
    // target, this second pass never runs.
    let median = median_length(&results.records);
    let mut effective = params.clone();
    if let Some(msg) = effective.apply_small_genome_rule(median) {
        sink.warn(msg.clone());
        notes.push(msg);
        let (records, contig_names, _, rescan_ms) =
            scan_all(&pool, &files, &effective, &sink, &controls)?;
        results.records = records;
        results.contig_names = contig_names;
        scan_ms_total = scan_ms_total.saturating_add(rescan_ms);
        params = effective.clone();
        results.params = params.clone();
    }
    results.notes = notes;
    results.timings.scan_ms = Some(scan_ms_total);
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
        let preview = summary_of(&results);
        crate::output::write_all(dir, &results, &preview)
            .map_err(|e| format!("writing output to {}: {e}", dir.display()))?;
        results.timings.output_ms = Some(t.elapsed().as_millis() as u64);
        sink.progress(stage::OUTPUT, 1, 1);
    }

    Ok(finish(results, &sink, started))
}

fn finish(mut results: RunResults, sink: &EventSink, started: Instant) -> RunResults {
    results.timings.total_ms = Some(started.elapsed().as_millis() as u64);
    let summary = summary_of(&results);
    sink.emit(Event::Counters(results.counters.clone()));
    sink.emit(Event::Finished(summary));
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
    lens[lens.len() / 2]
}

// ---------------------------------------------------------------------------
// Scan
// ---------------------------------------------------------------------------

fn scan_all(
    pool: &ThreadPool,
    files: &[PathBuf],
    params: &Parameters,
    sink: &EventSink,
    controls: &Controls,
) -> Result<
    (
        Vec<GenomeRecord>,
        HashMap<u32, Vec<String>>,
        Vec<String>,
        u64,
    ),
    String,
> {
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
    Ok((records, contig_names, notes, t.elapsed().as_millis() as u64))
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
        for want in prioritised.iter().rev() {
            if let Some(pos) = order
                .iter()
                .position(|&i| clustering.clusters[i].id == *want)
            {
                order.remove(pos);
                order.insert(0, pos);
            }
        }

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
        // "Clusters with a single genome need no index" (6.5).
        return Ok(ClusterResult {
            focus: 2,
            pivot_len: pivot_rec.stats.length,
            pivot_contig_names: results
                .contig_names
                .get(&pivot_id)
                .cloned()
                .unwrap_or_default(),
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
    let resolved: Vec<Option<MemberComparison>> = pool.install(|| {
        f1.par_iter()
            .map(|c1| {
                if controls.cancelled() {
                    return None;
                }
                let rec = &results.records[c1.genome_id as usize];
                Some(resolve_focus1(&pivot, rec, &c1.focus1, params, &cache))
            })
            .collect()
    });

    let mut comparisons: Vec<(u32, MemberComparison)> = Vec::new();
    for c in resolved.into_iter().flatten() {
        counters.merge(&c.counters);
        comparisons.push((c.genome_id, c));
    }
    comparisons.sort_by_key(|x| x.0);
    counters.cache_hits = cache.hits();

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

/// Turn a Focus 1 result into an exact one without repeating steps A-D.
///
/// A member whose difference pieces are all already in the cache is not read
/// from disk again at all, which is the memory rule of 6.10.
fn resolve_focus1(
    pivot: &PivotIndex,
    rec: &GenomeRecord,
    f1: &Focus1,
    params: &Parameters,
    cache: &PieceCache,
) -> MemberComparison {
    // Fresh counters: the Focus 1 pass already accounted for steps A-D and the
    // caller merged those. Counting them again here would double every number
    // the experiment of section 11 reports.
    let mut counters = Counters::default();
    let dp = DiffParams::from(params);

    let needs_data = f1
        .pieces
        .iter()
        .any(|p| cache.cached(&p.cache_key()).is_none());
    let data = if needs_data {
        match crate::fasta::read_genome(&rec.path) {
            Ok(d) => Some(d),
            Err(_) => None,
        }
    } else {
        None
    };

    let mut variants = Vec::new();
    let mut not_called: Vec<(u32, u64, u64)> = Vec::new();
    for p in &f1.pieces {
        let key = p.cache_key();
        let result = match cache.cached(&key) {
            Some((vs, ns, failed)) => {
                counters.cache_hits += 1;
                DiffResult {
                    variants: vs,
                    n_positions: ns,
                    align_failed: failed,
                    used_dp: false,
                }
            }
            None => {
                let alt = match &data {
                    Some(d) => match d.contigs.get(p.gcontig as usize) {
                        Some(c) => p.alt(&c.seq),
                        None => continue,
                    },
                    None => continue,
                };
                let ref_slice = &pivot.contig_seq(p.pcontig)[p.ps as usize..p.pe as usize];
                counters.bases_aligned += ref_slice.len().max(alt.len()) as u64;
                let r = crate::variant::diff(ref_slice, &alt, &dp);
                if r.used_dp {
                    counters.alignments_run += 1;
                }
                cache.insert(key, r.clone());
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
    variants.sort_by(|a, b| {
        (a.pivot_contig, a.pivot_start, a.pivot_end, a.kind as u8)
            .cmp(&(b.pivot_contig, b.pivot_start, b.pivot_end, b.kind as u8))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });
    variants.dedup();

    let mut focus1 = f1.clone();
    focus1.not_called.extend(not_called);
    let novel_bp = focus1.novel_bp;
    let cov = coverage(pivot, &focus1, &variants);

    MemberComparison {
        genome_id: rec.id,
        novel_bp,
        focus1,
        uncovered: cov.uncovered,
        duplicated: cov.duplicated,
        variants,
        counters,
    }
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
            graft(&mut full, pivot_name, sub.clone());
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
        let mut r = RunResults::default();
        r.params = Parameters::default();
        r.records = (0..2u32)
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
            .collect();
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
}
