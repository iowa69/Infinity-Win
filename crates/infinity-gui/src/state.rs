//! The UI model: a pure fold of the engine's event stream (section 7) into the
//! data every panel draws.
//!
//! This module knows nothing about egui. [`UiModel::apply`] is the only place an
//! [`Event`] is interpreted, which is what makes the whole desktop program
//! testable: a test drives a sequence of events through the fold and inspects
//! the model, exactly as the window would.
//!
//! # Two ways the stream is defended against
//!
//! 1. **Genome ids.** `scan_all` numbers `GenomeScanned` events by *input
//!    index*, then renumbers the surviving records densely, so when one file
//!    fails to read the ids in the events and the ids in every later event
//!    disagree. Every id from the stream is therefore resolved through the
//!    genome *name* first and only falls back to the raw id.
//! 2. **Focus regressions.** After a `SetClusterDist` command the engine emits a
//!    fresh Focus 0 tree for *every* cluster, including the ones whose members
//!    did not change and whose cached results it never re-emits. Accepting that
//!    would show an exact tree as blurry again, so a Focus 0 event is ignored
//!    when the cluster's member set is unchanged and a sharper tree is already
//!    on screen. Both are reported in the run notes rather than hidden.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use infinity_core::counters::Counters;
use infinity_core::events::{
    ClusterDetailInfo, ClusterFocusInfo, ClusterInfo, Event, GenomeScannedInfo, MapUpdateInfo,
    MemberNumbers, OutlierInfo,
};
use infinity_core::summary::RunSummary;
use infinity_core::tree::{parse_newick, Tree};

use crate::fmt;
use crate::layout::{self, DisplayNode, NodeKind};

/// One row of the genome list.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GenomeRow {
    pub id: u32,
    pub name: String,
    pub length: u64,
    pub contigs: u32,
    pub n50: u64,
    pub gc: f64,
    pub sketch_size: usize,
    /// `None` until clustering has happened.
    pub cluster: Option<u32>,
    pub is_pivot: bool,
    pub focus: u8,
    /// True while the numbers next to this genome are estimates.
    pub approximate: bool,
    pub numbers: Option<MemberNumbers>,
}

/// One cluster as the interface sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterView {
    pub id: u32,
    pub pivot: u32,
    pub pivot_name: String,
    /// Genome ids, in the model's own numbering.
    pub members: Vec<u32>,
    /// Names exactly as the engine labels them in trees and files.
    pub member_names: Vec<String>,
    pub focus: u8,
    pub approximate: bool,
    /// The latest tree for this cluster, parsed from the event's Newick.
    pub tree: Option<Tree>,
    pub numbers: HashMap<u32, MemberNumbers>,
    /// Member ids, kept so that a Focus 0 event can be recognised as a
    /// regression rather than as a genuine re-clustering.
    pub signature: Vec<u32>,
}

impl ClusterView {
    pub fn size(&self) -> usize {
        self.members.len()
    }
}

/// One segment of the status bar.
#[derive(Debug, Clone, PartialEq)]
pub struct StageState {
    pub name: String,
    pub done: u64,
    pub total: u64,
    /// Engine clock when the stage was first seen.
    pub started_ms: u64,
    /// Engine clock of the latest progress event for this stage.
    pub last_ms: u64,
    pub finished: bool,
}

impl StageState {
    /// How long the stage took, or has taken so far.
    pub fn elapsed_ms(&self) -> u64 {
        self.last_ms.saturating_sub(self.started_ms)
    }

    /// Fraction complete in `[0, 1]`, or `None` when the stage never reported a
    /// total. A stage without a total is not "0 %" and not "100 %".
    pub fn fraction(&self) -> Option<f64> {
        if self.total == 0 {
            return None;
        }
        Some((self.done as f64 / self.total as f64).clamp(0.0, 1.0))
    }
}

/// A warning from the engine, or from the interface itself.
#[derive(Debug, Clone, PartialEq)]
pub struct Note {
    pub t_ms: f64,
    pub text: String,
    /// True when the interface added it rather than the engine.
    pub local: bool,
}

/// Where the model's data came from.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Source {
    #[default]
    Empty,
    /// A run that is going now, or one that just finished in this window.
    Live,
    /// A finished output folder, opened without recomputing (section 8).
    Folder(PathBuf),
}

/// The unit of the backbone tree's branch lengths.
///
/// Today's `BackboneTree` event carries sketch distances in substitutions per
/// site, which the tree panel turns into "≈ SNPs" by multiplying with the median
/// genome length. The engine is about to emit a backbone in SNP units of its own
/// (a separate event); when it does, the fold sets [`BackboneUnit::Snps`] and
/// the conversion disappears.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BackboneUnit {
    /// Substitutions per site, from the sketches (Focus 0).
    #[default]
    SubsPerSite,
    /// SNPs already. `exact` is false while they are still estimates.
    Snps { exact: bool },
}

/// Everything the panels draw, folded out of the event stream.
#[derive(Debug, Clone, Default)]
pub struct UiModel {
    pub genomes: Vec<GenomeRow>,
    /// Genome id to position in `genomes`.
    index: HashMap<u32, usize>,
    /// Name to genome id. Only trustworthy when the name is not in `dup_names`.
    name_index: HashMap<String, u32>,
    /// Names that two genomes share, so a name cannot identify a genome.
    dup_names: HashSet<String>,

    pub clusters: Vec<ClusterView>,
    /// Cluster id to position in `clusters`.
    cluster_index: HashMap<u32, usize>,
    /// Pivot name to cluster id, for grafting clusters into the backbone.
    pivot_to_cluster: HashMap<String, u32>,

    pub backbone: Option<Tree>,
    pub backbone_pivots: Vec<String>,
    pub counters: Counters,
    pub stages: Vec<StageState>,
    pub notes: Vec<Note>,
    pub summary: Option<RunSummary>,
    /// Cluster id to genome id to map data.
    pub maps: HashMap<u32, HashMap<u32, MapUpdateInfo>>,

    /// What the backbone's branch lengths are measured in.
    pub backbone_unit: BackboneUnit,

    /// Genomes the engine flagged as divergent (another species,
    /// contamination, an incomplete assembly), in the engine's order.
    pub outliers: Vec<OutlierInfo>,
    /// Cluster id to the engine's distance matrix for it (`ClusterDetail`),
    /// which the multi-genome selection summary and the SNP groups read.
    pub details: HashMap<u32, ClusterDetailInfo>,
    /// Bumped whenever `details` changes, so the SNP groups know to
    /// recompute without the tree having to relayout.
    pub detail_version: u64,

    /// Cluster distance threshold in force, as last reported by the engine.
    pub threshold: f64,
    pub events_seen: u64,
    pub highest_seq: u64,
    pub finished: bool,
    /// Bumped whenever the *shape* of the tree could have changed, so that the
    /// tree panel knows when to relayout and start the move animation instead of
    /// rebuilding both every frame.
    pub version: u64,
    pub source: Source,
}

impl UiModel {
    pub fn new(source: Source) -> Self {
        Self {
            source,
            ..Default::default()
        }
    }

    // -- lookups ------------------------------------------------------------

    pub fn genome(&self, id: u32) -> Option<&GenomeRow> {
        self.index.get(&id).map(|&i| &self.genomes[i])
    }

    pub fn cluster(&self, id: u32) -> Option<&ClusterView> {
        self.cluster_index.get(&id).map(|&i| &self.clusters[i])
    }

    pub fn map(&self, cluster: u32, genome: u32) -> Option<&MapUpdateInfo> {
        self.maps.get(&cluster)?.get(&genome)
    }

    pub fn genome_count(&self) -> usize {
        self.genomes.len()
    }

    pub fn cluster_count(&self) -> usize {
        self.clusters.len()
    }

    pub fn total_length(&self) -> u64 {
        self.genomes.iter().map(|g| g.length).sum()
    }

    /// Median genome length, the factor that turns a sketch distance
    /// (substitutions per site) into an estimated SNP count. 0 before any
    /// genome is scanned.
    pub fn median_length(&self) -> u64 {
        let mut v: Vec<u64> = self.genomes.iter().map(|g| g.length).collect();
        if v.is_empty() {
            return 0;
        }
        let mid = v.len() / 2;
        let (_, m, _) = v.select_nth_unstable(mid);
        *m
    }

    /// The factor from the backbone's own unit to SNPs, and whether the result
    /// is an estimate.
    pub fn backbone_to_snps(&self) -> (f64, bool) {
        match self.backbone_unit {
            BackboneUnit::SubsPerSite => (self.median_length() as f64, true),
            BackboneUnit::Snps { exact } => (1.0, !exact),
        }
    }

    /// Clusters that have been refined at least to Focus 1.
    pub fn refined_clusters(&self) -> usize {
        self.clusters.iter().filter(|c| c.focus >= 1).count()
    }

    /// Clusters that reached Focus 2.
    pub fn exact_clusters(&self) -> usize {
        self.clusters.iter().filter(|c| c.focus >= 2).count()
    }

    /// Name to genome id, for turning tree leaf labels back into genomes. Names
    /// two genomes share are left out: guessing which of the two a label means
    /// would be worse than not linking the leaf at all.
    pub fn tree_name_index(&self) -> HashMap<String, u32> {
        self.name_index
            .iter()
            .filter(|(n, _)| !self.dup_names.contains(*n))
            .map(|(n, &id)| (n.clone(), id))
            .collect()
    }

    /// The first genome whose name contains `query`, case-insensitively. This is
    /// the tree panel's search box (section 8).
    pub fn search(&self, query: &str) -> Option<u32> {
        let q = query.trim().to_lowercase();
        if q.is_empty() {
            return None;
        }
        self.genomes
            .iter()
            .find(|g| g.name.to_lowercase().contains(&q))
            .map(|g| g.id)
    }

    // -- the fold -----------------------------------------------------------

    /// Apply one whole envelope, remembering its sequence number so that a
    /// dropped message can be told from a late one.
    pub fn apply_envelope(&mut self, env: infinity_core::events::Envelope) {
        self.highest_seq = self.highest_seq.max(env.seq);
        let t_ms = env.t_ms;
        self.apply(t_ms, env.event);
    }

    /// Apply one event, timestamped `t_ms` milliseconds into the run.
    pub fn apply(&mut self, t_ms: f64, event: Event) {
        self.events_seen += 1;
        match event {
            Event::GenomeScanned(g) => self.on_scanned(g),
            Event::Clustered {
                threshold,
                clusters,
            } => self.on_clustered(t_ms, threshold, clusters),
            Event::BackboneTree { newick, pivots } => {
                self.on_backbone(t_ms, newick, pivots);
                self.backbone_unit = BackboneUnit::SubsPerSite;
            }
            Event::ClusterFocus(info) => self.on_cluster_focus(t_ms, info),
            Event::MapUpdate(info) => self.on_map_update(*info),
            Event::Counters(c) => self.counters.merge(&c),
            Event::Progress {
                stage,
                done,
                total,
                elapsed_ms,
            } => self.on_progress(&stage, done, total, elapsed_ms),
            Event::Warning { text } => self.notes.push(Note {
                t_ms,
                text,
                local: false,
            }),
            Event::Finished(s) => self.on_finished(*s),
            // The exact backbone: pivot-to-pivot SNPs replace the sketch tree.
            Event::BackboneExact(b) => {
                // Pairs too far apart to align carry a sketch estimate and
                // nothing compared; the tree is exact only if no pair is.
                let n = b.pivots.len();
                let exact = (0..n).all(|i| (0..n).all(|j| i == j || b.comparable[i][j] > 0));
                self.on_backbone(t_ms, b.newick, b.pivots);
                self.backbone_unit = BackboneUnit::Snps { exact };
            }
            Event::Outliers(list) => {
                self.notes.push(Note {
                    t_ms,
                    text: format!(
                        "{} divergent assembl{} flagged (another species, contamination or an \
                         incomplete assembly): {}",
                        list.len(),
                        if list.len() == 1 { "y" } else { "ies" },
                        list.iter()
                            .map(|o| o.name.as_str())
                            .collect::<Vec<_>>()
                            .join(", ")
                    ),
                    local: false,
                });
                self.outliers = list;
            }
            // A later, sharper matrix replaces an earlier one.
            Event::ClusterDetail(d) => {
                let keep = self
                    .details
                    .get(&d.cluster)
                    .map_or(true, |old| old.focus <= d.focus || old.names != d.names);
                if keep {
                    self.details.insert(d.cluster, *d);
                    self.detail_version += 1;
                }
            }
            // Focus 1 maps and per-variant filters are written to the output
            // folder; the panels do not draw them yet.
            Event::Focus1Map(_) | Event::MapFilters { .. } => {}
        }
    }

    fn on_scanned(&mut self, g: GenomeScannedInfo) {
        let id = g.id;
        match self.index.get(&id).copied() {
            Some(i) => {
                // The small-genome rule of 6.1 rescans, so the same id arrives
                // twice; the second scan is the one that counts.
                let row = &mut self.genomes[i];
                row.name = g.name.clone();
                row.length = g.length;
                row.contigs = g.contigs;
                row.n50 = g.n50;
                row.gc = g.gc;
                row.sketch_size = g.sketch_size;
            }
            None => {
                self.genomes.push(GenomeRow {
                    id,
                    name: g.name.clone(),
                    length: g.length,
                    contigs: g.contigs,
                    n50: g.n50,
                    gc: g.gc,
                    sketch_size: g.sketch_size,
                    cluster: None,
                    is_pivot: false,
                    focus: fmt::FOCUS_SKETCH,
                    approximate: true,
                    numbers: None,
                });
                self.index.insert(id, self.genomes.len() - 1);
            }
        }
        self.note_name(&g.name, id, false);
    }

    /// Record a name to id binding. `authoritative` is true when the name comes
    /// from the engine's own de-duplicated list, which overrides a guess.
    fn note_name(&mut self, name: &str, id: u32, authoritative: bool) {
        if authoritative {
            self.name_index.insert(name.to_string(), id);
            self.dup_names.remove(name);
            return;
        }
        match self.name_index.get(name) {
            Some(&old) if old != id => {
                self.dup_names.insert(name.to_string());
            }
            Some(_) => {}
            None => {
                self.name_index.insert(name.to_string(), id);
            }
        }
    }

    /// Map an id from the stream onto this model's genome id, preferring the
    /// name (see the module documentation).
    fn resolve(&self, name: &str, id: u32) -> Option<u32> {
        if !self.dup_names.contains(name) {
            if let Some(&gid) = self.name_index.get(name) {
                if self.index.contains_key(&gid) {
                    return Some(gid);
                }
            }
        }
        if self.index.contains_key(&id) {
            return Some(id);
        }
        None
    }

    fn on_clustered(&mut self, t_ms: f64, threshold: f64, clusters: Vec<ClusterInfo>) {
        self.version += 1;
        self.threshold = threshold;
        // Results keyed by membership, exactly like the engine's own cache: a
        // cluster that survives a re-clustering keeps its tree and its focus.
        let mut old_by_sig: HashMap<Vec<u32>, ClusterView> = HashMap::new();
        for c in self.clusters.drain(..) {
            old_by_sig.insert(c.signature.clone(), c);
        }
        self.cluster_index.clear();
        self.pivot_to_cluster.clear();
        for g in &mut self.genomes {
            g.cluster = None;
            g.is_pivot = false;
            g.focus = fmt::FOCUS_SKETCH;
            g.approximate = true;
        }

        let mut new_clusters = 0usize;
        for info in clusters {
            let mut members: Vec<u32> = Vec::with_capacity(info.members.len());
            let mut names: Vec<String> = Vec::with_capacity(info.members.len());
            for (i, &raw) in info.members.iter().enumerate() {
                let name = info
                    .member_names
                    .get(i)
                    .cloned()
                    .unwrap_or_else(|| format!("genome {raw}"));
                if let Some(gid) = self.resolve(&name, raw) {
                    self.note_name(&name, gid, true);
                    members.push(gid);
                    names.push(name);
                }
            }
            let pivot_name = info.pivot_name.clone();
            let pivot = self.resolve(&pivot_name, info.pivot).unwrap_or(info.pivot);
            let signature = members.clone();

            let mut view = match old_by_sig.get(&signature) {
                // Same members, so the engine serves this cluster from its own
                // cache and never re-emits its result: keep what we have.
                Some(old) => ClusterView {
                    id: info.id,
                    pivot,
                    pivot_name: pivot_name.clone(),
                    member_names: names,
                    members: signature.clone(),
                    focus: old.focus,
                    approximate: old.approximate,
                    tree: old.tree.clone(),
                    numbers: old.numbers.clone(),
                    signature,
                },
                None => {
                    if !old_by_sig.is_empty() {
                        new_clusters += 1;
                    }
                    ClusterView {
                        id: info.id,
                        pivot,
                        pivot_name: pivot_name.clone(),
                        members: signature.clone(),
                        member_names: names,
                        focus: fmt::FOCUS_SKETCH,
                        approximate: true,
                        tree: None,
                        numbers: HashMap::new(),
                        signature,
                    }
                }
            };
            if view.members.len() <= 1 {
                // Nothing to compare, so the engine reports it as exact; see
                // `refine_one_cluster`, which returns before emitting anything.
                view.focus = fmt::FOCUS_EXACT;
                view.approximate = false;
            }
            self.pivot_to_cluster.insert(pivot_name, info.id);
            self.cluster_index.insert(info.id, self.clusters.len());
            self.clusters.push(view);
        }

        // Assign cluster membership now that the views exist.
        for ci in 0..self.clusters.len() {
            let id = self.clusters[ci].id;
            let pivot = self.clusters[ci].pivot;
            let focus = self.clusters[ci].focus;
            let approx = self.clusters[ci].approximate;
            let members = self.clusters[ci].members.clone();
            for &m in &members {
                if let Some(&gi) = self.index.get(&m) {
                    self.genomes[gi].cluster = Some(id);
                    self.genomes[gi].is_pivot = m == pivot;
                    self.genomes[gi].focus = focus;
                    self.genomes[gi].approximate = approx;
                }
            }
        }
        // A matrix belongs to a membership, not to a cluster number: after a
        // re-clustering only the ones whose cluster kept its members stay.
        let clusters = &self.clusters;
        let index = &self.cluster_index;
        self.detail_version += 1;
        self.details.retain(|cid, d| {
            index.get(cid).is_some_and(|&ci| {
                let mut a: Vec<&str> = clusters[ci]
                    .member_names
                    .iter()
                    .map(String::as_str)
                    .collect();
                let mut b: Vec<&str> = d.names.iter().map(String::as_str).collect();
                a.sort_unstable();
                b.sort_unstable();
                a == b
            })
        });
        if new_clusters > 0 {
            self.notes.push(Note {
                t_ms,
                text: format!(
                    "re-clustered at d<={threshold}: {new_clusters} clusters changed membership \
                     and start again from Focus 0"
                ),
                local: true,
            });
        }
    }

    fn on_backbone(&mut self, t_ms: f64, newick: String, pivots: Vec<String>) {
        match parse_newick(&newick) {
            Ok(t) => {
                self.version += 1;
                self.backbone = Some(t);
                self.backbone_pivots = pivots;
            }
            Err(e) => self.notes.push(Note {
                t_ms,
                text: format!("cannot read the backbone tree: {e}"),
                local: true,
            }),
        }
    }

    fn on_cluster_focus(&mut self, t_ms: f64, info: ClusterFocusInfo) {
        let Some(&ci) = self.cluster_index.get(&info.cluster) else {
            return;
        };
        if info.focus == 0 && self.clusters[ci].focus > 0 {
            // The engine re-emits a Focus 0 star for every cluster after a
            // re-clustering, including the ones it then serves from its cache
            // and never reports again. Trusting it would show an exact tree as
            // blurry. A membership change goes through `Clustered` first, which
            // has already reset this cluster to Focus 0.
            return;
        }
        let tree = match parse_newick(&info.newick) {
            Ok(t) => Some(t),
            Err(e) => {
                self.notes.push(Note {
                    t_ms,
                    text: format!("cannot read the tree of cluster {}: {e}", info.cluster),
                    local: true,
                });
                None
            }
        };
        let mut numbers = HashMap::with_capacity(info.members.len());
        for m in &info.members {
            if let Some(gid) = self.resolve(&m.name, m.id) {
                numbers.insert(gid, m.clone());
            }
        }
        let focus = info.focus;
        let approx = info.approximate || focus < fmt::FOCUS_EXACT;
        {
            let c = &mut self.clusters[ci];
            c.focus = focus;
            c.approximate = approx;
            if tree.is_some() {
                c.tree = tree;
            }
            if !numbers.is_empty() {
                c.numbers = numbers;
            }
        }
        self.version += 1;
        let members = self.clusters[ci].members.clone();
        for &gid in &members {
            if let Some(&gi) = self.index.get(&gid) {
                self.genomes[gi].focus = focus;
                self.genomes[gi].approximate = approx;
                if let Some(n) = self.clusters[ci].numbers.get(&gid) {
                    self.genomes[gi].numbers = Some(n.clone());
                }
            }
        }
    }

    fn on_map_update(&mut self, m: MapUpdateInfo) {
        let gid = match self.resolve(&m.genome_name, m.genome) {
            Some(g) => g,
            None => return,
        };
        let cluster = m.cluster;
        self.maps.entry(cluster).or_default().insert(gid, m);
    }

    fn on_progress(&mut self, stage: &str, done: u64, total: u64, elapsed_ms: u64) {
        let done = if stage == infinity_core::events::stage::REFINE {
            // After a re-clustering the engine restarts its refine counter from
            // zero and never counts the clusters it serves from its cache, so
            // its number can be lower than the truth on screen.
            done.max(self.refined_clusters() as u64)
        } else {
            done
        };
        match self.stages.iter_mut().find(|s| s.name == stage) {
            Some(s) => {
                s.done = done;
                s.total = s.total.max(total);
                s.last_ms = elapsed_ms;
                if total > 0 && done >= total {
                    s.finished = true;
                }
            }
            None => self.stages.push(StageState {
                name: stage.to_string(),
                done,
                total,
                started_ms: elapsed_ms,
                last_ms: elapsed_ms,
                finished: total > 0 && done >= total,
            }),
        }
    }

    fn on_finished(&mut self, s: RunSummary) {
        // The engine emits one `Counters` event per cluster and a final one that
        // is the total of all of them; accumulating the per-cluster ones and
        // then taking the summary's value is the only reading in which the
        // counters neither jump nor double-count.
        self.counters = s.counters.clone();
        self.summary = Some(s);
        for st in &mut self.stages {
            st.finished = true;
            if st.total > 0 {
                st.done = st.total;
            }
        }
        self.finished = true;
    }

    /// A note the interface itself wants to show in the same list.
    pub fn push_local_note(&mut self, text: String) {
        self.notes.push(Note {
            t_ms: 0.0,
            text,
            local: true,
        });
    }

    // -- derived views ------------------------------------------------------

    /// The tree to draw: the backbone with every cluster the user opened grafted
    /// in at its pivot leaf, and every other cluster collapsed to a triangle.
    ///
    /// Every branch length comes out in SNPs: the backbone (sketch distances)
    /// and Focus 0 cluster trees are multiplied by the median genome length and
    /// stay marked approximate, Focus 1 trees are piece counts (≈), Focus 2
    /// trees are the engine's SNP counts.
    pub fn display_tree(&self, expanded: &HashSet<u32>) -> Option<DisplayNode> {
        let bb = self.backbone.as_ref()?;
        let ids = self.tree_name_index();
        let (factor, approx) = self.backbone_to_snps();
        let mut top = layout::from_core(bb, fmt::FOCUS_SKETCH, approx, &ids);
        top.scale_lengths(factor);
        self.graft_clusters(&mut top, expanded, &ids);
        Some(if matches!(top.kind, NodeKind::Internal) {
            top.kind = NodeKind::Root;
            top.len = 0.0;
            top
        } else {
            // A run with one cluster has a single-leaf backbone, so the pivot
            // leaf *is* the root. Wrap it: overwriting its kind would throw away
            // the cluster's own node, its triangle and its length region.
            let mut root = DisplayNode::new(NodeKind::Root, 0.0, fmt::FOCUS_SKETCH, true);
            root.children = vec![top];
            root
        })
    }

    fn graft_clusters(
        &self,
        node: &mut DisplayNode,
        expanded: &HashSet<u32>,
        ids: &HashMap<String, u32>,
    ) {
        if let NodeKind::Genome { name, .. } = &node.kind {
            if let Some(&cid) = self.pivot_to_cluster.get(name.as_str()) {
                if let Some(c) = self.cluster(cid) {
                    let len = node.len;
                    let members = c.size();
                    let focus = c.focus;
                    let approximate = c.approximate;
                    let opened = expanded.contains(&cid);
                    let pivot_name = c.pivot_name.clone();
                    let single = members <= 1;
                    let mut sub = match c.tree.as_ref().filter(|_| opened) {
                        Some(t) => Some(self.cluster_subtree(t, focus, approximate, ids)),
                        // A cluster of one genome is trivially exact and the
                        // engine emits no tree for it, so draw the pivot itself.
                        None if opened && single => {
                            let mut s =
                                DisplayNode::new(NodeKind::Internal, 0.0, focus, approximate);
                            s.children = vec![DisplayNode::new(
                                NodeKind::Genome {
                                    id: c.pivot,
                                    name: pivot_name.clone(),
                                },
                                0.0,
                                focus,
                                false,
                            )];
                            Some(s)
                        }
                        None => None,
                    };
                    if let Some(s) = &mut sub {
                        s.kind = NodeKind::ClusterRoot {
                            id: cid,
                            members,
                            focus,
                        };
                        // The backbone edge stays on the cluster root; the root
                        // stub of a midpoint-rooted tree has no length of its
                        // own to lose.
                        s.len = len;
                        *node = s.clone();
                    } else {
                        node.kind = NodeKind::ClusterCollapsed {
                            id: cid,
                            members,
                            focus,
                        };
                    }
                }
            }
            return;
        }
        for c in &mut node.children {
            self.graft_clusters(c, expanded, ids);
        }
    }

    /// A cluster's own tree as display nodes, in SNPs: a Focus 0 star is a
    /// sketch distance and is scaled by the median genome length; Focus 1 and 2
    /// trees are already counts.
    fn cluster_subtree(
        &self,
        t: &Tree,
        focus: u8,
        approximate: bool,
        ids: &HashMap<String, u32>,
    ) -> DisplayNode {
        let mut d = layout::from_core(t, focus, approximate, ids);
        if focus == fmt::FOCUS_SKETCH {
            d.scale_lengths(self.median_length() as f64);
        }
        d
    }

    /// One cluster on its own, as a standard phylogram in SNPs: the "this
    /// cluster only" view of the tree panel. A cluster without a tree yet is a
    /// collapsed triangle; a one-genome cluster is its pivot.
    pub fn cluster_display_tree(&self, cid: u32) -> Option<DisplayNode> {
        let c = self.cluster(cid)?;
        let ids = self.tree_name_index();
        let kind = NodeKind::ClusterRoot {
            id: cid,
            members: c.size(),
            focus: c.focus,
        };
        let sub = match &c.tree {
            Some(t) => {
                let mut s = self.cluster_subtree(t, c.focus, c.approximate, &ids);
                s.kind = kind;
                s.len = 0.0;
                s
            }
            None if c.size() <= 1 => {
                let mut s = DisplayNode::new(kind, 0.0, c.focus, c.approximate);
                s.children = vec![DisplayNode::new(
                    NodeKind::Genome {
                        id: c.pivot,
                        name: c.pivot_name.clone(),
                    },
                    0.0,
                    c.focus,
                    false,
                )];
                s
            }
            None => DisplayNode::new(
                NodeKind::ClusterCollapsed {
                    id: cid,
                    members: c.size(),
                    focus: c.focus,
                },
                0.0,
                c.focus,
                c.approximate,
            ),
        };
        let mut root = DisplayNode::new(NodeKind::Root, 0.0, c.focus, c.approximate);
        root.children = vec![sub];
        Some(root)
    }

    /// The cluster a genome belongs to, or `None` before clustering.
    pub fn cluster_of(&self, genome: u32) -> Option<u32> {
        self.genome(genome)?.cluster
    }

    /// Members of a cluster that have map data, in id order.
    pub fn mapped_members(&self, cluster: u32) -> Vec<u32> {
        let mut v: Vec<u32> = self
            .maps
            .get(&cluster)
            .map(|m| m.keys().copied().collect())
            .unwrap_or_default();
        v.sort_unstable();
        v
    }
}

// ---------------------------------------------------------------------------
// Sorting and filtering the genome list
// ---------------------------------------------------------------------------

/// A column of the genome list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortColumn {
    Id,
    Name,
    Cluster,
    Length,
    Contigs,
    N50,
    Gc,
    Focus,
}

/// Which column the genome list is ordered by, and in which direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SortSpec {
    pub column: SortColumn,
    pub ascending: bool,
}

impl Default for SortSpec {
    fn default() -> Self {
        Self {
            column: SortColumn::Id,
            ascending: true,
        }
    }
}

impl SortSpec {
    /// Flip the direction, as clicking a sorted header does.
    pub fn toggled(self, column: SortColumn) -> Self {
        if self.column == column {
            Self {
                column,
                ascending: !self.ascending,
            }
        } else {
            // A fresh column starts ascending, except for the sizes, where the
            // interesting end is the large one.
            Self {
                column,
                ascending: !matches!(
                    column,
                    SortColumn::Length | SortColumn::Contigs | SortColumn::N50
                ),
            }
        }
    }
}

/// True when a row passes the filter. The empty filter passes everything.
pub fn matches_filter(g: &GenomeRow, filter: &str) -> bool {
    let f = filter.trim();
    if f.is_empty() {
        return true;
    }
    g.name.to_lowercase().contains(&f.to_lowercase())
}

/// Indices into `genomes` of the rows that pass `filter`, ordered by `spec`.
///
/// Ties are always broken by row index, so the list never reorders itself
/// between two frames while the engine is still adding rows.
pub fn sorted_indices(genomes: &[GenomeRow], filter: &str, spec: &SortSpec) -> Vec<usize> {
    let mut idx: Vec<usize> = genomes
        .iter()
        .enumerate()
        .filter(|(_, g)| matches_filter(g, filter))
        .map(|(i, _)| i)
        .collect();
    idx.sort_by(|&a, &b| {
        let (ga, gb) = (&genomes[a], &genomes[b]);
        let ord = match spec.column {
            SortColumn::Id => ga.id.cmp(&gb.id),
            SortColumn::Name => ga.name.to_lowercase().cmp(&gb.name.to_lowercase()),
            // Ungrouped genomes sort last in both directions.
            SortColumn::Cluster => ga
                .cluster
                .map(|c| (0u8, c))
                .unwrap_or((1, 0))
                .cmp(&gb.cluster.map(|c| (0u8, c)).unwrap_or((1, 0))),
            SortColumn::Length => ga.length.cmp(&gb.length),
            SortColumn::Contigs => ga.contigs.cmp(&gb.contigs),
            SortColumn::N50 => ga.n50.cmp(&gb.n50),
            SortColumn::Gc => ga
                .gc
                .partial_cmp(&gb.gc)
                .unwrap_or(std::cmp::Ordering::Equal),
            SortColumn::Focus => ga.focus.cmp(&gb.focus),
        };
        let ord = if spec.ascending { ord } else { ord.reverse() };
        ord.then(a.cmp(&b))
    });
    idx
}

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

/// A `GenomeScanned` event, so that any module's tests can drive the fold.
#[cfg(test)]
pub fn scan_event(id: u32, name: &str, length: u64, contigs: u32) -> Event {
    Event::GenomeScanned(infinity_core::events::GenomeScannedInfo {
        id,
        name: name.to_string(),
        length,
        contigs,
        n50: length,
        gc: 0.5,
        sketch_size: 100,
    })
}

/// A `Clustered` event. Each entry is `(cluster id, pivot id, member ids)`; the
/// member names are `g<id>`, matching [`scan_event`].
#[cfg(test)]
pub fn clustered_event(specs: &[(u32, u32, &[u32])]) -> Event {
    let clusters = specs
        .iter()
        .map(|&(id, pivot, members)| infinity_core::events::ClusterInfo {
            id,
            pivot,
            pivot_name: format!("g{pivot}"),
            members: members.to_vec(),
            member_names: members.iter().map(|m| format!("g{m}")).collect(),
            focus: 0,
            cluster_dist: 0.005,
        })
        .collect();
    Event::Clustered {
        threshold: 0.005,
        clusters,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::{stage, PieceView};
    use infinity_core::params::Parameters;
    use infinity_core::variant::{Variant, VariantKind};
    fn scanned(id: u32, name: &str, length: u64, contigs: u32) -> Event {
        Event::GenomeScanned(GenomeScannedInfo {
            id,
            name: name.to_string(),
            length,
            contigs,
            n50: length,
            gc: 0.5,
            sketch_size: 100,
        })
    }

    fn cluster(id: u32, pivot: u32, members: &[u32], names: &[&str]) -> ClusterInfo {
        ClusterInfo {
            id,
            pivot,
            pivot_name: names[0].to_string(),
            members: members.to_vec(),
            member_names: names.iter().map(|s| s.to_string()).collect(),
            focus: 0,
            cluster_dist: 0.005,
        }
    }

    /// Three genomes in two clusters, folded the way a real run arrives.
    fn three_genomes() -> UiModel {
        let mut m = UiModel::new(Source::Live);
        m.apply(0.0, scanned(0, "aaa", 1000, 1));
        m.apply(0.0, scanned(1, "bbb", 2000, 2));
        m.apply(0.0, scanned(2, "ccc", 3000, 3));
        m.apply(
            1.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![
                    cluster(0, 0, &[0, 1], &["aaa", "bbb"]),
                    cluster(1, 2, &[2], &["ccc"]),
                ],
            },
        );
        m.apply(
            2.0,
            Event::BackboneTree {
                newick: "(aaa:0.001,ccc:0.002);".to_string(),
                pivots: vec!["aaa".into(), "ccc".into()],
            },
        );
        m
    }

    #[test]
    fn scanned_genomes_appear_as_rows() {
        let m = three_genomes();
        assert_eq!(m.genome_count(), 3);
        assert_eq!(m.total_length(), 6000);
        assert_eq!(m.genome(1).unwrap().contigs, 2);
    }

    #[test]
    fn a_rescan_replaces_the_row_instead_of_duplicating_it() {
        let mut m = UiModel::new(Source::Live);
        m.apply(0.0, scanned(0, "a", 100, 1));
        m.apply(0.0, scanned(0, "a", 400, 4));
        assert_eq!(m.genome_count(), 1);
        assert_eq!(m.genome(0).unwrap().length, 400);
    }

    #[test]
    fn clustering_assigns_members_pivots_and_colours() {
        let m = three_genomes();
        assert_eq!(m.cluster_count(), 2);
        assert_eq!(m.threshold, 0.005);
        assert_eq!(m.cluster_of(0), Some(0));
        assert_eq!(m.cluster_of(2), Some(1));
        assert!(m.genome(0).unwrap().is_pivot);
        assert!(!m.genome(1).unwrap().is_pivot);
        assert_eq!(m.cluster(0).unwrap().pivot_name, "aaa");
    }

    #[test]
    fn a_singleton_cluster_is_exact_because_there_is_nothing_to_compare() {
        let m = three_genomes();
        assert_eq!(m.cluster(1).unwrap().focus, fmt::FOCUS_EXACT);
        assert!(!m.cluster(1).unwrap().approximate);
        assert_eq!(m.exact_clusters(), 1);
        assert_eq!(m.genome(2).unwrap().focus, fmt::FOCUS_EXACT);
    }

    #[test]
    fn focus_events_raise_a_cluster_and_its_genomes() {
        let mut m = three_genomes();
        m.apply(
            3.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 1,
                approximate: true,
                newick: "(aaa:1,bbb:2);".to_string(),
                members: vec![
                    MemberNumbers {
                        id: 0,
                        name: "aaa".into(),
                        diff_pieces: 3,
                        snps: 0,
                        indels: 0,
                        uncovered_bp: 0,
                        novel_bp: 0,
                    },
                    MemberNumbers {
                        id: 1,
                        name: "bbb".into(),
                        diff_pieces: 5,
                        snps: 0,
                        indels: 0,
                        uncovered_bp: 10,
                        novel_bp: 2,
                    },
                ],
            }),
        );
        assert_eq!(m.cluster(0).unwrap().focus, 1);
        assert!(m.cluster(0).unwrap().approximate);
        assert!(m.cluster(0).unwrap().tree.is_some());
        assert_eq!(
            m.genome(1).unwrap().numbers.as_ref().unwrap().diff_pieces,
            5
        );
        assert_eq!(m.genome(1).unwrap().focus, 1);
        // Cluster 1 is a singleton, which is exact from the start, so it counts
        // as refined too.
        assert_eq!(m.refined_clusters(), 2);
        assert_eq!(m.exact_clusters(), 1);
    }

    #[test]
    fn a_focus_zero_event_never_downgrades_an_unchanged_cluster() {
        let mut m = three_genomes();
        let sharp = Event::ClusterFocus(ClusterFocusInfo {
            cluster: 0,
            focus: 2,
            approximate: false,
            newick: "(aaa:1,bbb:2);".to_string(),
            members: vec![],
        });
        m.apply(3.0, sharp);
        assert_eq!(m.cluster(0).unwrap().focus, 2);
        // What the engine emits for every cluster after a re-clustering.
        let blurry = Event::ClusterFocus(ClusterFocusInfo {
            cluster: 0,
            focus: 0,
            approximate: true,
            newick: "(aaa:0,bbb:0);".to_string(),
            members: vec![],
        });
        m.apply(4.0, blurry);
        assert_eq!(
            m.cluster(0).unwrap().focus,
            2,
            "an exact tree must not be shown as blurry again"
        );
        assert!(m.cluster(0).unwrap().tree.is_some());
        // A genuinely new membership does replace the tree.
        m.apply(
            5.0,
            Event::Clustered {
                threshold: 0.02,
                clusters: vec![cluster(0, 0, &[0, 1, 2], &["aaa", "bbb", "ccc"])],
            },
        );
        m.apply(
            6.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 0,
                approximate: true,
                newick: "(aaa:0,(bbb:0,ccc:0));".to_string(),
                members: vec![],
            }),
        );
        assert_eq!(m.cluster(0).unwrap().focus, 0);
    }

    #[test]
    fn recluster_keeps_results_for_an_unchanged_member_set() {
        let mut m = three_genomes();
        m.apply(
            3.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(aaa:1,bbb:2);".to_string(),
                members: vec![],
            }),
        );
        // Same two clusters, same members, new threshold.
        m.apply(
            4.0,
            Event::Clustered {
                threshold: 0.01,
                clusters: vec![
                    cluster(0, 0, &[0, 1], &["aaa", "bbb"]),
                    cluster(1, 2, &[2], &["ccc"]),
                ],
            },
        );
        assert_eq!(m.threshold, 0.01);
        assert_eq!(
            m.cluster(0).unwrap().focus,
            2,
            "cached by the engine, so kept"
        );
        assert!(m.cluster(0).unwrap().tree.is_some());
    }

    #[test]
    fn counters_accumulate_and_the_summary_overrides_them() {
        let mut m = UiModel::new(Source::Live);
        let c1 = Counters {
            pieces_compared: 10,
            ..Default::default()
        };
        let c2 = Counters {
            pieces_compared: 5,
            ..Default::default()
        };
        m.apply(0.0, Event::Counters(c1));
        m.apply(0.0, Event::Counters(c2));
        // The engine emits one Counters event per cluster, so the interface sums
        // them; the summary then carries the authoritative total.
        assert_eq!(m.counters.pieces_compared, 15);
        let total = Counters {
            pieces_compared: 15,
            ..Default::default()
        };
        let s = RunSummary {
            counters: total,
            ..RunSummary::new("test", Parameters::default())
        };
        m.apply(1.0, Event::Finished(Box::new(s)));
        assert_eq!(m.counters.pieces_compared, 15, "not 30");
        assert!(m.finished);
    }

    #[test]
    fn stages_are_tracked_in_arrival_order_with_their_own_elapsed_time() {
        let mut m = UiModel::new(Source::Live);
        m.apply(0.0, progress(stage::SCAN, 0, 10, 100));
        m.apply(0.0, progress(stage::SCAN, 5, 10, 600));
        m.apply(0.0, progress(stage::SCAN, 10, 10, 1100));
        m.apply(0.0, progress(stage::CLUSTER, 1, 1, 1200));
        assert_eq!(m.stages.len(), 2);
        assert_eq!(m.stages[0].name, "scan");
        assert!(m.stages[0].finished);
        assert_eq!(m.stages[0].elapsed_ms(), 1000);
        assert_eq!(m.stages[0].fraction(), Some(1.0));
        assert_eq!(m.stages[1].fraction(), Some(1.0));
        // A stage that never reports a total has no fraction, not zero.
        m.apply(0.0, progress(stage::DISCOVER, 0, 0, 10));
        assert_eq!(m.stages[2].fraction(), None);
    }

    fn progress(s: &str, done: u64, total: u64, elapsed_ms: u64) -> Event {
        Event::Progress {
            stage: s.to_string(),
            done,
            total,
            elapsed_ms,
        }
    }

    #[test]
    fn refine_progress_never_reports_fewer_clusters_than_are_on_screen() {
        let mut m = three_genomes();
        m.apply(
            3.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(aaa:1,bbb:2);".to_string(),
                members: vec![],
            }),
        );
        // After a re-clustering the engine's counter restarts at zero.
        m.apply(0.0, progress(stage::REFINE, 0, 2, 500));
        let s = m.stages.iter().find(|s| s.name == stage::REFINE).unwrap();
        assert_eq!(
            s.done,
            m.refined_clusters() as u64,
            "never fewer than the clusters already refined on screen"
        );
        assert_eq!(s.done, 2, "cluster 0 was refined, the singleton counts too");
    }

    #[test]
    fn a_bad_newick_string_is_reported_and_does_not_break_the_fold() {
        let mut m = UiModel::new(Source::Live);
        m.apply(0.0, scanned(0, "a", 10, 1));
        m.apply(
            0.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![cluster(0, 0, &[0], &["a"])],
            },
        );
        m.apply(
            1.0,
            Event::BackboneTree {
                newick: "((((".to_string(),
                pivots: vec![],
            },
        );
        assert!(m.backbone.is_none());
        assert!(m.notes.iter().any(|n| n.text.contains("backbone")));
    }

    #[test]
    fn genome_ids_are_resolved_by_name_when_the_engine_renumbered_them() {
        let mut m = UiModel::new(Source::Live);
        // File 1 failed to read, so the stream ids skip it while the engine's
        // own ids do not.
        m.apply(0.0, scanned(0, "aaa", 10, 1));
        m.apply(0.0, scanned(2, "ccc", 30, 1));
        m.apply(
            1.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![cluster(0, 0, &[0, 1], &["aaa", "ccc"])],
            },
        );
        let c = m.cluster(0).unwrap();
        assert_eq!(c.members, vec![0, 2], "id 1 was mapped back to ccc (id 2)");
        assert_eq!(m.genome(2).unwrap().cluster, Some(0));
    }

    #[test]
    fn duplicate_names_are_not_guessed_at() {
        let mut m = UiModel::new(Source::Live);
        m.apply(0.0, scanned(0, "same", 10, 1));
        m.apply(0.0, scanned(1, "same", 20, 1));
        assert!(m.dup_names.contains("same"));
        assert!(m.tree_name_index().is_empty());
        // The engine's own de-duplicated labels restore the link.
        m.apply(
            1.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![cluster(0, 0, &[0, 1], &["same", "same#1"])],
            },
        );
        let idx = m.tree_name_index();
        assert_eq!(idx.get("same#1"), Some(&1));
        assert_eq!(idx.get("same"), Some(&0));
    }

    #[test]
    fn map_updates_are_stored_per_cluster_and_genome() {
        let mut m = three_genomes();
        m.apply(
            4.0,
            Event::MapUpdate(Box::new(MapUpdateInfo {
                cluster: 0,
                genome: 1,
                genome_name: "bbb".into(),
                pivot_name: "aaa".into(),
                pieces: vec![PieceView {
                    pivot_contig: 0,
                    pivot_contig_name: "c0".into(),
                    start: 10,
                    end: 40,
                    resolved: true,
                }],
                variants: vec![snp(0, 20)],
                uncovered: vec![(900, 950)],
                pivot_contigs: vec!["c0".into()],
                pivot_contig_lens: vec![1000],
            })),
        );
        let mu = m.map(0, 1).expect("map data");
        assert_eq!(mu.variants.len(), 1);
        assert_eq!(m.mapped_members(0), vec![1]);
    }

    fn snp(contig: u32, pos: u64) -> Variant {
        Variant {
            kind: VariantKind::Snp,
            pivot_contig: contig,
            pivot_start: pos,
            pivot_end: pos + 1,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            genome_contig: 0,
            genome_start: pos,
            genome_end: pos + 1,
            strand: '+',
            flags: 0,
        }
    }

    #[test]
    fn the_display_tree_collapses_clusters_by_default_and_opens_them_on_request() {
        let m = three_genomes();
        let collapsed = m.display_tree(&HashSet::new()).expect("a backbone exists");
        assert_eq!(collapsed.children.len(), 2);
        assert!(matches!(
            collapsed.children[0].kind,
            NodeKind::ClusterCollapsed {
                id: 0,
                members: 2,
                ..
            }
        ));
        // The singleton cluster is opened, so its genome becomes a leaf.
        let mut open = HashSet::new();
        open.insert(1u32);
        let t = m.display_tree(&open).unwrap();
        let c = t
            .children
            .iter()
            .find(|c| matches!(c.kind, NodeKind::ClusterRoot { id: 1, .. }))
            .unwrap();
        assert_eq!(c.children.len(), 1);
        assert_eq!(c.children[0].kind.name(), Some("ccc"));
        // The cluster's own focus, not the backbone's, travels with the subtree.
        assert_eq!(c.children[0].focus, fmt::FOCUS_EXACT);
        assert_eq!(collapsed.children[0].focus, fmt::FOCUS_SKETCH);
    }

    #[test]
    fn every_tree_is_drawn_in_snps() {
        let mut m = three_genomes();
        // Median genome length is 2000 (lengths 1000, 2000, 3000).
        assert_eq!(m.median_length(), 2000);
        let (factor, approx) = m.backbone_to_snps();
        assert_eq!(factor, 2000.0);
        assert!(approx);
        // The backbone edge of 0.001 subs/site becomes ≈ 2 SNPs.
        let t = m.display_tree(&HashSet::new()).unwrap();
        assert!(
            (t.children[0].len - 2.0).abs() < 1e-9,
            "{}",
            t.children[0].len
        );
        assert!(t.children[0].approximate);
        // A Focus 2 cluster tree keeps its SNP counts as they are.
        m.apply(
            3.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(aaa:1,bbb:4);".to_string(),
                members: vec![],
            }),
        );
        let mut open = HashSet::new();
        open.insert(0u32);
        let t = m.display_tree(&open).unwrap();
        let c = &t.children[0];
        assert!(matches!(c.kind, NodeKind::ClusterRoot { id: 0, .. }));
        assert!((c.len - 2.0).abs() < 1e-9, "the backbone edge stays ≈ SNPs");
        assert_eq!(c.children[1].len, 4.0);
        assert!(!c.children[1].approximate);
        // A Focus 0 star is a sketch distance, so it is scaled like the backbone.
        m.apply(
            4.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![cluster(0, 0, &[0, 1, 2], &["aaa", "bbb", "ccc"])],
            },
        );
        m.apply(
            5.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 0,
                approximate: true,
                newick: "(aaa:0,bbb:0.002,ccc:0.001);".to_string(),
                members: vec![],
            }),
        );
        let t = m.cluster_display_tree(0).unwrap();
        let c = &t.children[0];
        assert!((c.children[1].len - 4.0).abs() < 1e-9, "0.002 × 2000");
        assert!(c.children[1].approximate);
    }

    #[test]
    fn a_cluster_on_its_own_is_a_phylogram_that_fills_the_width() {
        let mut m = three_genomes();
        m.apply(
            3.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(aaa:1,bbb:3);".to_string(),
                members: vec![],
            }),
        );
        let t = m.cluster_display_tree(0).unwrap();
        let l = layout::layout(&t);
        assert_eq!(l.rows, 2.0);
        assert_eq!(l.regions.len(), 2);
        // No empty backbone width is reserved in front of the cluster.
        assert!((l.x_max - layout::CLUSTER_W).abs() < 1e-9, "{}", l.x_max);
        // A cluster with no tree yet is its triangle; a singleton its pivot.
        let mut m2 = three_genomes();
        m2.clusters[0].tree = None;
        let t = m2.cluster_display_tree(0).unwrap();
        assert!(matches!(
            t.children[0].kind,
            NodeKind::ClusterCollapsed { id: 0, .. }
        ));
        let t = m2.cluster_display_tree(1).unwrap();
        assert_eq!(t.children[0].children[0].kind.name(), Some("ccc"));
        assert!(m2.cluster_display_tree(99).is_none());
    }

    #[test]
    fn the_median_of_nothing_is_zero() {
        assert_eq!(UiModel::default().median_length(), 0);
    }

    #[test]
    fn without_a_backbone_there_is_no_tree_to_draw() {
        let m = UiModel::new(Source::Empty);
        assert!(m.display_tree(&HashSet::new()).is_none());
    }

    #[test]
    fn a_single_cluster_backbone_keeps_its_cluster_node() {
        // One cluster means a one-leaf backbone, so the pivot leaf is the root.
        // It must still come out as a cluster node, or the tree loses its
        // triangle and its own length region.
        let mut m = UiModel::new(Source::Live);
        for i in 0..3u32 {
            m.apply(0.0, scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(0.0, clustered_event(&[(0, 0, &[0, 1, 2])]));
        m.apply(
            0.0,
            Event::BackboneTree {
                newick: "g0;".to_string(),
                pivots: vec!["g0".into()],
            },
        );
        m.apply(
            0.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,(g1:2,g2:3));".to_string(),
                members: vec![],
            }),
        );

        let collapsed = m.display_tree(&HashSet::new()).unwrap();
        assert_eq!(collapsed.kind, NodeKind::Root);
        assert_eq!(collapsed.children.len(), 1);
        assert!(matches!(
            collapsed.children[0].kind,
            NodeKind::ClusterCollapsed {
                id: 0,
                members: 3,
                focus: 2
            }
        ));
        let l = layout::layout(&collapsed);
        assert_eq!(l.rows, 1.0);

        let mut open = HashSet::new();
        open.insert(0u32);
        let opened = m.display_tree(&open).unwrap();
        assert!(matches!(
            opened.children[0].kind,
            NodeKind::ClusterRoot {
                id: 0,
                members: 3,
                focus: 2
            }
        ));
        let l = layout::layout(&opened);
        assert_eq!(l.rows, 3.0, "all three genomes become leaves");
        assert_eq!(l.regions.len(), 2, "the backbone and the opened cluster");
        assert_eq!(l.regions[1].unit, "SNPs");
        assert_eq!(l.regions[1].cluster, Some(0));
    }

    #[test]
    fn filtering_is_case_insensitive_and_substring_based() {
        let m = three_genomes();
        assert_eq!(
            sorted_indices(&m.genomes, "", &SortSpec::default()).len(),
            3
        );
        assert_eq!(
            sorted_indices(&m.genomes, "BB", &SortSpec::default()).len(),
            1
        );
        assert_eq!(
            sorted_indices(&m.genomes, "zz", &SortSpec::default()).len(),
            0
        );
    }

    #[test]
    fn sorting_orders_by_the_chosen_column_and_breaks_ties_by_id() {
        let m = three_genomes();
        let by_len = SortSpec {
            column: SortColumn::Length,
            ascending: false,
        };
        let ids: Vec<u32> = sorted_indices(&m.genomes, "", &by_len)
            .iter()
            .map(|&i| m.genomes[i].id)
            .collect();
        assert_eq!(ids, vec![2, 1, 0]);
        let by_len_asc = SortSpec {
            column: SortColumn::Length,
            ascending: true,
        };
        let ids: Vec<u32> = sorted_indices(&m.genomes, "", &by_len_asc)
            .iter()
            .map(|&i| m.genomes[i].id)
            .collect();
        assert_eq!(ids, vec![0, 1, 2]);
        // Everything has the same contigs count? No: 1, 2, 3, so name order and
        // id order agree here; check the tie-break explicitly instead.
        let mut rows = m.genomes.clone();
        rows[0].contigs = 5;
        rows[1].contigs = 5;
        let spec = SortSpec {
            column: SortColumn::Contigs,
            ascending: true,
        };
        let ids: Vec<u32> = sorted_indices(&rows, "", &spec)
            .iter()
            .map(|&i| rows[i].id)
            .collect();
        assert_eq!(ids, vec![2, 0, 1]);
    }

    #[test]
    fn clicking_a_sorted_header_flips_it_and_a_new_one_starts_sensibly() {
        let s = SortSpec::default();
        assert!(s.ascending && s.column == SortColumn::Id);
        // Clicking the column that is already sorted flips it.
        assert!(!s.toggled(SortColumn::Id).ascending);
        assert!(s.toggled(SortColumn::Id).toggled(SortColumn::Id).ascending);
        // A size column starts descending: the big genomes are the interesting
        // ones. Everything else starts ascending.
        assert!(!s.toggled(SortColumn::Length).ascending);
        assert!(!s.toggled(SortColumn::Contigs).ascending);
        assert!(!s.toggled(SortColumn::N50).ascending);
        assert!(s.toggled(SortColumn::Name).ascending);
        assert!(s.toggled(SortColumn::Focus).ascending);
        // Clicking the same new column again flips it.
        assert!(
            s.toggled(SortColumn::Length)
                .toggled(SortColumn::Length)
                .ascending
        );
    }

    #[test]
    fn search_finds_the_first_matching_genome() {
        let m = three_genomes();
        assert_eq!(m.search("bb"), Some(1));
        assert_eq!(m.search("B"), Some(1));
        assert_eq!(m.search("  cc "), Some(2));
        assert_eq!(m.search(""), None);
        assert_eq!(m.search("zzz"), None);
    }
}
