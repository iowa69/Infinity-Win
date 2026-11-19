//! The window itself: one `eframe::App` that is a thin renderer over
//! [`UiModel`].
//!
//! Three rules shape this file (section 8, and the owner's priorities):
//!
//! 1. **Nothing is computed here.** The engine runs on its own threads, an
//!    opened folder is read on its own thread, and the only work the update loop
//!    does is drain a bounded batch of events and draw.
//! 2. **The window repaints because something asked it to.** A relay thread
//!    wakes the interface when an event arrives; otherwise a heartbeat runs only
//!    while a run is going, and the window is idle (zero CPU) when nothing is.
//! 3. **Focus is never ambiguous.** Every number that is an estimate is drawn
//!    with the "≈" marker and in the colour of its focus level.

use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Duration;

use eframe::egui;
use egui::Context;

use infinity_core::engine::RunResults;

use crate::layout::{self, Positions, TreeLayout, ANIM_SECS};
use crate::open::OpenHandle;
use crate::pairwise::{self, PairView};
use crate::prefs::{Prefs, MAX_RECENT};
use crate::run::{RunHandle, BATCH};
use crate::selection::{ClickMods, Selection};
use crate::state::{SortSpec, Source, UiModel};
use crate::{export, open, prefs, run};

/// How often the interface wakes itself while a run is going, even when no event
/// has arrived: the status bar shows elapsed seconds, and the relay thread's
/// wake-up is the fast path this only backs up.
const HEARTBEAT: Duration = Duration::from_millis(100);

/// The desktop program.
pub struct InfinityApp {
    /// Everything the panels draw, folded out of the event stream.
    pub(crate) model: UiModel,
    pub(crate) prefs: Prefs,

    // --- inputs -----------------------------------------------------------
    pub(crate) inputs: Vec<PathBuf>,
    pub(crate) input_text: String,
    pub(crate) outdir_text: String,
    /// What the command line asked for, done once on the first frame.
    pub(crate) startup: Startup,
    pub(crate) open_text: String,
    pub(crate) export_text: String,

    // --- run state --------------------------------------------------------
    pub(crate) run_handle: Option<RunHandle>,
    pub(crate) results: Option<RunResults>,
    pub(crate) opening: Option<OpenHandle>,
    /// Result of an asynchronous PNG capture.
    png_result: Arc<Mutex<Option<Result<PathBuf, String>>>>,

    // --- selection --------------------------------------------------------
    /// Any number of genomes; see [`crate::selection`] for the click rules.
    pub(crate) selection: Selection,
    pub(crate) selected_cluster: Option<u32>,
    pub(crate) expanded: HashSet<u32>,

    // --- genome list ------------------------------------------------------
    pub(crate) filter: String,
    pub(crate) sort: SortSpec,

    // --- tree view --------------------------------------------------------
    pub(crate) tree: TreeView,
    pub(crate) tree_layout: Option<TreeLayout>,
    pub(crate) positions: Vec<(f64, f64)>,
    anim_from: Positions,
    anim_start: Option<f64>,
    laid_out_version: u64,
    laid_out_expanded: u64,
    /// Screen rectangle of the tree canvas and the pixels-per-point in force,
    /// remembered so a screenshot can be cropped to the tree alone.
    pub(crate) tree_rect: Option<(egui::Rect, f32)>,

    // --- map view ---------------------------------------------------------
    pub(crate) map_view: crate::map_panel::MapView,

    // --- misc -------------------------------------------------------------
    pub(crate) status: String,
    pub(crate) params_error: Option<String>,
    pub(crate) show_notes: bool,
    /// Whether the "add genomes" input area under the toolbar is open.
    pub(crate) show_inputs: bool,
    /// Whether the right-hand panel is open.
    pub(crate) show_side: bool,
    /// Whether the Settings window is open.
    pub(crate) show_settings: bool,
    /// The SNP threshold being typed in Settings, and why it was refused.
    pub(crate) threshold_text: String,
    pub(crate) threshold_error: Option<String>,
    /// Set once the theme has been installed into the egui context.
    theme_installed: bool,
    /// Whether the last frame was drawn in the dark theme (for exports that
    /// draw "like the screen" without a `Ui` at hand).
    pub(crate) dark: bool,
    /// The cluster distance the slider last sent, so dragging it does not queue
    /// a re-clustering per pixel.
    last_dist_sent: f64,

    // --- divergent assemblies ----------------------------------------------
    /// The warning dialog, while it is open.
    pub(crate) outlier_dialog: Option<crate::outlier_dialog::OutlierDialog>,
    /// The dialog was offered for this run already (it opens once per run).
    pub(crate) outliers_prompted: bool,
    /// "Continue anyway" was chosen: the Run panel keeps a warning.
    pub(crate) outliers_acknowledged: bool,
    /// A run was cancelled to start again without some genomes; it starts as
    /// soon as the old one has stopped.
    pub(crate) restart_pending: bool,
}

/// Pan, zoom and decoration state of the tree panel.
#[derive(Debug, Clone)]
pub struct TreeView {
    /// Pixels per leaf row.
    pub px_per_row: f32,
    /// Pixels per horizontal layout unit.
    pub px_per_unit: f32,
    /// Pan offset in pixels.
    pub offset: egui::Vec2,
    pub show_legend: bool,
    pub search: String,
    /// Set when the view must be refitted to the tree, consumed by the panel.
    pub needs_fit: bool,
    /// True once the user has panned or zoomed by hand: from then on a growing
    /// tree no longer yanks the view back to fit, until "Fit" is pressed.
    pub user_moved: bool,
    /// Genome the search box found, drawn with a marker.
    pub search_hit: Option<u32>,
    /// A zoom step asked for by a header button, applied to the canvas centre.
    pub zoom_request: Option<f32>,
    /// A layout point the next frame must put in the middle of the panel.
    pub centre: Option<(f64, f64)>,
    /// The whole picture, or the selected cluster on its own.
    pub mode: TreeMode,
    /// SNP counts written on the branches.
    pub branch_labels: bool,
    /// Tip labels lined up in one column, with dotted leaders to the tips.
    pub aligned_tips: bool,
    /// A strip of cluster colour next to every tip label.
    pub color_strip: bool,
}

/// What the tree panel shows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TreeMode {
    /// The backbone with clusters collapsed or opened.
    All,
    /// The selected cluster alone, as a standard phylogram in SNPs.
    Cluster,
}

impl Default for TreeView {
    fn default() -> Self {
        Self {
            px_per_row: 16.0,
            px_per_unit: 420.0,
            offset: egui::Vec2::ZERO,
            show_legend: true,
            search: String::new(),
            needs_fit: true,
            user_moved: false,
            search_hit: None,
            zoom_request: None,
            centre: None,
            mode: TreeMode::All,
            branch_labels: false,
            aligned_tips: false,
            color_strip: true,
        }
    }
}

impl Default for InfinityApp {
    fn default() -> Self {
        Self::new()
    }
}

impl InfinityApp {
    pub fn new() -> Self {
        let loaded = prefs::load();
        let params = loaded.params.clone();
        let outdir_text = if loaded.recent_outputs.is_empty() {
            default_outdir()
        } else {
            loaded.recent_outputs[0].clone()
        };
        let export_text = if loaded.export_dir.is_empty() {
            outdir_text.clone()
        } else {
            loaded.export_dir.clone()
        };
        let mut app = Self {
            model: UiModel::new(Source::Empty),
            outdir_text,
            startup: Startup::default(),
            export_text,
            open_text: String::new(),
            prefs: loaded,
            inputs: Vec::new(),
            input_text: String::new(),
            run_handle: None,
            results: None,
            opening: None,
            png_result: Arc::new(Mutex::new(None)),
            selection: Selection::default(),
            selected_cluster: None,
            expanded: HashSet::new(),
            filter: String::new(),
            sort: SortSpec::default(),
            tree: TreeView::default(),
            tree_layout: None,
            positions: Vec::new(),
            anim_from: Positions::new(),
            anim_start: None,
            laid_out_version: u64::MAX,
            laid_out_expanded: u64::MAX,
            tree_rect: None,
            map_view: crate::map_panel::MapView::default(),
            status: "Add genome files or a folder, then press Run.".to_string(),
            params_error: None,
            show_notes: false,
            // Open at start: the first thing a user has to do is get genomes in,
            // and there is no file dialog to discover.
            show_inputs: true,
            show_side: true,
            show_settings: false,
            threshold_text: String::new(),
            threshold_error: None,
            theme_installed: false,
            dark: true,
            last_dist_sent: params.cluster_dist,
            outlier_dialog: None,
            outliers_prompted: false,
            outliers_acknowledged: false,
            restart_pending: false,
        };
        app.validate_params();
        app
    }

    pub(crate) fn running(&self) -> bool {
        self.run_handle.as_ref().is_some_and(|r| r.is_running())
    }

    /// Remember a path so the next session offers it.
    pub(crate) fn remember_input(&mut self, path: &str) {
        prefs::push_recent(&mut self.prefs.recent_inputs, path, MAX_RECENT);
    }

    pub(crate) fn remember_output(&mut self, path: &str) {
        prefs::push_recent(&mut self.prefs.recent_outputs, path, MAX_RECENT);
    }

    /// Persist the preferences; a failure is only worth a status line.
    pub(crate) fn save_prefs(&mut self) {
        if !self.export_text.trim().is_empty() {
            self.prefs.export_dir = self.export_text.trim().to_string();
        }
        if let Err(e) = prefs::save(&self.prefs) {
            self.status = format!("cannot save the preferences: {e}");
        }
    }

    pub(crate) fn validate_params(&mut self) {
        self.params_error = self
            .prefs
            .params
            .validate()
            .and_then(|()| self.prefs.resolution.validate())
            .err();
    }

    /// The parameter set a run starts with: the edited parameters, with the
    /// Resolution and SNP-rule settings written in.
    pub(crate) fn run_params(&self) -> infinity_core::params::Parameters {
        let mut params = self.prefs.params.clone();
        self.prefs.resolution.apply_to(&mut params);
        params
    }

    // --- selection ---------------------------------------------------------

    /// Select a genome. `additive` adds or removes it (Ctrl-click) instead of
    /// replacing the selection.
    pub(crate) fn select_genome(&mut self, id: u32, additive: bool) {
        let mods = if additive {
            ClickMods::TOGGLE
        } else {
            ClickMods::PLAIN
        };
        self.click_genome(id, mods, &[]);
    }

    /// A click on a genome in a view whose rows are in `order` (needed for a
    /// Shift-click range). The clicked genome's cluster becomes the one the
    /// map shows, when the click left it selected.
    pub(crate) fn click_genome(&mut self, id: u32, mods: ClickMods, order: &[u32]) {
        self.selection.click(id, mods, order);
        if self.selection.contains(id) {
            if let Some(cid) = self.model.cluster_of(id) {
                self.set_cluster(cid, false);
            }
        }
    }

    pub(crate) fn clear_selection(&mut self) {
        self.selection.clear();
    }

    /// The genome ids of the list, in the order and with the filter it shows.
    pub(crate) fn list_order(&self) -> Vec<u32> {
        crate::state::sorted_indices(&self.model.genomes, &self.filter, &self.sort)
            .into_iter()
            .map(|i| self.model.genomes[i].id)
            .collect()
    }

    /// Ctrl+A: every genome the list shows (after the filter).
    pub(crate) fn select_all_visible(&mut self) {
        let order = self.list_order();
        self.selection.select_all(&order);
    }

    /// Selected genomes the tree cannot show because their cluster is drawn
    /// collapsed: `(cluster, how many)`, in cluster order.
    pub(crate) fn hidden_selected(&self) -> Vec<(u32, usize)> {
        if self.selection.is_empty() || self.single_cluster().is_some() {
            return Vec::new();
        }
        let mut out: Vec<(u32, usize)> = Vec::new();
        for c in &self.model.clusters {
            let opened = self.expanded.contains(&c.id) && (c.tree.is_some() || c.size() <= 1);
            if opened {
                continue;
            }
            let n = c
                .members
                .iter()
                .filter(|&&g| self.selection.contains(g))
                .count();
            if n > 0 {
                out.push((c.id, n));
            }
        }
        out
    }

    /// Open every collapsed cluster that holds a selected genome.
    pub(crate) fn expand_selected(&mut self) {
        for (c, _) in self.hidden_selected() {
            self.expanded.insert(c);
        }
    }

    /// Make a cluster the one the tree, the map and the parameters talk about.
    /// `prioritize` also asks the engine to refine it next (6.11.3).
    pub(crate) fn set_cluster(&mut self, id: u32, prioritize: bool) {
        if self.selected_cluster != Some(id) {
            self.selected_cluster = Some(id);
            self.map_view.needs_fit = true;
        }
        if prioritize {
            if let Some(r) = &self.run_handle {
                if r.is_running() {
                    r.prioritize(id);
                    self.status = format!("cluster {id} moved to the front of the queue");
                }
            }
        }
    }

    /// Open or close a cluster in the tree, and ask the engine to refine it.
    pub(crate) fn toggle_cluster(&mut self, id: u32) {
        if self.expanded.contains(&id) {
            self.expanded.remove(&id);
        } else {
            self.expanded.insert(id);
        }
        self.set_cluster(id, true);
    }

    /// Put the tree and map views back to "fit everything", which is what a new
    /// run or a newly opened folder needs.
    pub(crate) fn reset_view(&mut self) {
        self.tree.needs_fit = true;
        self.tree.user_moved = false;
        self.tree.offset = egui::Vec2::ZERO;
        self.tree.search_hit = None;
        self.tree.zoom_request = None;
        self.tree.centre = None;
        self.map_view.needs_fit = true;
        self.map_view.user_moved = false;
        self.map_view.zoom_request = None;
    }

    // --- run control -------------------------------------------------------

    pub(crate) fn start_run(&mut self, ctx: &Context) {
        if self.running() {
            return;
        }
        self.validate_params();
        if self.params_error.is_some() {
            self.status = "fix the parameters before running".to_string();
            return;
        }
        let outdir = PathBuf::from(self.outdir_text.trim());
        if self.outdir_text.trim().is_empty() {
            self.status = "choose an output folder first".to_string();
            return;
        }
        let params = self.run_params();
        match run::start(self.inputs.clone(), params, outdir.clone(), ctx.clone()) {
            Ok(handle) => {
                if let Some(note) = &handle.startup_note {
                    self.model.push_local_note(note.clone());
                    self.show_notes = true;
                }
                self.remember_output(&outdir.display().to_string());
                self.save_prefs();
                self.model = UiModel::new(Source::Live);
                self.results = None;
                self.outlier_dialog = None;
                self.outliers_prompted = false;
                self.outliers_acknowledged = false;
                self.selection.clear();
                self.expanded.clear();
                self.selected_cluster = None;
                self.reset_view();
                self.last_dist_sent = self.prefs.params.cluster_dist;
                self.run_handle = Some(handle);
                self.status = format!("running, writing to {}", outdir.display());
                ctx.request_repaint();
            }
            Err(e) => {
                self.status = e;
                self.show_notes = true;
            }
        }
    }

    pub(crate) fn stop_run(&mut self) {
        if let Some(r) = &self.run_handle {
            r.cancel();
            self.status =
                "cancellation requested; the engine stops after the current genome".to_string();
        }
    }

    /// Leave `names` out from now on and start the same run again: cancel the
    /// one that is going, and start the new one as soon as it has stopped
    /// (both write to the same folder, so they must not overlap).
    pub(crate) fn exclude_and_restart(&mut self, names: &[String], ctx: &Context) {
        let added = crate::outlier_dialog::add_exclusions(&mut self.prefs.params.exclude, names);
        self.save_prefs();
        if self.running() {
            self.stop_run();
            self.restart_pending = true;
            self.status = format!(
                "excluding {added} genome(s): restarting as soon as the current run has stopped"
            );
            ctx.request_repaint();
        } else {
            self.start_run(ctx);
        }
    }

    /// Open the divergent-assembly dialog once per live run, as soon as the
    /// engine reports flagged genomes.
    fn poll_outliers(&mut self) {
        if self.outliers_prompted || self.model.outliers.is_empty() {
            return;
        }
        self.outliers_prompted = true;
        if self.model.source == Source::Live && !self.restart_pending {
            self.outlier_dialog = Some(crate::outlier_dialog::OutlierDialog::new(
                self.model.outliers.clone(),
            ));
        }
    }

    /// Start the run a "restart without these genomes" is waiting for, once
    /// the cancelled one has really stopped.
    fn poll_restart(&mut self, ctx: &Context) {
        if self.restart_pending && !self.running() && self.opening.is_none() {
            self.restart_pending = false;
            self.start_run(ctx);
        }
    }

    /// Send a new cluster distance straight to the engine (section 8: the slider
    /// "takes effect immediately").
    pub(crate) fn push_cluster_dist(&mut self) {
        let v = self.prefs.params.cluster_dist;
        // Only send on a real change: each command makes the engine re-cluster.
        if (v - self.last_dist_sent).abs() < f64::EPSILON {
            return;
        }
        self.last_dist_sent = v;
        match &self.run_handle {
            Some(r) if r.is_running() => {
                r.set_cluster_dist(v);
                self.status = format!("cluster distance {v}: re-clustering from the sketches");
            }
            _ => {
                self.status = format!("cluster distance {v} will be used by the next run");
            }
        }
    }

    pub(crate) fn open_folder(&mut self, dir: &str, ctx: &Context) {
        let dir = dir.trim();
        if dir.is_empty() {
            self.status = "type the path of a finished output folder".to_string();
            return;
        }
        if self.opening.is_some() {
            return;
        }
        match open::start(PathBuf::from(dir), ctx.clone()) {
            Ok(handle) => {
                self.remember_input(dir);
                self.save_prefs();
                self.status = format!("reading {} ...", dir);
                self.opening = Some(handle);
                ctx.request_repaint();
            }
            Err(e) => {
                self.status = e;
            }
        }
    }

    // --- export ------------------------------------------------------------

    pub(crate) fn export_files(&mut self, item: Option<export::Item>) {
        let dir = PathBuf::from(self.export_text.trim());
        if self.export_text.trim().is_empty() {
            self.status = "type an export folder first".to_string();
            return;
        }
        self.prefs.export_dir = self.export_text.clone();
        self.save_prefs();
        let mut notes = Vec::new();
        let run_dir = self.run_dir();
        let ctx = export::ExportCtx {
            run_dir: run_dir.as_deref(),
            layout: self.tree_layout.as_ref(),
            svg: self.svg_options(),
            ..export::ExportCtx::new(&self.model, self.results.as_ref())
        };
        let report = match item {
            Some(i) => {
                let mut r = export::Report::default();
                if self.model.genome_count() == 0 {
                    r.errors.push("nothing to export yet".to_string());
                    r
                } else if let Err(e) = std::fs::create_dir_all(&dir) {
                    r.errors
                        .push(format!("cannot create {}: {e}", dir.display()));
                    r
                } else {
                    match export::write_item(i, &dir, &ctx, &mut notes) {
                        Ok(paths) => r.written.extend(paths),
                        Err(e) => r.errors.push(format!("{}: {e}", i.label())),
                    }
                    r
                }
            }
            None => export::write_all_with(&dir, &ctx),
        };
        for n in &report.notes {
            self.model.push_local_note(n.clone());
        }
        for n in &notes {
            self.model.push_local_note(n.clone());
        }
        self.status = report.summary();
        if !report.ok() {
            self.model.push_local_note(report.errors.join("; "));
            self.show_notes = true;
        }
    }

    /// Carry out an export asked for by the menu or the Run panel. "Everything"
    /// writes every file and then captures the PNG as well.
    pub(crate) fn run_export(&mut self, action: crate::chrome::ExportAction, ctx: &Context) {
        use crate::chrome::ExportAction;
        match action {
            ExportAction::Item(export::Item::Png) => self.export_png(ctx),
            ExportAction::Item(item) => self.export_files(Some(item)),
            ExportAction::All => {
                self.export_files(None);
                let summary = self.status.clone();
                if self.tree_rect.is_some() && self.tree_layout.is_some() {
                    self.export_png(ctx);
                    self.status = format!("{summary}; capturing tree.png ...");
                }
            }
        }
    }

    /// The folder the model came from: where a live run writes, or the output
    /// folder that was opened. Engine-written exports are copied from here.
    pub(crate) fn run_dir(&self) -> Option<PathBuf> {
        match &self.model.source {
            crate::state::Source::Folder(p) => Some(p.clone()),
            crate::state::Source::Live => {
                self.run_handle.as_ref().map(|h| h.outdir().to_path_buf())
            }
            crate::state::Source::Empty => None,
        }
    }

    /// The SVG looks like the tree panel: same theme, same decorations, a row
    /// height that is readable whatever the zoom on screen.
    pub(crate) fn svg_options(&self) -> crate::svg::SvgOptions {
        let dark = self.dark;
        let width = self
            .tree_rect
            .map(|(r, _)| (r.width() * 0.68).max(480.0))
            .unwrap_or(760.0);
        crate::svg::SvgOptions {
            dark,
            legend: self.tree.show_legend,
            branch_labels: self.tree.branch_labels,
            aligned_tips: self.tree.aligned_tips,
            color_strip: self.tree.color_strip,
            tree_width: width,
            row_px: self.tree.px_per_row.clamp(14.0, 24.0),
            selected_cluster: self.selected_cluster,
        }
    }

    /// Ask for a screenshot and write the tree panel to `path` when it arrives.
    pub(crate) fn export_png(&mut self, ctx: &Context) {
        let Some((rect, ppp)) = self.tree_rect else {
            self.status = "the tree panel has not been drawn yet".to_string();
            return;
        };
        let dir = PathBuf::from(self.export_text.trim());
        if self.export_text.trim().is_empty() {
            self.status = "type an export folder first".to_string();
            return;
        }
        let path = dir.join("tree.png");
        let slot = Arc::clone(&self.png_result);
        *slot.lock().unwrap_or_else(|e| e.into_inner()) = None;
        ctx.request_screenshot(move |img| {
            let res = export::write_png(&path, &img, Some((rect, ppp)))
                .map(|()| path.clone())
                .map_err(|e| e.to_string());
            if let Ok(mut g) = slot.lock() {
                *g = Some(res);
            }
        });
        self.status = "capturing the tree ...".to_string();
    }

    /// Collect the result of an earlier screenshot request.
    fn poll_png(&mut self) {
        let taken = self.png_result.lock().ok().and_then(|mut g| g.take());
        if let Some(res) = taken {
            match res {
                Ok(p) => self.status = format!("wrote {}", p.display()),
                Err(e) => {
                    self.status = format!("cannot write the picture: {e}");
                    self.model.push_local_note(self.status.clone());
                }
            }
        }
    }

    // --- pairwise ----------------------------------------------------------

    /// The pairwise readout for exactly two selected genomes.
    pub(crate) fn pair_view(&self) -> Option<PairView> {
        if self.selection.len() != 2 {
            return None;
        }
        let (a, b) = (self.selection.ids()[0], self.selection.ids()[1]);
        let ca = self.model.cluster_of(a)?;
        if Some(ca) != self.model.cluster_of(b) {
            return None;
        }
        if let Some(r) = &self.results {
            if let Some(m) = r.matrices.get(&ca) {
                let ia = m.ids.iter().position(|&x| x == a);
                let ib = m.ids.iter().position(|&x| x == b);
                if let (Some(ia), Some(ib)) = (ia, ib) {
                    let focus = self.model.cluster(ca).map(|c| c.focus).unwrap_or(2);
                    return Some(pairwise::pair_from_matrix(m, ia, ib, a, b, focus));
                }
            }
        }
        let (Some(ma), Some(mb)) = (self.model.map(ca, a), self.model.map(ca, b)) else {
            return None;
        };
        let focus = self
            .model
            .genome(a)
            .map(|g| g.focus)
            .unwrap_or(0)
            .min(self.model.genome(b).map(|g| g.focus).unwrap_or(0));
        Some(pairwise::pair_from_maps(ma, mb, a, b, focus))
    }

    /// What the header says about more than two selected genomes: clusters
    /// involved and the SNP range among those that share a cluster.
    pub(crate) fn selection_summary(&self) -> crate::selection::Summary {
        crate::selection::summarise(&self.model, self.results.as_ref(), self.selection.ids())
    }

    /// Why the pairwise readout is empty, in words.
    pub(crate) fn pair_hint(&self) -> &'static str {
        match self.selection.len() {
            0 => {
                "Select two genomes (Ctrl-click) to see their distance; Shift-click selects a \
                 range, Ctrl+A all."
            }
            1 => "Select one more genome (Ctrl-click) to see a pairwise distance.",
            _ => match (
                self.model.cluster_of(self.selection.ids()[0]),
                self.model.cluster_of(self.selection.ids()[1]),
            ) {
                (Some(a), Some(b)) if a != b => {
                    "The two genomes are in different clusters: no pairwise number exists."
                }
                _ => "Not measured yet: both genomes must reach Focus 2.",
            },
        }
    }

    // --- layout and animation ----------------------------------------------

    /// A cheap signature of what the tree panel shows: the set of opened
    /// clusters, or, in single-cluster mode, which cluster.
    pub(crate) fn expanded_signature(&self) -> u64 {
        let mut ids: Vec<u32> = self.expanded.iter().copied().collect();
        ids.sort_unstable();
        let mut h = 0x9E37_79B9_7F4A_7C15u64;
        for id in ids {
            h = h.rotate_left(5) ^ u64::from(id).wrapping_mul(0x100_0000_01b3);
        }
        if let Some(c) = self.single_cluster() {
            h = h.rotate_left(17) ^ (u64::from(c) + 1).wrapping_mul(0xff51_afd7_ed55_8ccd);
        }
        h
    }

    /// The cluster drawn on its own, when the tree panel is in that mode.
    pub(crate) fn single_cluster(&self) -> Option<u32> {
        match self.tree.mode {
            TreeMode::Cluster => self.selected_cluster,
            TreeMode::All => None,
        }
    }

    /// The display tree for the current mode.
    pub(crate) fn current_display_tree(&self) -> Option<layout::DisplayNode> {
        match self.single_cluster() {
            Some(c) => self.model.cluster_display_tree(c),
            None => self.model.display_tree(&self.expanded),
        }
    }

    /// Rebuild the tree layout only when the model or the opened set changed,
    /// and start the ~300 ms move animation from wherever the nodes are now.
    pub(crate) fn refresh_layout(&mut self, now: f64) {
        let sig = self.expanded_signature();
        if self.tree_layout.is_some()
            && self.laid_out_version == self.model.version
            && self.laid_out_expanded == sig
        {
            return;
        }
        let Some(display) = self.current_display_tree() else {
            self.tree_layout = None;
            self.positions.clear();
            self.laid_out_version = self.model.version;
            self.laid_out_expanded = sig;
            return;
        };
        let next = layout::layout(&display);
        // A tree that grew or shrank needs a refit; one that merely sharpened
        // keeps the view the user is already looking at.
        let rows_changed = match &self.tree_layout {
            Some(old) => (old.rows - next.rows).abs() > 1e-9,
            None => true,
        };
        // Chain an interrupted animation: start from the positions on screen,
        // not from the previous target.
        let from: Positions = match &self.tree_layout {
            Some(old) if !self.positions.is_empty() => {
                if self.anim_start.is_none() {
                    layout::snapshot(old)
                } else {
                    layout::positions_map(old, &self.positions)
                }
            }
            _ => Positions::new(),
        };
        let animate = !from.is_empty();
        self.positions = if animate {
            layout::animated(&next, &from, 0.0)
        } else {
            next.nodes.iter().map(|n| (n.x, n.y)).collect()
        };
        self.anim_from = from;
        self.anim_start = if animate { Some(now) } else { None };
        self.tree_layout = Some(next);
        self.laid_out_version = self.model.version;
        self.laid_out_expanded = sig;
        if rows_changed {
            self.tree.needs_fit = true;
        }
    }

    /// Advance the animation and keep the interpolated positions current.
    pub(crate) fn step_animation(&mut self, now: f64) {
        let Some(start) = self.anim_start else {
            return;
        };
        let Some(l) = &self.tree_layout else {
            self.anim_start = None;
            return;
        };
        let t = (now - start) / ANIM_SECS;
        if t >= 1.0 {
            self.positions = l.nodes.iter().map(|n| (n.x, n.y)).collect();
            self.anim_start = None;
        } else {
            self.positions = layout::animated(l, &self.anim_from, t);
        }
    }

    pub(crate) fn animating(&self) -> bool {
        self.anim_start.is_some()
    }

    // --- event plumbing ----------------------------------------------------

    /// Drain at most one batch of events into the model. Returning `true` means
    /// more are waiting and another frame should follow at once.
    pub(crate) fn drain_events(&mut self) -> bool {
        // Cloning the `Arc` rather than holding a reference, so that folding the
        // events can borrow the model mutably.
        let Some(pump) = self.run_handle.as_ref().map(|r| r.pump().clone()) else {
            return false;
        };
        for env in pump.drain(BATCH) {
            self.model.apply_envelope(env);
        }
        !pump.is_empty()
    }

    /// Collect the engine's results once it has finished.
    fn poll_run(&mut self) {
        let Some(handle) = &self.run_handle else {
            return;
        };
        if let Some(res) = handle.take_results() {
            match res {
                Ok(results) => {
                    self.status = if results.cancelled {
                        format!(
                            "stopped after {} of {} clusters",
                            self.model.exact_clusters(),
                            self.model.cluster_count()
                        )
                    } else {
                        format!(
                            "finished: {} genomes, {} clusters, output in {}",
                            self.model.genome_count(),
                            self.model.cluster_count(),
                            handle.outdir().display()
                        )
                    };
                    self.results = Some(results);
                }
                Err(e) => {
                    self.status = format!("the run failed: {e}");
                    self.model.push_local_note(e);
                    self.show_notes = true;
                }
            }
        }
    }

    /// Collect a folder that finished loading.
    fn poll_open(&mut self) {
        let Some(handle) = &self.opening else {
            return;
        };
        if !handle.is_done() {
            if let Some(p) = handle.progress() {
                self.status = format!("reading {} … {:.0} %", handle.dir().display(), p * 100.0);
            }
            return;
        }
        match handle.done() {
            Some(Ok(model)) => {
                let n = model.genome_count();
                let c = model.cluster_count();
                self.status = format!(
                    "opened {}: {n} genomes in {c} clusters",
                    handle.dir().display()
                );
                self.model = model;
                self.results = None;
                self.run_handle = None;
                self.selection.clear();
                self.expanded.clear();
                self.selected_cluster = None;
                self.reset_view();
                self.opening = None;
            }
            Some(Err(e)) => {
                self.status = e;
                self.opening = None;
            }
            None => {}
        }
    }

    /// Genomes dropped onto the window (section 8: no native file dialog).
    fn handle_dropped_files(&mut self, ctx: &Context) {
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if dropped.is_empty() {
            return;
        }
        let mut added = 0usize;
        for f in &dropped {
            let p = f.path().to_path_buf();
            if self.add_input_path(&p) {
                added += 1;
            }
        }
        if added > 0 {
            self.status = format!("added {added} path(s); {} in total", self.inputs.len());
        }
    }

    /// Add a file or a folder to the input list. Returns false when it was
    /// already there, so dropping the same folder twice does not double it.
    pub(crate) fn add_input_path(&mut self, path: &PathBuf) -> bool {
        if self.inputs.iter().any(|p| p == path) {
            return false;
        }
        self.inputs.push(path.clone());
        self.remember_input(&path.display().to_string());
        true
    }

    pub(crate) fn remove_input(&mut self, i: usize) {
        if i < self.inputs.len() {
            self.inputs.remove(i);
        }
    }
}

/// Command-line requests: open a results folder, or load genomes (and
/// optionally start the run) so that the window comes up ready.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Startup {
    pub open: Option<String>,
    pub inputs: Vec<PathBuf>,
    pub out: Option<String>,
    pub run: bool,
}

impl Startup {
    pub fn parse(args: &[String]) -> Result<Self, String> {
        let mut s = Startup::default();
        let mut it = args.iter();
        while let Some(a) = it.next() {
            match a.as_str() {
                "--open" => {
                    s.open = Some(it.next().ok_or("--open needs a folder")?.clone());
                }
                "--out" | "-o" => {
                    s.out = Some(it.next().ok_or("--out needs a folder")?.clone());
                }
                "--run" => s.run = true,
                x if x.starts_with('-') => return Err(format!("unknown option {x}")),
                x => s.inputs.push(PathBuf::from(x)),
            }
        }
        if s.open.is_some() && !s.inputs.is_empty() {
            return Err("give either --open or genome paths, not both".to_string());
        }
        if s.run && s.inputs.is_empty() {
            return Err("--run needs genome files or folders".to_string());
        }
        Ok(s)
    }
}

impl InfinityApp {
    /// Ctrl+A selects every genome the list shows, Escape clears the
    /// selection; neither while a text field has the keyboard.
    fn selection_keys(&mut self, ctx: &Context) {
        if ctx.egui_wants_keyboard_input() || self.outlier_dialog.is_some() {
            return;
        }
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::A)) {
            self.select_all_visible();
            self.status = format!("{} genomes selected", self.selection.len());
        }
        if !self.selection.is_empty()
            && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape))
        {
            self.clear_selection();
        }
    }

    fn apply_startup(&mut self, ctx: &Context) {
        if self.startup == Startup::default() {
            return;
        }
        let s = std::mem::take(&mut self.startup);
        if let Some(dir) = s.open {
            self.open_folder(&dir, ctx);
            return;
        }
        for p in &s.inputs {
            self.add_input_path(p);
        }
        if let Some(o) = s.out {
            self.outdir_text = o;
        }
        if s.run {
            self.start_run(ctx);
        }
    }
}

#[cfg(test)]
mod startup_tests {
    use super::Startup;

    fn v(a: &[&str]) -> Vec<String> {
        a.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn command_lines_parse() {
        assert_eq!(Startup::parse(&[]).unwrap(), Startup::default());
        let s = Startup::parse(&v(&["--run", "-o", "out", "a", "b"])).unwrap();
        assert!(s.run);
        assert_eq!(s.out.as_deref(), Some("out"));
        assert_eq!(s.inputs.len(), 2);
        assert_eq!(
            Startup::parse(&v(&["--open", "x"]))
                .unwrap()
                .open
                .as_deref(),
            Some("x")
        );
        assert!(Startup::parse(&v(&["--run"])).is_err());
        assert!(Startup::parse(&v(&["--open", "x", "a"])).is_err());
        assert!(Startup::parse(&v(&["--bogus"])).is_err());
    }
}

impl eframe::App for InfinityApp {
    /// Everything that is not drawing: draining the event stream, collecting
    /// finished work, and rebuilding the tree layout.
    ///
    /// eframe calls this before every [`Self::ui`], and also while the window is
    /// hidden, so a long run keeps being folded into the model even when the
    /// user has minimised the window.
    fn logic(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        self.apply_startup(ctx);
        if !self.theme_installed {
            crate::theme::install(ctx);
            ctx.set_theme(crate::theme::preference(&self.prefs.theme));
            self.theme_installed = true;
        }
        let now = ctx.input(|i| i.time);
        let more = self.drain_events();
        self.poll_outliers();
        self.poll_run();
        self.poll_restart(ctx);
        self.poll_open();
        self.poll_png();
        self.handle_dropped_files(ctx);
        self.refresh_layout(now);
        self.step_animation(now);

        // Repaint policy: immediately when events are still queued or the tree
        // is moving, on a heartbeat while a run is going (the relay thread is
        // the fast path), and not at all when idle.
        if more || self.animating() {
            ctx.request_repaint();
        } else if self.running() || self.opening.is_some() {
            ctx.request_repaint_after(HEARTBEAT);
        }
    }

    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        // Panels are shown in the order that carves the window up the way the
        // mock-up of section 8 lays it out: the full-width bars first, then the
        // two sides, then the map across what is left, and the tree in the
        // middle.
        //
        // The side panels are carved straight out of the root `Ui`, and the
        // tree and the map live inside a `CentralPanel`. That panel's child
        // `Ui` is bounded *and clipped* to what the side panels left over, so
        // a header that is wider than the space it gets can never spill over
        // a neighbouring panel: it wraps or is clipped instead. (Before, all
        // five were shown inside one `CentralPanel`, so the tree header and
        // canvas could widen past the right-hand panel and paint over it.)
        self.dark = ui.visuals().dark_mode;
        self.selection_keys(ui.ctx());
        crate::chrome::toolbar(ui, self);
        crate::chrome::status_bar(ui, self);
        let total_w = ui.available_width();
        crate::genome_panel::show(ui, self, total_w);
        crate::side_panel::show(ui, self, total_w);
        let p = crate::theme::palette(ui);
        egui::CentralPanel::default()
            .frame(
                egui::Frame::central_panel(ui.style())
                    .inner_margin(egui::Margin::same(6))
                    .fill(p.panel),
            )
            .show(ui, |ui| {
                crate::map_panel::show(ui, self);
                egui::CentralPanel::default()
                    .frame(egui::Frame::NONE)
                    .show(ui, |ui| crate::tree_panel::show(ui, self));
            });
        let ctx = ui.ctx().clone();
        crate::settings_panel::show(&ctx, self);
        crate::outlier_dialog::show(&ctx, self);
    }
}

/// Where a run writes by default: a folder next to the program's working
/// directory, because there is no file dialog to ask with.
fn default_outdir() -> String {
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    format!("infinity-out-{stamp}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout::NodeKind;
    use infinity_core::events::{ClusterInfo, Event, GenomeScannedInfo};

    fn app_with_two_clusters() -> InfinityApp {
        let mut app = InfinityApp::default();
        for i in 0..4u32 {
            app.model.apply(
                0.0,
                Event::GenomeScanned(GenomeScannedInfo {
                    id: i,
                    name: format!("g{i}"),
                    length: 1000,
                    contigs: 1,
                    n50: 1000,
                    gc: 0.5,
                    sketch_size: 10,
                }),
            );
        }
        app.model.apply(
            0.0,
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![
                    ClusterInfo {
                        id: 0,
                        pivot: 0,
                        pivot_name: "g0".into(),
                        members: vec![0, 1],
                        member_names: vec!["g0".into(), "g1".into()],
                        focus: 0,
                        cluster_dist: 0.005,
                    },
                    ClusterInfo {
                        id: 1,
                        pivot: 2,
                        pivot_name: "g2".into(),
                        members: vec![2, 3],
                        member_names: vec!["g2".into(), "g3".into()],
                        focus: 0,
                        cluster_dist: 0.005,
                    },
                ],
            },
        );
        app.model.apply(
            0.0,
            Event::BackboneTree {
                newick: "(g0:0.001,g2:0.002);".into(),
                pivots: vec!["g0".into(), "g2".into()],
            },
        );
        app
    }

    #[test]
    fn a_fresh_app_has_valid_default_parameters_and_nothing_to_run() {
        let app = InfinityApp::new();
        assert!(app.params_error.is_none());
        assert!(!app.running());
        assert!(app.inputs.is_empty());
        assert!(app.pair_view().is_none());
        assert!(app.pair_hint().contains("Select two"));
    }

    #[test]
    fn selection_holds_any_number_of_genomes() {
        let mut app = app_with_two_clusters();
        app.select_genome(0, false);
        app.select_genome(1, true);
        assert_eq!(app.selection, vec![0, 1]);
        app.select_genome(2, true);
        app.select_genome(3, true);
        assert_eq!(app.selection, vec![0, 1, 2, 3], "nothing drops off");
        // Clicking a selected genome again deselects it.
        app.select_genome(2, true);
        assert_eq!(app.selection, vec![0, 1, 3]);
        // A plain click starts again.
        app.select_genome(2, false);
        assert_eq!(app.selection, vec![2]);
        app.clear_selection();
        assert!(app.selection.is_empty());
    }

    #[test]
    fn shift_click_in_the_list_follows_the_lists_order_and_filter() {
        let mut app = app_with_two_clusters();
        // Sorted by name descending: g3, g2, g1, g0.
        app.sort = SortSpec {
            column: crate::state::SortColumn::Name,
            ascending: false,
        };
        let order = app.list_order();
        assert_eq!(order, vec![3, 2, 1, 0]);
        app.click_genome(2, ClickMods::PLAIN, &order);
        app.click_genome(0, ClickMods::RANGE, &order);
        assert_eq!(app.selection, vec![2, 1, 0]);
        // Ctrl+A takes what the filter shows, nothing more.
        app.filter = "g3".into();
        app.select_all_visible();
        assert_eq!(app.selection, vec![3]);
        app.filter.clear();
        app.select_all_visible();
        assert_eq!(app.selection.len(), 4);
    }

    #[test]
    fn selected_genomes_in_collapsed_clusters_are_counted_and_can_be_opened() {
        let mut app = app_with_two_clusters();
        app.model.apply(
            0.0,
            Event::ClusterFocus(infinity_core::events::ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,g1:2);".into(),
                members: vec![],
            }),
        );
        app.select_genome(0, false);
        app.select_genome(1, true);
        app.select_genome(3, true);
        assert_eq!(app.hidden_selected(), vec![(0, 2), (1, 1)]);
        app.expand_selected();
        // Cluster 0 has a tree and opens; cluster 1 has none yet and stays a
        // triangle, so its selected genome is still counted as hidden.
        assert_eq!(app.hidden_selected(), vec![(1, 1)]);
        assert!(app.expanded.contains(&0));
    }

    fn outliers_event() -> Event {
        Event::Outliers(vec![infinity_core::events::OutlierInfo {
            id: 3,
            name: "g3".into(),
            length: 1000,
            gc: 0.65,
            nearest_dist: 1.0,
            nearest_name: String::new(),
            reasons: vec!["no close relative".into()],
        }])
    }

    #[test]
    fn flagged_genomes_open_the_dialog_once_per_live_run() {
        let mut app = app_with_two_clusters();
        app.model.source = Source::Live;
        app.poll_outliers();
        assert!(app.outlier_dialog.is_none(), "nothing flagged yet");
        app.model.apply(0.0, outliers_event());
        app.poll_outliers();
        let d = app.outlier_dialog.as_ref().expect("the dialog opens");
        assert_eq!(d.ticked_names(), vec!["g3"]);
        // Closed ("Continue anyway"), it does not come back by itself.
        app.outlier_dialog = None;
        app.outliers_acknowledged = true;
        app.poll_outliers();
        assert!(app.outlier_dialog.is_none());
        // The warning is in the notes, where the Run panel lists it.
        assert!(app
            .model
            .notes
            .iter()
            .any(|n| n.text.contains("divergent assembl") && n.text.contains("g3")));
    }

    #[test]
    fn an_opened_folder_never_pops_the_dialog() {
        let mut app = app_with_two_clusters();
        app.model.source = Source::Folder(PathBuf::from("/results"));
        app.model.apply(0.0, outliers_event());
        app.poll_outliers();
        assert!(app.outlier_dialog.is_none());
        assert_eq!(app.model.outliers.len(), 1, "but the warning stays");
    }

    #[test]
    fn a_restart_waits_for_the_cancelled_run_and_nothing_else() {
        let mut app = app_with_two_clusters();
        // Nothing is running, so a pending restart starts at once; with no
        // genomes to compare it reports that instead of running.
        app.restart_pending = true;
        app.poll_restart(&Context::default());
        assert!(!app.restart_pending);
        assert!(app.status.contains("nothing to compare") || app.params_error.is_some());
    }

    #[test]
    fn more_than_two_genomes_get_a_summary_instead_of_a_pair() {
        let mut app = app_with_two_clusters();
        for g in [0, 1, 2] {
            app.select_genome(g, true);
        }
        assert!(app.pair_view().is_none());
        let s = app.selection_summary();
        assert_eq!(s.n, 3);
        assert_eq!(s.clusters, vec![(0, 2), (1, 1)]);
        assert_eq!(s.pairs, 1);
    }

    #[test]
    fn selecting_a_genome_selects_its_cluster_for_the_map() {
        let mut app = app_with_two_clusters();
        app.select_genome(3, false);
        assert_eq!(app.selected_cluster, Some(1));
    }

    #[test]
    fn opening_a_cluster_marks_it_for_a_layout_rebuild() {
        let mut app = app_with_two_clusters();
        app.refresh_layout(0.0);
        let sig = app.expanded_signature();
        app.toggle_cluster(0);
        assert_ne!(app.expanded_signature(), sig);
        // Toggling twice is back to the same signature.
        app.toggle_cluster(0);
        assert_eq!(app.expanded_signature(), sig);
        assert!(app.expanded.is_empty());
    }

    #[test]
    fn the_layout_is_rebuilt_only_when_something_changed() {
        let mut app = app_with_two_clusters();
        app.refresh_layout(0.0);
        assert!(app.tree_layout.is_some());
        assert!(app.anim_start.is_none(), "a first layout does not animate");
        let first = app.tree_layout.as_ref().unwrap().nodes.len();
        app.refresh_layout(1.0);
        assert_eq!(app.tree_layout.as_ref().unwrap().nodes.len(), first);
        assert!(app.anim_start.is_none());
        // A new tree for one cluster changes the model version, and now the
        // nodes move.
        app.model.apply(
            2.0,
            Event::ClusterFocus(infinity_core::events::ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,g1:2);".into(),
                members: vec![],
            }),
        );
        app.expanded.insert(0);
        app.refresh_layout(3.0);
        assert!(app.anim_start.is_some(), "a changed tree animates");
        app.step_animation(3.0 + ANIM_SECS + 0.01);
        assert!(app.anim_start.is_none(), "and stops once it has arrived");
        let l = app.tree_layout.as_ref().unwrap();
        for (n, p) in l.nodes.iter().zip(app.positions.iter()) {
            assert!((n.x - p.0).abs() < 1e-9 && (n.y - p.1).abs() < 1e-9);
        }
    }

    #[test]
    fn the_animation_moves_monotonically_towards_the_target() {
        let mut app = app_with_two_clusters();
        app.refresh_layout(0.0);
        app.expanded.insert(0);
        app.model.apply(
            0.0,
            Event::ClusterFocus(infinity_core::events::ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,g1:2);".into(),
                members: vec![],
            }),
        );
        app.refresh_layout(1.0);
        assert!(app.animating());
        let l = app.tree_layout.clone().unwrap();
        let mut last = app.positions.iter().map(|p| p.1).collect::<Vec<_>>();
        for i in 1..=10 {
            app.step_animation(1.0 + ANIM_SECS * i as f64 / 10.0);
            let cur: Vec<f64> = app.positions.iter().map(|p| p.1).collect();
            for (j, (c, target)) in cur.iter().zip(l.nodes.iter().map(|n| n.y)).enumerate() {
                let before = last[j];
                // Each step is closer to (or at) the target than the one before.
                assert!(
                    (c - target).abs() <= (before - target).abs() + 1e-9,
                    "row {j} moved away from its target"
                );
            }
            last = cur;
        }
        assert!(!app.animating());
    }

    #[test]
    fn a_dropped_path_is_added_once() {
        let mut app = InfinityApp::default();
        let p = PathBuf::from("/data/genomes");
        assert!(app.add_input_path(&p));
        assert!(!app.add_input_path(&p));
        assert_eq!(app.inputs, vec![p]);
        app.remove_input(0);
        assert!(app.inputs.is_empty());
        // Removing something that is not there is not an error.
        app.remove_input(5);
        assert!(app.inputs.is_empty());
    }

    #[test]
    fn the_display_tree_follows_the_expanded_set() {
        let app = app_with_two_clusters();
        let d = app.model.display_tree(&HashSet::new()).unwrap();
        assert!(d
            .children
            .iter()
            .all(|c| matches!(c.kind, NodeKind::ClusterCollapsed { .. })));
        let mut open = HashSet::new();
        open.insert(0u32);
        // Cluster 0 has no tree of its own yet, so it stays collapsed.
        let d2 = app.model.display_tree(&open).unwrap();
        assert!(d2
            .children
            .iter()
            .all(|c| matches!(c.kind, NodeKind::ClusterCollapsed { .. })));
    }

    #[test]
    fn pair_view_is_empty_until_two_genomes_of_one_cluster_are_chosen() {
        let mut app = app_with_two_clusters();
        assert!(app.pair_view().is_none());
        app.select_genome(0, false);
        assert!(app.pair_view().is_none());
        app.select_genome(2, true);
        assert!(app.pair_view().is_none(), "different clusters");
        assert!(app.pair_hint().contains("different clusters"));
        app.select_genome(2, true);
        app.select_genome(1, true);
        assert!(app.pair_view().is_none(), "no map data yet");
        assert!(app.pair_hint().contains("Focus 2"));
    }
}
