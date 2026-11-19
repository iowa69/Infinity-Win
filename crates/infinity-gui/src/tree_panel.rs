//! The tree panel (section 8): the tree as a clean rectangular phylogram, with
//! pan and zoom, the SNPs of each genome's own branch after its name, SNP
//! groups shaded at a threshold the user slides, collapsed clusters that open
//! on click, a scale bar and a legend of what is on screen, hover numbers and
//! a search box.
//!
//! The picture itself is built by [`crate::pubtree`], the same builder the
//! tree window and the SVG export use; this file turns the pan and zoom into a
//! [`Geometry`], paints the scene, and turns the scene's click targets back
//! into selections. The transform is a small pure type ([`View`]) so that the
//! mapping can be tested without a window.

use eframe::egui;
use egui::{Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Vec2};

use crate::app::{InfinityApp, TreeMode};
use crate::layout::{NodeKind, TreeLayout};
use crate::pairwise::PairSource;
use crate::pubtree::{self, Geometry, Hit, Ink};
use crate::selection::ClickMods;
use crate::state::UiModel;
use crate::{colors, fmt, groups, theme};

const MARGIN: f32 = 12.0;
const MIN_PX_PER_ROW: f32 = 0.35;
const MAX_PX_PER_ROW: f32 = 240.0;
const MIN_PX_PER_UNIT: f32 = 8.0;
const MAX_PX_PER_UNIT: f32 = 40_000.0;
/// Rows are never fitted taller than this: a three-genome tree should not
/// fill the panel with 200-pixel rows.
const MAX_FIT_ROW: f32 = 24.0;

/// Layout coordinates to screen coordinates and back.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct View {
    pub origin: Pos2,
    pub px_per_row: f32,
    pub px_per_unit: f32,
}

impl View {
    pub fn screen(&self, p: (f64, f64)) -> Pos2 {
        Pos2::new(
            self.origin.x + p.0 as f32 * self.px_per_unit,
            self.origin.y + p.1 as f32 * self.px_per_row,
        )
    }

    pub fn layout(&self, s: Pos2) -> (f64, f64) {
        (
            f64::from((s.x - self.origin.x) / self.px_per_unit),
            f64::from((s.y - self.origin.y) / self.px_per_row),
        )
    }

    /// Scale both axes by `factor` about a fixed screen point, clamped so that a
    /// runaway wheel event cannot lose the tree.
    pub fn zoomed_about(self, pointer: Pos2, factor: f32) -> Self {
        let here = self.layout(pointer);
        let mut v = self;
        v.px_per_row = (v.px_per_row * factor).clamp(MIN_PX_PER_ROW, MAX_PX_PER_ROW);
        v.px_per_unit = (v.px_per_unit * factor).clamp(MIN_PX_PER_UNIT, MAX_PX_PER_UNIT);
        let now = v.screen(here);
        v.origin += pointer - now;
        v
    }

    fn geometry(&self, clip: Rect) -> Geometry {
        Geometry {
            origin: self.origin,
            px_per_row: self.px_per_row,
            px_per_unit: self.px_per_unit,
            clip,
        }
    }
}

/// The genome ids in the tree's row order, for a Shift-click range. A
/// collapsed cluster stands for all of its members, in the order the cluster
/// lists them, so a range drawn across a wedge includes what it hides.
pub(crate) fn row_order(l: &TreeLayout, model: &UiModel) -> Vec<u32> {
    let mut out = Vec::with_capacity(l.leaves.len());
    for &i in &l.leaves {
        match &l.nodes[i].kind {
            NodeKind::Genome { id, .. } if *id != u32::MAX => out.push(*id),
            NodeKind::ClusterCollapsed { id, .. } => {
                if let Some(c) = model.cluster(*id) {
                    out.extend(c.members.iter().copied());
                }
            }
            _ => {}
        }
    }
    out
}

pub fn show(ui: &mut egui::Ui, app: &mut InfinityApp) {
    header(ui, app);
    canvas(ui, app);
    let ctx = ui.ctx().clone();
    help_window(&ctx, app);
}

// ---------------------------------------------------------------------------
// Header
// ---------------------------------------------------------------------------

fn header(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let mut fit_now = false;
    let mut zoom_by: Option<f32> = None;
    let mut picked: Option<u32> = None;
    let mut open_window = false;

    // Wrapped rather than one long row: at the smallest window the controls
    // flow onto a second line instead of running under the side panel.
    theme::panel_header(ui, "Tree", |ui| {
        let sel = app.selected_cluster;
        let label = match sel.and_then(|c| app.model.cluster(c)) {
            Some(c) => format!(
                "c{} · {} genomes · {}",
                c.id,
                c.size(),
                state_word(c.focus, c.approximate)
            ),
            None => format!("{} clusters", app.model.cluster_count()),
        };
        let ids: Vec<u32> = app.model.clusters.iter().map(|c| c.id).collect();
        egui::ComboBox::from_id_salt("cluster-select")
            .selected_text(label)
            .show_ui(ui, |ui| {
                for id in ids {
                    let text = app
                        .model
                        .cluster(id)
                        .map(|c| {
                            format!(
                                "cluster {} — {} genomes, {}",
                                c.id,
                                c.size(),
                                state_word(c.focus, c.approximate)
                            )
                        })
                        .unwrap_or_default();
                    if ui.selectable_label(sel == Some(id), text).clicked() {
                        picked = Some(id);
                    }
                }
            });

        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.tree.search)
                .desired_width(120.0)
                .hint_text("🔍 find a genome"),
        );
        if resp.changed() {
            app.tree.search_hit = app.model.search(&app.tree.search);
            app.filter = app.tree.search.clone();
        }
        if let Some(id) = app.tree.search_hit {
            let name = app
                .model
                .genome(id)
                .map(|g| g.name.clone())
                .unwrap_or_default();
            let go = ui
                .button("go to")
                .on_hover_text(format!("select {name} and centre the tree on it"))
                .clicked()
                || (resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)));
            if go {
                app.select_genome(id, false);
                centre_on(app, id);
            }
        }

        ui.separator();
        // All clusters, or the selected one alone as a standard phylogram.
        let can_single = app.selected_cluster.is_some();
        let before = app.tree.mode;
        ui.selectable_value(&mut app.tree.mode, TreeMode::All, "All")
            .on_hover_text("the backbone with every cluster, collapsed or opened");
        ui.add_enabled_ui(can_single, |ui| {
            ui.selectable_value(&mut app.tree.mode, TreeMode::Cluster, "Cluster")
                .on_hover_text("the selected cluster on its own, branch lengths in SNPs")
                .on_disabled_hover_text("select a cluster first");
        });
        if app.tree.mode != before {
            app.tree.needs_fit = true;
            app.tree.user_moved = false;
        }

        ui.separator();
        if ui.button("−").on_hover_text("zoom out").clicked() {
            zoom_by = Some(1.0 / 1.25);
        }
        if ui.button("+").on_hover_text("zoom in").clicked() {
            zoom_by = Some(1.25);
        }
        if ui
            .button("Fit")
            .on_hover_text("fit the whole tree")
            .clicked()
        {
            fit_now = true;
        }
        ui.menu_button("View ▾", |ui| {
            let o = &mut app.tree.opts;
            ui.checkbox(&mut o.tip_snps, "SNPs at the tips (+n)")
                .on_hover_text("after each name, the SNPs on that genome's own branch");
            ui.checkbox(&mut o.branch_labels, "SNP counts on internal branches");
            ui.checkbox(&mut o.cluster_brackets, "Cluster brackets");
            ui.checkbox(&mut o.aligned_tips, "Aligned tip labels");
            ui.checkbox(&mut o.color_strip, "Cluster colour strip");
            ui.checkbox(&mut o.legend, "Legend");
            ui.separator();
            ui.checkbox(&mut app.tree.paper, "White background")
                .on_hover_text("draw the tree on white, like the exported figure");
        });
        if ui
            .button("⛶ Tree window")
            .on_hover_text(
                "the whole tree in a window of its own, on white, ready for a figure: \
                 pan, zoom, save PNG or SVG",
            )
            .clicked()
        {
            open_window = true;
        }
        if ui
            .selectable_label(app.tree.show_help, "？ How to read")
            .on_hover_text("what the branches, numbers and shaded groups mean")
            .clicked()
        {
            app.tree.show_help = !app.tree.show_help;
        }
    });
    if open_window {
        let ctx = ui.ctx().clone();
        crate::tree_window::open(app, &ctx);
    }
    group_controls(ui, app);
    // The pairwise readout gets a line of its own: it is the longest thing in
    // the header and the one that changes most.
    let mut open_hidden = false;
    ui.horizontal_wrapped(|ui| {
        pair_readout(ui, app);
        let hidden = app.hidden_selected();
        if !hidden.is_empty() {
            let n: usize = hidden.iter().map(|h| h.1).sum();
            let p = theme::palette(ui);
            if ui
                .small_button(
                    RichText::new(format!(
                        "⤢ show {n} selected in {} collapsed cluster{}",
                        hidden.len(),
                        if hidden.len() == 1 { "" } else { "s" }
                    ))
                    .color(p.accent),
                )
                .on_hover_text("open the clusters that hold selected genomes")
                .clicked()
            {
                open_hidden = true;
            }
        }
    });
    if open_hidden {
        app.expand_selected();
        for &(c, _) in &app.hidden_selected() {
            // Clusters without a tree yet cannot open; ask for them first.
            app.set_cluster(c, true);
        }
    }

    if let Some(f) = zoom_by {
        app.tree.zoom_request = Some(f);
    }
    if fit_now {
        app.tree.needs_fit = true;
        app.tree.user_moved = false;
    }
    if let Some(id) = picked {
        app.set_cluster(id, true);
    }
}

/// "exact" or "≈ estimate".
fn state_word(focus: u8, approximate: bool) -> &'static str {
    if focus >= fmt::FOCUS_EXACT && !approximate {
        "exact"
    } else {
        "≈ estimate"
    }
}

/// The SNP-threshold row: on/off, a slider that is fine at the low end, a
/// number box, quick picks, and what it found.
fn group_controls(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let p = theme::palette(ui);
    let mut changed = false;
    ui.horizontal_wrapped(|ui| {
        changed |= ui
            .checkbox(&mut app.tree.groups_on, "SNP groups")
            .on_hover_text(
                "Shade the genomes linked by chains of at most this many SNPs \
                 (single linkage), within each exact cluster.",
            )
            .changed();
        ui.add_enabled_ui(app.tree.groups_on, |ui| {
            ui.label(RichText::new("≤").color(p.text_weak));
            let mut t = groups::snps_to_slider(app.tree.group_n);
            ui.spacing_mut().slider_width = 150.0;
            if ui
                .add(egui::Slider::new(&mut t, 0.0..=1.0).show_value(false))
                .on_hover_text(
                    "0–20 SNPs one by one on the left part of the track, up to 100 on the right",
                )
                .changed()
            {
                app.tree.group_n = groups::slider_to_snps(t);
                changed = true;
            }
            changed |= ui
                .add(
                    egui::DragValue::new(&mut app.tree.group_n)
                        .range(0..=100)
                        .suffix(" SNPs"),
                )
                .changed();
            for q in [5u32, 10, 20] {
                if ui
                    .selectable_label(app.tree.group_n == q, q.to_string())
                    .on_hover_text(format!("{q} SNPs"))
                    .clicked()
                {
                    app.tree.group_n = q;
                    changed = true;
                }
            }
        });
        if changed {
            app.refresh_groups();
        }
        if app.tree.groups_on {
            let s = &app.groups.set;
            let text = if s.covered.is_empty() {
                if app.model.exact_clusters() == 0 || app.model.genome_count() == 0 {
                    "groups appear when a cluster becomes exact".to_string()
                } else {
                    "no SNP matrix for the exact clusters yet".to_string()
                }
            } else if s.groups.is_empty() {
                format!("no two genomes within {} SNPs", s.threshold)
            } else {
                format!(
                    "{} group{} · {} genomes",
                    s.groups.len(),
                    if s.groups.len() == 1 { "" } else { "s" },
                    s.grouped_genomes()
                )
            };
            ui.label(RichText::new(text).small().color(p.text_weak))
                .on_hover_text("the list is in the Run panel; click a group label to select it");
        }
    });
}

/// The pairwise readout of section 8: distance, SNP count and focus of the two
/// selected leaves.
fn pair_readout(ui: &mut egui::Ui, app: &InfinityApp) {
    match app.pair_view() {
        Some(p) => {
            let na = app
                .model
                .genome(p.a)
                .map(|g| g.name.as_str())
                .unwrap_or("?");
            let nb = app
                .model
                .genome(p.b)
                .map(|g| g.name.as_str())
                .unwrap_or("?");
            let nums = format!(
                "{} SNPs | indels {}",
                p.snps_text(),
                if p.available() {
                    fmt::count(p.indels)
                } else {
                    "n/a".to_string()
                }
            );
            ui.label(
                RichText::new(format!(
                    "{na} × {nb}:  {nums} | d = {} | {}",
                    p.dist_text(),
                    fmt::focus_label(p.focus, p.approximate)
                ))
                .small()
                .color(colors::focus_color(p.focus, ui.visuals().dark_mode)),
            )
            .on_hover_text(match p.source {
                PairSource::Matrix => "taken from the engine's own distance matrix",
                PairSource::Variants => {
                    "recomputed from the two members' variant lists against the pivot (6.9)"
                }
                PairSource::Unavailable => app.pair_hint(),
            });
        }
        None if app.selection.len() > 2 => {
            let s = app.selection_summary();
            let approx = s.snps.is_some_and(|r| r.approximate);
            let focus = if approx {
                fmt::FOCUS_BLOCKS
            } else {
                fmt::FOCUS_EXACT
            };
            let colour = if s.snps.is_some() {
                colors::focus_color(focus, ui.visuals().dark_mode)
            } else {
                theme::palette(ui).text
            };
            ui.label(RichText::new(s.text()).small().color(colour))
                .on_hover_text(
                    "SNP distances are counted only between genomes of the same cluster; \
                     select exactly two for their full pairwise numbers.",
                );
        }
        None => {
            theme::hint(ui, app.pair_hint());
        }
    }
}

/// The "How to read this tree" panel.
fn help_window(ctx: &egui::Context, app: &mut InfinityApp) {
    if !app.tree.show_help {
        return;
    }
    let n = app.tree.group_n;
    let mut open = true;
    egui::Window::new("How to read this tree")
        .open(&mut open)
        .collapsible(false)
        .resizable(false)
        .default_width(400.0)
        .pivot(Align2::RIGHT_TOP)
        .default_pos(ctx.content_rect().right_top() + Vec2::new(-290.0, 110.0))
        .show(ctx, |ui| {
            let p = theme::palette(ui);
            let item = |ui: &mut egui::Ui, head: &str, body: String| {
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new(head).strong().color(p.text_strong));
                    ui.label(RichText::new(body).color(p.text));
                });
                ui.add_space(3.0);
            };
            item(
                ui,
                "Branch length = SNPs.",
                "The further right a genome sits, the more SNPs separate it from the others; \
                 the scale bar at the bottom gives the length of 1, 2, 5 or 10 SNPs."
                    .into(),
            );
            item(
                ui,
                "+n after a name",
                "is the number of SNPs on that genome's own branch: differences it has that its \
                 closest relatives in the tree do not. +0 means none of its own."
                    .into(),
            );
            item(
                ui,
                "Shaded groups (G1, G2 …)",
                format!(
                    "are genomes linked by chains of at most {n} SNPs (single linkage), counted \
                     from the exact pairwise SNP matrix. Move the slider to see what groups at \
                     5, 10 or 20 SNPs; click a group's label to select its genomes."
                ),
            );
            item(
                ui,
                "Dark branches",
                "are exact SNP counts inside a cluster. Grey branches are estimates: the \
                 distances between clusters, or a cluster still being refined (≈)."
                    .into(),
            );
            item(
                ui,
                "Triangles",
                "are collapsed clusters: click one to open it. A bracket on the right names \
                 each opened cluster."
                    .into(),
            );
            ui.add_space(2.0);
            theme::hint(
                ui,
                "SNPs are counted only in the core genome every member of a cluster shares, \
                 after the SNP rules (Settings → SNP rules).",
            );
        });
    if !open {
        app.tree.show_help = false;
    }
}

/// Move the view so that a genome's leaf ends up in the middle of the panel. If
/// the genome is inside a collapsed cluster, that cluster is opened instead and
/// the centring happens on the next frame, once it has a leaf.
fn centre_on(app: &mut InfinityApp, genome: u32) {
    let Some(name) = app.model.genome(genome).map(|g| g.name.clone()) else {
        return;
    };
    let coords = app
        .tree_layout
        .as_ref()
        .and_then(|l| l.find(&name).map(|i| (l.nodes[i].x, l.nodes[i].y)));
    match coords {
        Some(c) => app.tree.centre = Some(c),
        None => {
            if let Some(cid) = app.model.cluster_of(genome) {
                app.expanded.insert(cid);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Canvas
// ---------------------------------------------------------------------------

fn canvas(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let rect = remaining_rect(ui);
    let resp = ui.allocate_rect(rect, Sense::click_and_drag());
    app.tree_rect = Some((rect, ui.ctx().pixels_per_point()));
    let pal = theme::palette(ui);
    let ink = if app.tree.paper {
        Ink::PAPER
    } else {
        Ink::screen(&pal)
    };
    let painter = ui.painter().with_clip_rect(rect);
    painter.rect_filled(rect, 6, ink.paper);
    painter.rect_stroke(
        rect,
        6,
        Stroke::new(1.0, pal.border),
        egui::StrokeKind::Inside,
    );

    // Pointer state is read once, up front: `Response` is consumed by the
    // tooltip call at the end of the frame.
    let pointer = resp.hover_pos();
    let dragged_by = if resp.dragged() {
        Some(resp.drag_delta())
    } else {
        None
    };
    let clicked = resp.clicked();
    let mods = ui.ctx().input(|i| ClickMods::from_modifiers(&i.modifiers));
    let scroll = ui.ctx().input(|i| i.smooth_scroll_delta);

    // Taken out for the frame (and put back below) rather than cloned: the
    // layout of a thousand genomes is not small.
    let Some(l) = app.tree_layout.take() else {
        hint(
            &painter,
            rect,
            if app.model.genome_count() == 0 {
                "No genomes yet.\nAdd files or a folder, then press Run."
            } else if app.model.backbone.is_none() {
                "Scanning genomes …\nThe first tree appears as soon as clustering is done."
            } else {
                "No tree to draw."
            },
            ink.text_weak,
        );
        return;
    };
    let opts = app.tree_opts();
    let marks = app.marks(&l, true);

    // --- interaction first, so this frame draws the result ------------------
    if app.tree.needs_fit && !app.tree.user_moved {
        fit_view(app, rect, &l, &marks, &opts);
        app.tree.needs_fit = false;
    }
    if let Some(d) = dragged_by {
        app.tree.offset += d;
        app.tree.user_moved = true;
    }
    if scroll.y.abs() > 0.01 {
        let anchor = pointer.unwrap_or_else(|| rect.center());
        let v = current_view(app, rect).zoomed_about(anchor, (scroll.y * 0.004).exp());
        apply_view(app, rect, v);
        app.tree.user_moved = true;
    }
    if scroll.x.abs() > 0.01 {
        app.tree.offset.x += scroll.x;
        app.tree.user_moved = true;
    }
    if let Some(f) = app.tree.zoom_request.take() {
        let v = current_view(app, rect).zoomed_about(rect.center(), f);
        apply_view(app, rect, v);
        app.tree.user_moved = true;
    }
    if let Some(c) = app.tree.centre.take() {
        let mut v = current_view(app, rect);
        v.origin = Pos2::new(
            rect.center().x - c.0 as f32 * v.px_per_unit,
            rect.center().y - c.1 as f32 * v.px_per_row,
        );
        apply_view(app, rect, v);
        app.tree.user_moved = true;
    }

    // --- draw ----------------------------------------------------------------
    let g = current_view(app, rect).geometry(rect.shrink(1.0));
    let positions = (app.positions.len() == l.nodes.len()).then_some(app.positions.as_slice());
    let sc = pubtree::scene(&l, &marks, &g, &opts, &ink, positions);
    pubtree::paint(&painter, &sc.prims);
    let fh = pubtree::footer_height(&l, &marks, &g, &sc.used, opts.legend, rect.width());
    let footer = Rect::from_min_max(
        Pos2::new(rect.left() + 1.0, rect.bottom() - fh - 4.0),
        Pos2::new(rect.right() - 1.0, rect.bottom() - 2.0),
    );
    if fh > 0.0 {
        pubtree::paint(
            &painter,
            &pubtree::footer(&l, &marks, &g, &sc.used, opts.legend, &ink, footer, true),
        );
    }

    // --- hover, tooltip, click -------------------------------------------------
    let hover = pointer
        .filter(|p| fh <= 0.0 || !footer.contains(*p))
        .and_then(|p| {
            sc.hits
                .iter()
                .rev()
                .find(|(_, r)| r.contains(p))
                .map(|(h, _)| *h)
        });
    let tooltip: Vec<String> = match hover {
        Some(Hit::Genome(id)) => genome_tooltip(&app.model, id),
        Some(Hit::Cluster(id)) => cluster_tooltip(&app.model, id, app.expanded.contains(&id)),
        Some(Hit::Badge(id)) => vec![
            format!(
                "{} selected genome(s) inside cluster {id}",
                marks.selected_in.get(&id).copied().unwrap_or(0)
            ),
            "click to open the cluster and show them".to_string(),
        ],
        Some(Hit::Group(i)) => group_tooltip(app, i),
        None => Vec::new(),
    };
    if clicked {
        match hover {
            Some(Hit::Genome(id)) if id != u32::MAX => {
                let order = row_order(&l, &app.model);
                app.click_genome(id, mods, &order);
            }
            Some(Hit::Cluster(id)) => app.toggle_cluster(id),
            Some(Hit::Badge(id)) => {
                app.expanded.insert(id);
                app.set_cluster(id, true);
            }
            Some(Hit::Group(i)) => app.select_group(i, mods.toggle),
            // A click on empty canvas keeps the selection: a long hand-picked
            // selection is too easy to lose that way. Escape clears it.
            _ => {}
        }
    }
    ui.ctx().output_mut(|o| {
        o.cursor_icon = if hover.is_some() {
            egui::CursorIcon::PointingHand
        } else if resp.dragged() {
            egui::CursorIcon::Grabbing
        } else {
            egui::CursorIcon::Default
        };
    });
    app.tree_layout = Some(l);
    if !tooltip.is_empty() {
        resp.on_hover_ui(|ui| {
            for line in tooltip {
                ui.label(line);
            }
        });
    }
}

/// What is left of the panel below the header, never wider than the panel
/// itself. `available_size_before_wrap` can report more than the panel when a
/// header row overflowed, and a canvas sized from it painted over the
/// neighbouring panel.
pub(crate) fn remaining_rect(ui: &egui::Ui) -> Rect {
    let avail = ui.available_rect_before_wrap();
    let bound = ui.max_rect();
    let r = Rect::from_min_max(
        avail.min,
        Pos2::new(
            avail.right().min(bound.right()),
            avail.bottom().min(bound.bottom()),
        ),
    );
    if r.width() < 40.0 || r.height() < 40.0 {
        Rect::from_min_size(r.min, r.size().max(Vec2::splat(40.0)))
    } else {
        r
    }
}

fn hint(painter: &egui::Painter, rect: Rect, msg: &str, colour: Color32) {
    painter.text(
        rect.center(),
        Align2::CENTER_CENTER,
        msg,
        FontId::proportional(15.0),
        colour,
    );
}

/// The view: the first row's centre sits half a row below the top margin, so
/// its label is never cut.
fn current_view(app: &InfinityApp, rect: Rect) -> View {
    View {
        origin: rect.min + Vec2::new(MARGIN, MARGIN + app.tree.px_per_row * 0.5) + app.tree.offset,
        px_per_row: app.tree.px_per_row,
        px_per_unit: app.tree.px_per_unit,
    }
}

/// Store a view back as a pan offset from the panel corner, so that resizing the
/// window keeps the tree where it was.
fn apply_view(app: &mut InfinityApp, rect: Rect, v: View) {
    app.tree.offset = v.origin - rect.min - Vec2::new(MARGIN, MARGIN + v.px_per_row * 0.5);
    app.tree.px_per_row = v.px_per_row;
    app.tree.px_per_unit = v.px_per_unit;
}

/// Fit the whole tree into the panel: the rows fill the height above the
/// footer, and the tree leaves exactly the room its labels, group labels and
/// cluster brackets need on the right.
fn fit_view(
    app: &mut InfinityApp,
    rect: Rect,
    l: &TreeLayout,
    marks: &pubtree::Marks,
    opts: &pubtree::PubOptions,
) {
    let guess = View {
        origin: rect.min,
        px_per_row: app.tree.px_per_row,
        px_per_unit: app.tree.px_per_unit,
    };
    let used = pubtree::Used {
        exact: true,
        estimate: true,
        groups: opts.groups && !marks.groups.is_empty(),
        selection: !marks.selection.is_empty(),
        wedges: true,
        tip_snps: opts.tip_snps,
    };
    let footer = pubtree::footer_height(
        l,
        marks,
        &guess.geometry(rect),
        &used,
        opts.legend,
        rect.width(),
    );
    let h = (rect.height() - 2.0 * MARGIN - footer - 4.0).max(10.0);
    let ppr = (h / l.rows.max(1.0) as f32).clamp(MIN_PX_PER_ROW, MAX_FIT_ROW);
    let room = pubtree::label_room(l, marks, opts, ppr);
    let w = (rect.width() - 2.0 * MARGIN - room).max(rect.width() * 0.35);
    app.tree.px_per_row = ppr;
    app.tree.px_per_unit = (w / l.x_max.max(1e-6) as f32).clamp(MIN_PX_PER_UNIT, MAX_PX_PER_UNIT);
    app.tree.offset = Vec2::ZERO;
}

fn genome_tooltip(model: &UiModel, id: u32) -> Vec<String> {
    let Some(g) = model.genome(id) else {
        return vec!["unknown genome".to_string()];
    };
    let mut out = vec![
        g.name.clone(),
        format!(
            "length {} | contigs {} | N50 {} | GC {}",
            fmt::bp(g.length),
            fmt::count(u64::from(g.contigs)),
            fmt::bp(g.n50),
            fmt::pct(Some(g.gc))
        ),
        format!(
            "cluster {} | {} | {}",
            g.cluster
                .map(|c| c.to_string())
                .unwrap_or_else(|| "-".to_string()),
            if g.is_pivot { "pivot" } else { "member" },
            fmt::focus_label(g.focus, g.approximate)
        ),
    ];
    match &g.numbers {
        Some(n) => {
            let m = fmt::approx_mark(g.approximate);
            out.push(format!(
                "against the pivot: {m}SNPs {} | indels {} | diff pieces {}",
                fmt::count(n.snps),
                fmt::count(n.indels),
                fmt::count(n.diff_pieces)
            ));
            out.push(format!(
                "uncovered {} | novel {}",
                fmt::bp(n.uncovered_bp),
                fmt::bp(n.novel_bp)
            ));
        }
        None => out.push("no numbers yet".to_string()),
    }
    out.push("click to select · Ctrl-click to add or remove · Shift-click for a range".to_string());
    out
}

fn cluster_tooltip(model: &UiModel, id: u32, opened: bool) -> Vec<String> {
    let Some(c) = model.cluster(id) else {
        return vec!["unknown cluster".to_string()];
    };
    vec![
        format!("cluster {} — pivot {}", c.id, c.pivot_name),
        format!(
            "{} genomes · {}",
            c.size(),
            fmt::focus_label(c.focus, c.approximate)
        ),
        format!(
            "branch lengths in {} (from {})",
            fmt::marked_unit(c.focus),
            fmt::source_unit(c.focus)
        ),
        if opened {
            "click to close it".to_string()
        } else {
            "click to open it (and to refine it next)".to_string()
        },
    ]
}

fn group_tooltip(app: &InfinityApp, i: usize) -> Vec<String> {
    let set = &app.groups.set;
    let Some(g) = set.groups.get(i) else {
        return Vec::new();
    };
    vec![
        format!(
            "G{} — {} genomes in cluster {}",
            g.number,
            g.members.len(),
            g.cluster
        ),
        format!(
            "linked by chains of ≤{} SNPs; the widest pair is {} SNPs apart",
            set.threshold,
            fmt::snps(g.diameter)
        ),
        "click to select them · Ctrl-click to add them".to_string(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::selection::Selection;

    fn view() -> View {
        View {
            origin: Pos2::new(100.0, 50.0),
            px_per_row: 16.0,
            px_per_unit: 400.0,
        }
    }

    #[test]
    fn screen_and_layout_are_inverses() {
        let v = view();
        for p in [(0.0, 0.0), (1.0, 3.0), (0.25, 12.5)] {
            let s = v.screen(p);
            let back = v.layout(s);
            assert!((back.0 - p.0).abs() < 1e-5, "{p:?} -> {s:?} -> {back:?}");
            assert!((back.1 - p.1).abs() < 1e-5);
        }
        assert_eq!(v.screen((0.0, 0.0)), Pos2::new(100.0, 50.0));
        assert_eq!(v.screen((1.0, 2.0)), Pos2::new(500.0, 82.0));
    }

    #[test]
    fn zooming_keeps_the_point_under_the_pointer_fixed() {
        let v = view();
        let pointer = Pos2::new(233.0, 91.0);
        let before = v.layout(pointer);
        for f in [0.5f32, 1.25, 3.0, 0.1] {
            let z = v.zoomed_about(pointer, f);
            let after = z.layout(pointer);
            assert!(
                (after.0 - before.0).abs() < 1e-3 && (after.1 - before.1).abs() < 1e-3,
                "factor {f} moved {before:?} to {after:?}"
            );
        }
    }

    #[test]
    fn zoom_is_clamped_so_the_tree_cannot_be_lost() {
        let v = view();
        let pointer = Pos2::new(10.0, 10.0);
        let tiny = (0..80).fold(v, |acc, _| acc.zoomed_about(pointer, 0.2));
        assert!(tiny.px_per_row >= MIN_PX_PER_ROW);
        assert!(tiny.px_per_unit >= MIN_PX_PER_UNIT);
        let huge = (0..80).fold(v, |acc, _| acc.zoomed_about(pointer, 5.0));
        assert!(huge.px_per_row <= MAX_PX_PER_ROW);
        assert!(huge.px_per_unit <= MAX_PX_PER_UNIT);
        assert!(huge.px_per_row.is_finite() && huge.origin.x.is_finite());
    }

    #[test]
    fn the_row_order_expands_collapsed_clusters_into_their_members() {
        let mut m = UiModel::default();
        for i in 0..5u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2]), (1, 3, &[3, 4])]),
        );
        m.apply(
            0.0,
            infinity_core::events::Event::BackboneTree {
                newick: "(g0:0.001,g3:0.002);".into(),
                pivots: vec!["g0".into(), "g3".into()],
            },
        );
        m.apply(
            0.0,
            infinity_core::events::Event::ClusterFocus(infinity_core::events::ClusterFocusInfo {
                cluster: 1,
                focus: 2,
                approximate: false,
                newick: "(g4:1,g3:2);".into(),
                members: vec![],
            }),
        );
        let open: std::collections::HashSet<u32> = [1u32].into_iter().collect();
        let l = crate::layout::layout(&m.display_tree(&open).unwrap());
        let order = row_order(&l, &m);
        assert_eq!(order.len(), 5, "{order:?}");
        // The collapsed cluster 0 contributes its members as one block, and
        // the opened cluster 1 its leaves in row order.
        let pos0 = order.iter().position(|&g| g == 0).unwrap();
        assert_eq!(&order[pos0..pos0 + 3], &[0, 1, 2]);
        assert!(order.contains(&3) && order.contains(&4));
        let mut sel = Selection::default();
        sel.extend([1, 2, 4]);
        assert_eq!(sel.len(), 3);
    }

    #[test]
    fn a_view_round_trips_through_the_pan_offset() {
        let rect = Rect::from_min_size(Pos2::new(20.0, 30.0), Vec2::new(400.0, 300.0));
        let mut app = InfinityApp::default();
        app.tree.offset = Vec2::new(13.0, -7.0);
        let v = current_view(&app, rect);
        // Storing it back must give exactly the same offset, also after a zoom.
        apply_view(&mut app, rect, v);
        assert_eq!(app.tree.offset, Vec2::new(13.0, -7.0));
        let z = v.zoomed_about(Pos2::new(200.0, 100.0), 1.7);
        apply_view(&mut app, rect, z);
        assert_eq!(current_view(&app, rect), z);
    }

    #[test]
    fn cluster_states_read_as_words() {
        assert_eq!(state_word(2, false), "exact");
        assert_eq!(state_word(2, true), "≈ estimate");
        assert_eq!(state_word(0, true), "≈ estimate");
    }
}
