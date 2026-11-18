//! The tree panel (section 8): a rectangular tree drawn with `egui::Painter`,
//! with pan and zoom, a scale bar, collapsed clusters that open on click, focus
//! styling, a legend, hover numbers and a search box.
//!
//! egui has no tree widget, so the geometry comes from [`crate::layout`] and this
//! file only turns layout coordinates into pixels and pixels back into hits. The
//! transform is a small pure type ([`View`]) so that the mapping can be tested
//! without a window.
//!
//! The canvas is drawn in one read-only pass over the model; clicks are collected
//! and applied afterwards, so drawing never has to fight the borrow checker.

use eframe::egui;
use egui::{Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};

use crate::app::{InfinityApp, TreeMode, TreeView};
use crate::layout::{self, NodeKind, TreeLayout};
use crate::pairwise::PairSource;
use crate::state::UiModel;
use crate::{colors, fmt, theme};

/// Fraction of the panel width kept free for leaf labels when fitting.
const LABEL_FRACTION: f32 = 0.32;
const MARGIN: f32 = 10.0;
const MIN_PX_PER_ROW: f32 = 0.35;
const MAX_PX_PER_ROW: f32 = 240.0;
const MIN_PX_PER_UNIT: f32 = 8.0;
const MAX_PX_PER_UNIT: f32 = 40_000.0;
/// Below this row height labels would overlap, so they are dropped.
const MIN_ROW_FOR_LABELS: f32 = 7.0;
/// Below this row height the "≈" markers on Focus 1 branches are dropped, so a
/// zoomed-out tree does not turn into a field of glyphs. The legend still says
/// what Focus 1 means.
const MIN_ROW_FOR_APPROX: f32 = 12.0;

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
}

/// What the pointer is over.
#[derive(Debug, Clone, PartialEq)]
enum Hover {
    Genome(u32),
    Cluster(u32),
}

pub fn show(ui: &mut egui::Ui, app: &mut InfinityApp) {
    header(ui, app);
    canvas(ui, app);
}

// ---------------------------------------------------------------------------
// Header
// ---------------------------------------------------------------------------

fn header(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let mut fit_now = false;
    let mut zoom_by: Option<f32> = None;
    let mut picked: Option<u32> = None;

    // Wrapped rather than one long row: at the smallest window the controls
    // flow onto a second line instead of running under the side panel.
    theme::panel_header(ui, "Tree", |ui| {
        let sel = app.selected_cluster;
        let label = match sel.and_then(|c| app.model.cluster(c)) {
            Some(c) => format!(
                "c{} · {} genomes · {}{}",
                c.id,
                c.size(),
                fmt::focus_word(c.focus),
                if c.approximate { " ≈" } else { "" }
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
                                "cluster {} — {} genomes, {}{}",
                                c.id,
                                c.size(),
                                fmt::focus_word(c.focus),
                                if c.approximate { " ≈" } else { "" }
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
            ui.checkbox(&mut app.tree.show_legend, "Legend");
            ui.checkbox(&mut app.tree.branch_labels, "SNP counts on branches");
            ui.checkbox(&mut app.tree.aligned_tips, "Aligned tip labels");
            ui.checkbox(&mut app.tree.color_strip, "Cluster colour strip");
        });
    });
    // The pairwise readout gets a line of its own: it is the longest thing in
    // the header and the one that changes most.
    ui.horizontal_wrapped(|ui| pair_readout(ui, app));

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
                egui::RichText::new(format!(
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
        None => {
            theme::hint(ui, app.pair_hint());
        }
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

/// Everything the drawing pass needs, gathered once so that the draw functions
/// stay small.
struct Draw<'a> {
    painter: egui::Painter,
    view: View,
    rect: Rect,
    dark: bool,
    weak: Color32,
    strong: Color32,
    text: Color32,
    model: &'a UiModel,
    selection: &'a [u32],
    search_hit: Option<u32>,
    /// Node positions in layout coordinates, parallel to `TreeLayout::nodes`,
    /// mid-animation when the tree is moving.
    positions: &'a [(f64, f64)],
    /// Branch labels, aligned tips, colour strip.
    opts: &'a TreeView,
    selected_cluster: Option<u32>,
}

impl Draw<'_> {
    /// Screen x of the aligned tip-label column: just right of the deepest tip.
    fn tip_column(&self, l: &TreeLayout) -> f32 {
        self.view.screen((l.x_max, 0.0)).x + 10.0
    }
}

fn canvas(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let rect = remaining_rect(ui);
    let resp = ui.allocate_rect(rect, Sense::click_and_drag());
    app.tree_rect = Some((rect, ui.ctx().pixels_per_point()));
    let dark = ui.visuals().dark_mode;
    let pal = theme::palette(ui);
    // The canvas is the darkest (or lightest) surface of the window, framed by
    // a hairline: the "stage" the tree is drawn on.
    let bg = pal.canvas;
    let border = Stroke::new(1.0, pal.border);
    let weak = pal.text_weak;
    let strong = pal.text_strong;
    let text_col = pal.text;
    let painter = ui.painter().with_clip_rect(rect);
    painter.rect_filled(rect, 6, bg);
    painter.rect_stroke(rect, 6, border, egui::StrokeKind::Inside);

    // Pointer state is read once, up front: `Response` is consumed by the
    // tooltip call at the end of the frame.
    let pointer = resp.hover_pos();
    let dragged_by = if resp.dragged() {
        Some(resp.drag_delta())
    } else {
        None
    };
    let clicked = resp.clicked();
    let additive = ui.ctx().input(|i| i.modifiers.command);
    let scroll = ui.ctx().input(|i| i.smooth_scroll_delta);

    // --- interaction first, so this frame draws the result ------------------
    let fit_to = app.tree_layout.as_ref().map(|l| (l.rows, l.x_max));
    if app.tree.needs_fit && !app.tree.user_moved {
        if let Some((rows, x_max)) = fit_to {
            fit_view(app, rect, rows, x_max);
        }
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

    if app.tree_layout.is_none() {
        hint(
            &painter,
            rect,
            if app.model.genome_count() == 0 {
                "No genomes yet.\nAdd a file or a folder above and press Run."
            } else if app.model.backbone.is_none() {
                "Scanning genomes …\nThe blurry tree appears as soon as clustering is done."
            } else {
                "No tree to draw."
            },
            dark,
        );
        return;
    }

    let view = current_view(app, rect);

    // --- read-only drawing pass --------------------------------------------
    let hits: Vec<(Hover, Rect)>;
    let mut hover: Option<Hover> = None;
    let mut tooltip: Vec<String> = Vec::new();
    {
        let l = app.tree_layout.as_ref().expect("checked above");
        let d = Draw {
            painter: painter.clone(),
            view,
            rect,
            dark,
            weak,
            strong,
            text: text_col,
            model: &app.model,
            selection: &app.selection,
            search_hit: app.tree.search_hit,
            positions: &app.positions,
            opts: &app.tree,
            selected_cluster: app.selected_cluster,
        };
        draw_edges(&d, l);
        if d.opts.branch_labels {
            draw_branch_labels(&d, l);
        }
        hits = draw_nodes(&d, l);
        if let Some(p) = pointer {
            hover = hits
                .iter()
                .rev()
                .find(|(_, r)| r.contains(p))
                .map(|(h, _)| h.clone());
        }
        draw_scale_bars(&d, l, pointer);
        if let Some(h) = &hover {
            tooltip = match h {
                Hover::Genome(id) => genome_tooltip(d.model, *id),
                Hover::Cluster(id) => cluster_tooltip(d.model, *id),
            };
        }
    }

    if app.tree.show_legend {
        draw_legend(&painter, rect, dark, bg, border);
    }

    // --- clicks, applied after the read-only pass ---------------------------
    let mut click: Option<(Hover, bool)> = None;
    if clicked {
        if let Some(p) = pointer {
            click = hits
                .iter()
                .rev()
                .find(|(_, r)| r.contains(p))
                .map(|(h, _)| (h.clone(), additive));
        }
    }
    if let Some((h, add)) = click {
        match h {
            Hover::Genome(id) => app.select_genome(id, add),
            Hover::Cluster(id) => app.toggle_cluster(id),
        }
    }
    ui.ctx().output_mut(|o| {
        o.cursor_icon = if hover.is_some() {
            egui::CursorIcon::PointingHand
        } else {
            egui::CursorIcon::Default
        };
    });
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

fn hint(painter: &egui::Painter, rect: Rect, msg: &str, dark: bool) {
    painter.text(
        rect.center(),
        Align2::CENTER_CENTER,
        msg,
        FontId::proportional(15.0),
        if dark {
            Color32::from_gray(140)
        } else {
            Color32::from_gray(110)
        },
    );
}

fn current_view(app: &InfinityApp, rect: Rect) -> View {
    View {
        origin: rect.min + Vec2::splat(MARGIN) + app.tree.offset,
        px_per_row: app.tree.px_per_row,
        px_per_unit: app.tree.px_per_unit,
    }
}

/// Store a view back as a pan offset from the panel corner, so that resizing the
/// window keeps the tree where it was.
fn apply_view(app: &mut InfinityApp, rect: Rect, v: View) {
    app.tree.offset = v.origin - rect.min - Vec2::splat(MARGIN);
    app.tree.px_per_row = v.px_per_row;
    app.tree.px_per_unit = v.px_per_unit;
}

/// Fit the whole tree into the panel: the rows fill the height, and the
/// horizontal extent leaves [`LABEL_FRACTION`] of the width free for the labels.
fn fit_view(app: &mut InfinityApp, rect: Rect, rows: f64, x_max: f64) {
    let h = (rect.height() - 2.0 * MARGIN).max(10.0);
    let w = ((rect.width() - 2.0 * MARGIN) * (1.0 - LABEL_FRACTION)).max(10.0);
    app.tree.px_per_row = (h / rows.max(1.0) as f32).clamp(MIN_PX_PER_ROW, MAX_PX_PER_ROW);
    app.tree.px_per_unit = (w / x_max.max(1e-6) as f32).clamp(MIN_PX_PER_UNIT, MAX_PX_PER_UNIT);
    app.tree.offset = Vec2::ZERO;
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

fn draw_edges(d: &Draw<'_>, l: &TreeLayout) {
    let positions = d.positions;
    for n in &l.nodes {
        let Some(&p) = n.parent.and_then(|p| positions.get(p)) else {
            continue;
        };
        let Some(&c) = positions.get(n.idx) else {
            continue;
        };
        let style = colors::focus_style(n.focus, d.dark);
        let child = d.view.screen(c);
        let parent = d.view.screen(p);
        // Cull branches that are entirely off screen: a thousand-genome tree
        // shows a few dozen rows at a time.
        if (child.y < d.rect.top() - 2.0 && parent.y < d.rect.top() - 2.0)
            || (child.y > d.rect.bottom() + 2.0 && parent.y > d.rect.bottom() + 2.0)
        {
            continue;
        }
        // A rectangular tree: the branch runs horizontally at the child's height.
        segment(&d.painter, Pos2::new(parent.x, child.y), child, &style);
    }
    for n in &l.nodes {
        if n.leaf_count < 2 {
            continue;
        }
        let Some(&c) = positions.get(n.idx) else {
            continue;
        };
        let top = d.view.screen((c.0, n.row_min));
        let bottom = d.view.screen((c.0, n.row_max));
        if bottom.y < d.rect.top() - 2.0 || top.y > d.rect.bottom() + 2.0 {
            continue;
        }
        let style = colors::focus_style(n.focus, d.dark);
        segment(&d.painter, top, bottom, &style);
    }
}

/// One branch: dashed at Focus 0, solid above it (section 8).
fn segment(painter: &egui::Painter, a: Pos2, b: Pos2, style: &colors::FocusStyle) {
    if !style.dashed {
        painter.line_segment([a, b], style.stroke());
        return;
    }
    let delta = b - a;
    let len = delta.length();
    if len < 0.5 {
        return;
    }
    let step = style.dash.max(1.0);
    let unit = delta / len;
    let mut t = 0.0f32;
    while t < len {
        let end = (t + step).min(len);
        painter.line_segment([a + unit * t, a + unit * end], style.stroke());
        t = end + step;
    }
}

fn draw_nodes(d: &Draw<'_>, l: &TreeLayout) -> Vec<(Hover, Rect)> {
    let mut hits: Vec<(Hover, Rect)> = Vec::new();
    let label_font = FontId::proportional((d.view.px_per_row * 0.72).clamp(8.5, 14.0));
    let small_font = FontId::proportional((d.view.px_per_row * 0.62).clamp(8.0, 12.0));
    let show_labels = d.view.px_per_row >= MIN_ROW_FOR_LABELS;
    let show_approx = d.view.px_per_row >= MIN_ROW_FOR_APPROX;
    let weak = d.weak;
    let strong = d.strong;
    let text_col = d.text;

    for n in &l.nodes {
        let Some(&p) = d.positions.get(n.idx) else {
            continue;
        };
        let s = d.view.screen(p);
        let y = s.y;
        if y < d.rect.top() - 20.0 || y > d.rect.bottom() + 20.0 {
            continue;
        }
        match &n.kind {
            NodeKind::Genome { id, name } => {
                if !show_labels {
                    d.painter
                        .circle_filled(s, 1.0, colors::focus_color(n.focus, d.dark));
                    if d.selection.contains(id) {
                        d.painter.circle_filled(s, 2.2, colors::SELECTED);
                    }
                    hits.push((
                        Hover::Genome(*id),
                        Rect::from_center_size(s, Vec2::splat(6.0)),
                    ));
                    continue;
                }
                let selected = d.selection.contains(id);
                let colour = d
                    .model
                    .genome(*id)
                    .and_then(|g| g.cluster)
                    .map(|c| colors::cluster_color(c, d.dark))
                    .unwrap_or(weak);
                if selected {
                    let r = Rect::from_min_size(
                        Pos2::new(d.rect.left(), y - d.view.px_per_row * 0.5),
                        Vec2::new(d.rect.width(), d.view.px_per_row),
                    );
                    d.painter
                        .rect_filled(r, 0, colors::wash(colors::SELECTED, 46));
                }
                d.painter.circle_filled(s, 2.0, colour);
                if Some(*id) == d.search_hit {
                    d.painter
                        .circle_stroke(s, 5.0, Stroke::new(1.4, colors::HOVERED));
                }
                // The marker travels with the name, so an approximate branch is
                // never readable as an exact one (section 0 rule 8).
                let text = if n.approximate {
                    format!("≈ {name}")
                } else {
                    name.clone()
                };
                // Aligned tips start in one column, joined to the tip by a
                // dotted leader; otherwise the label sits at the tip.
                let mut lx = if d.opts.aligned_tips {
                    let col = d.tip_column(l);
                    if col - s.x > 8.0 {
                        dotted(&d.painter, s.x + 4.0, col - 3.0, y, colors::wash(weak, 110));
                    }
                    col
                } else {
                    s.x + 5.0
                };
                if d.opts.color_strip {
                    let h = (d.view.px_per_row * 0.86).max(2.0);
                    let strip = Rect::from_min_size(Pos2::new(lx, y - h / 2.0), Vec2::new(7.0, h));
                    d.painter.rect_filled(strip, 1, colour);
                    lx += 11.0;
                }
                let tr = d.painter.text(
                    Pos2::new(lx, y),
                    Align2::LEFT_CENTER,
                    text,
                    label_font.clone(),
                    if selected { strong } else { text_col },
                );
                hits.push((
                    Hover::Genome(*id),
                    Rect::from_min_max(s - Vec2::splat(3.0), tr.max + Vec2::new(2.0, 3.0)),
                ));
            }
            NodeKind::ClusterCollapsed { id, members, .. } => {
                // The triangle geometry lives in `layout`, where it is tested.
                let pts = layout::triangle(p.0, p.1, layout::collapsed_half_rows(*members));
                let tri = [
                    d.view.screen(pts[0]),
                    d.view.screen(pts[1]),
                    d.view.screen(pts[2]),
                ];
                let base_x = tri[1].x;
                let colour = colors::cluster_color(*id, d.dark);
                d.painter.add(egui::Shape::convex_polygon(
                    tri.to_vec(),
                    colors::wash(colour, 110),
                    Stroke::new(1.2, colour),
                ));
                let (focus, approx) = d
                    .model
                    .cluster(*id)
                    .map(|c| (c.focus, c.approximate))
                    .unwrap_or((n.focus, n.approximate));
                let tri_rect = Rect::from_points(&tri[..]);
                if show_labels {
                    let label = format!(
                        "c{id} ({members} genomes) {}{}",
                        fmt::focus_word(focus),
                        if approx { " ≈" } else { "" }
                    );
                    let tr = d.painter.text(
                        Pos2::new(base_x + 5.0, y),
                        Align2::LEFT_CENTER,
                        label,
                        label_font.clone(),
                        colour,
                    );
                    hits.push((Hover::Cluster(*id), tri_rect.union(tr)));
                } else {
                    hits.push((Hover::Cluster(*id), tri_rect));
                }
            }
            NodeKind::ClusterRoot { id, members, focus } => {
                let colour = colors::cluster_color(*id, d.dark);
                d.painter.circle_filled(s, 3.0, colour);
                d.painter
                    .circle_stroke(s, 4.6, Stroke::new(1.0, colors::wash(colour, 150)));
                if show_labels {
                    let label = format!(
                        "c{id} · {members} genomes · {}{}",
                        fmt::focus_word(*focus),
                        if n.approximate { " ≈" } else { "" }
                    );
                    let tr = d.painter.text(
                        Pos2::new(s.x + 7.0, y - d.view.px_per_row * 0.75),
                        Align2::LEFT_CENTER,
                        label,
                        small_font.clone(),
                        colour,
                    );
                    hits.push((Hover::Cluster(*id), tr.expand(2.0)));
                } else {
                    hits.push((
                        Hover::Cluster(*id),
                        Rect::from_center_size(s, Vec2::splat(8.0)),
                    ));
                }
            }
            NodeKind::Internal | NodeKind::Root => {
                if show_approx && n.approximate && n.focus == fmt::FOCUS_BLOCKS && n.leaf_count > 1
                {
                    d.painter.text(
                        Pos2::new(s.x + 3.0, y - d.view.px_per_row * 0.45),
                        Align2::LEFT_BOTTOM,
                        "≈",
                        small_font.clone(),
                        colors::focus_color(fmt::FOCUS_BLOCKS, d.dark),
                    );
                }
            }
        }
    }
    hits
}

/// The scale bars of section 8 ("|-----| 5 SNPs"), in SNPs. One for the
/// opened cluster the pointer is over (or the selected one, or the first one),
/// and one for the backbone between clusters, whose SNPs are estimates and say
/// so. Each bar is drawn in its region's focus colour.
fn draw_scale_bars(d: &Draw<'_>, l: &TreeLayout, pointer: Option<Pos2>) {
    let mut bars: Vec<(&layout::Region, &str)> = Vec::new();
    let under_pointer = pointer
        .map(|p| d.view.layout(p).0)
        .and_then(|lx| l.region_at(lx))
        .filter(|r| r.cluster.is_some());
    let cluster_region = under_pointer
        .or_else(|| {
            l.regions
                .iter()
                .find(|r| r.cluster.is_some() && r.cluster == d.selected_cluster)
        })
        .or_else(|| l.regions.iter().find(|r| r.cluster.is_some()));
    if let Some(r) = cluster_region {
        bars.push((r, "in cluster"));
    }
    if let Some(bb) = l.regions.first().filter(|r| r.raw_max > 0.0) {
        bars.push((bb, "between clusters"));
    }
    let mut y = d.rect.bottom() - MARGIN - 10.0;
    for (region, what) in bars {
        let Some(bar) = layout::scale_bar(region, 0.18) else {
            continue;
        };
        let px = bar.layout_width as f32 * d.view.px_per_unit;
        if !px.is_finite() || px <= 4.0 {
            continue;
        }
        let x0 = d.rect.left() + MARGIN;
        let style = colors::focus_style(region.focus, d.dark);
        let stroke = Stroke::new(1.6, style.color);
        segment(&d.painter, Pos2::new(x0, y), Pos2::new(x0 + px, y), &style);
        d.painter
            .line_segment([Pos2::new(x0, y - 4.0), Pos2::new(x0, y + 4.0)], stroke);
        d.painter.line_segment(
            [Pos2::new(x0 + px, y - 4.0), Pos2::new(x0 + px, y + 4.0)],
            stroke,
        );
        let label = match region.cluster {
            Some(c) => format!("{}  ({what} c{c})", bar.label),
            None => format!("{}  ({what})", bar.label),
        };
        d.painter.text(
            Pos2::new(x0 + px + 6.0, y),
            Align2::LEFT_CENTER,
            label,
            FontId::proportional(11.0),
            stroke.color,
        );
        y -= 18.0;
    }
}

/// SNP counts written above the branches, where there is room for them. The
/// "≈" marker travels with every estimated count.
fn draw_branch_labels(d: &Draw<'_>, l: &TreeLayout) {
    if d.view.px_per_row < 11.0 {
        return;
    }
    let font = FontId::proportional((d.view.px_per_row * 0.55).clamp(8.5, 11.0));
    let char_w = font.size * 0.56;
    for n in &l.nodes {
        if n.len <= 0.0 || matches!(n.kind, NodeKind::Root) {
            continue;
        }
        let (Some(&c), Some(&p)) = (
            d.positions.get(n.idx),
            n.parent.and_then(|p| d.positions.get(p)),
        ) else {
            continue;
        };
        let child = d.view.screen(c);
        let parent = d.view.screen(p);
        if child.y < d.rect.top() || child.y > d.rect.bottom() {
            continue;
        }
        let text = format!(
            "{}{}",
            if n.approximate { "≈" } else { "" },
            fmt::snps(n.len)
        );
        let w = text.chars().count() as f32 * char_w;
        if child.x - parent.x < w + 6.0 {
            continue;
        }
        d.painter.text(
            Pos2::new((parent.x + child.x) / 2.0, child.y - 1.5),
            Align2::CENTER_BOTTOM,
            text,
            font.clone(),
            colors::wash(colors::focus_color(n.focus, d.dark), 210),
        );
    }
}

/// A dotted horizontal leader from a tip to its aligned label.
fn dotted(painter: &egui::Painter, x0: f32, x1: f32, y: f32, colour: Color32) {
    painter.extend(egui::Shape::dashed_line(
        &[Pos2::new(x0, y), Pos2::new(x1, y)],
        Stroke::new(0.8, colour),
        1.2,
        2.8,
    ));
}

fn draw_legend(painter: &egui::Painter, rect: Rect, dark: bool, bg: Color32, border: Stroke) {
    let box_rect = Rect::from_min_size(rect.left_top() + Vec2::splat(8.0), Vec2::new(214.0, 62.0));
    // Half way between the canvas and a panel, slightly see-through, so the
    // legend reads as floating above the tree rather than as part of it.
    let fill = theme::mix(bg, theme::palette_for(dark).panel, 0.6);
    painter.rect_filled(box_rect, 6, theme::with_alpha(fill, 225));
    painter.rect_stroke(box_rect, 6, border, egui::StrokeKind::Inside);
    let font = FontId::proportional(11.0);
    for (i, text) in [
        "Focus 0 — sketch only, ≈ SNPs",
        "Focus 1 — blocks compared, ≈ SNPs",
        "Focus 2 — exact SNPs",
    ]
    .iter()
    .enumerate()
    {
        let y = box_rect.top() + 14.0 + i as f32 * 17.0;
        let style = colors::focus_style(i as u8, dark);
        let a = Pos2::new(box_rect.left() + 10.0, y);
        let b = Pos2::new(box_rect.left() + 40.0, y);
        segment(painter, a, b, &style);
        painter.text(
            Pos2::new(box_rect.left() + 46.0, y),
            Align2::LEFT_CENTER,
            *text,
            font.clone(),
            style.color,
        );
    }
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
                "{m}diff pieces {} | SNPs {} | indels {}",
                fmt::count(n.diff_pieces),
                fmt::count(n.snps),
                fmt::count(n.indels)
            ));
            out.push(format!(
                "uncovered {} | novel {}",
                fmt::bp(n.uncovered_bp),
                fmt::bp(n.novel_bp)
            ));
        }
        None => out.push("no numbers yet".to_string()),
    }
    out.push("click to select · ctrl-click for a pair".to_string());
    out
}

fn cluster_tooltip(model: &UiModel, id: u32) -> Vec<String> {
    let Some(c) = model.cluster(id) else {
        return vec!["unknown cluster".to_string()];
    };
    vec![
        format!("cluster {} — pivot {}", c.id, c.pivot_name),
        format!("{} genomes", c.size()),
        fmt::focus_label(c.focus, c.approximate),
        format!(
            "branch lengths in {} (from {})",
            fmt::marked_unit(c.focus),
            fmt::source_unit(c.focus)
        ),
        "click to open or close, and to refine this cluster next".to_string(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn a_view_round_trips_through_the_pan_offset() {
        let rect = Rect::from_min_size(Pos2::new(20.0, 30.0), Vec2::new(400.0, 300.0));
        let tree = crate::app::TreeView {
            offset: Vec2::new(13.0, -7.0),
            ..Default::default()
        };
        let v = View {
            origin: rect.min + Vec2::splat(MARGIN) + tree.offset,
            px_per_row: tree.px_per_row,
            px_per_unit: tree.px_per_unit,
        };
        // Storing it back must give exactly the same offset.
        let offset = v.origin - rect.min - Vec2::splat(MARGIN);
        assert_eq!(offset, tree.offset);
    }
}
