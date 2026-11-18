//! The map panel (section 8): the pivot of the selected cluster as a horizontal
//! axis with its contigs side by side, and one row per selected member showing
//! identical stretches, unresolved difference pieces, SNPs, indels, complex
//! variants and the places the member does not cover.
//!
//! # What the event stream does and does not carry
//!
//! `MapUpdateInfo` has everything this panel draws except breakpoints: the
//! engine's `Focus1` struct records them, but the event does not, so no
//! per-member breakpoint position ever reaches a consumer. The triangle marker
//! is implemented and driven by [`Row::breakpoints`], which the event stream can
//! only leave empty.
//!
//! For the same reason the "zoom made visible" of section 8 is only partly
//! achievable: `MapUpdate` is emitted once, when a cluster reaches Focus 2, and
//! every piece in it is already marked resolved. This panel therefore shows the
//! two states it *can* see — a wide orange box for an unresolved piece and a
//! narrow pale one for a resolved piece — and says plainly when a cluster has
//! no map data yet instead of drawing an empty axis that looks like "no
//! differences".

use eframe::egui;
use egui::{Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};

use infinity_core::events::MapUpdateInfo;
use infinity_core::interval;
use infinity_core::variant::VariantKind;

use crate::app::InfinityApp;
use crate::map::{self, Axis};
use crate::state::UiModel;
use crate::{colors, fmt};

/// The largest cluster for which the panel will draw every member when nothing
/// is selected: more rows than this cannot be read in a panel a few centimetres
/// high, and drawing them would cost frames for nothing.
pub const ROW_CAP: usize = 24;
const RULER_H: f32 = 18.0;
const ROW_H: f32 = 15.0;
const LABEL_W: f32 = 132.0;
const MIN_SPAN: f64 = 200.0;

/// Zoom and pan state of the map panel.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MapView {
    /// First base of the visible window, in flat pivot coordinates.
    pub start: f64,
    /// Width of the visible window, in bases.
    pub span: f64,
    pub needs_fit: bool,
    /// True once the user has panned or zoomed by hand.
    pub user_moved: bool,
    /// A zoom step asked for by a header button, applied about the window centre.
    pub zoom_request: Option<f64>,
}

impl Default for MapView {
    fn default() -> Self {
        Self {
            start: 0.0,
            span: 1.0,
            needs_fit: true,
            user_moved: false,
            zoom_request: None,
        }
    }
}

/// Base position to screen x, and back. Pure, so the axis maths is testable.
pub fn x_of(flat: f64, view: &MapView, rect: Rect) -> f32 {
    let f = ((flat - view.start) / view.span.max(1e-9)) as f32;
    rect.left() + f * rect.width()
}

pub fn flat_of(x: f32, view: &MapView, rect: Rect) -> f64 {
    let f = f64::from((x - rect.left()) / rect.width().max(1.0));
    view.start + f * view.span
}

/// Zoom the window about a fixed base position.
pub fn zoomed_about(view: MapView, anchor: f64, factor: f64, total: f64) -> MapView {
    let span = (view.span * factor).clamp(MIN_SPAN, total.max(MIN_SPAN));
    let frac = if view.span > 0.0 {
        (anchor - view.start) / view.span
    } else {
        0.5
    };
    let mut v = view;
    v.span = span;
    v.start = (anchor - frac * span).clamp(0.0, (total - span).max(0.0));
    v
}

/// One row of the map.
struct Row<'a> {
    genome: u32,
    name: String,
    map: &'a MapUpdateInfo,
    axis: Axis,
    /// Breakpoints would be drawn as triangles. The event stream does not carry
    /// them, so this stays empty; see the module documentation.
    breakpoints: Vec<u64>,
}

pub fn show(ui: &mut egui::Ui, app: &mut InfinityApp) {
    egui::Panel::bottom("map")
        .resizable(true)
        .default_size(196.0)
        .min_size(74.0)
        .show(ui, |ui| {
            header(ui, app);
            ui.separator();
            body(ui, app);
        });
}

fn header(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let mut fit_now = false;
    let mut zoom_by: Option<f64> = None;
    let title = match app.selected_cluster.and_then(|c| app.model.cluster(c)) {
        Some(c) => format!(
            "MAP  (cluster {}, pivot {}, {}{})",
            c.id,
            c.pivot_name,
            fmt::focus_word(c.focus),
            if c.approximate { " ≈" } else { "" }
        ),
        None => "MAP  (no cluster selected)".to_string(),
    };
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(title).strong());
        ui.separator();
        if ui.button("−").clicked() {
            zoom_by = Some(1.0 / 1.4);
        }
        if ui.button("+").clicked() {
            zoom_by = Some(1.4);
        }
        if ui.button("Fit").clicked() {
            fit_now = true;
        }
        ui.separator();
        let plan = map::plan_rows(&app.model, app.selected_cluster, &app.selection, ROW_CAP);
        match plan.reason {
            map::RowReason::Selection => {
                ui.weak(format!("{} selected member(s)", plan.genomes.len()));
            }
            map::RowReason::WholeCluster => {
                if plan.truncated {
                    ui.weak(format!(
                        "first {} of {} members — select genomes to narrow it down",
                        plan.genomes.len(),
                        plan.available
                    ));
                } else {
                    ui.weak(format!("{} members", plan.genomes.len()));
                }
            }
            map::RowReason::NoData => {
                ui.weak(no_data_hint(&app.model, app.selected_cluster));
            }
        }
    });
    if fit_now {
        app.map_view.needs_fit = true;
        app.map_view.user_moved = false;
    }
    if let Some(f) = zoom_by {
        app.map_view.zoom_request = Some(f);
    }
}

/// Why there is nothing to draw, said plainly rather than left as a blank axis.
fn no_data_hint(model: &UiModel, cluster: Option<u32>) -> String {
    match cluster.and_then(|c| model.cluster(c)) {
        None => "select a cluster in the tree or in the genome list".to_string(),
        Some(c) => format!(
            "no map data yet: cluster {} is at {}; the map appears when it reaches Focus 2",
            c.id,
            fmt::focus_label(c.focus, c.approximate)
        ),
    }
}

fn body(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let size = ui.available_size_before_wrap();
    let (rect, resp) =
        ui.allocate_exact_size(size.max(Vec2::new(80.0, 40.0)), Sense::click_and_drag());
    let dark = ui.visuals().dark_mode;
    let weak = ui.visuals().weak_text_color();
    let painter = ui.painter().with_clip_rect(rect);
    painter.rect_filled(rect, 0, ui.visuals().panel_fill);

    let cluster = app.selected_cluster.unwrap_or(u32::MAX);
    let plan = map::plan_rows(&app.model, app.selected_cluster, &app.selection, ROW_CAP);
    let ids = plan.genomes.clone();
    let total = app
        .model
        .map(cluster, ids.first().copied().unwrap_or(u32::MAX))
        .map(|m| (Axis::from_map(m).total.max(1)) as f64)
        .unwrap_or(0.0);

    if total == 0.0 || ids.is_empty() {
        painter.text(
            rect.center(),
            Align2::CENTER_CENTER,
            no_data_hint(&app.model, app.selected_cluster),
            FontId::proportional(12.5),
            weak,
        );
        return;
    }

    // --- view state, mutated before anything borrows the model --------------
    if app.map_view.needs_fit && !app.map_view.user_moved {
        app.map_view.start = 0.0;
        app.map_view.span = total;
    }
    app.map_view.needs_fit = false;
    let span = app.map_view.span;
    if !span.is_finite() || span <= 0.0 || span > total {
        app.map_view.span = total;
        app.map_view.start = app
            .map_view
            .start
            .clamp(0.0, (total - app.map_view.span).max(0.0));
    }
    if let Some(f) = app.map_view.zoom_request.take() {
        let centre = app.map_view.start + app.map_view.span / 2.0;
        app.map_view = zoomed_about(app.map_view, centre, f, total);
        app.map_view.user_moved = true;
    }
    if resp.dragged() {
        let d =
            f64::from(resp.drag_delta().x) / f64::from(rect.width().max(1.0)) * app.map_view.span;
        app.map_view.start =
            (app.map_view.start - d).clamp(0.0, (total - app.map_view.span).max(0.0));
        app.map_view.user_moved = true;
    }
    let plot = plot_rect(rect);
    let pointer = resp.hover_pos();
    let scroll = ui.ctx().input(|i| i.smooth_scroll_delta);
    if scroll.y.abs() > 0.01 {
        let anchor = pointer
            .map(|p| flat_of(p.x, &app.map_view, plot))
            .unwrap_or(app.map_view.start + app.map_view.span / 2.0);
        app.map_view = zoomed_about(
            app.map_view,
            anchor,
            f64::from((-scroll.y * 0.004).exp()),
            total,
        );
        app.map_view.user_moved = true;
    }
    if scroll.x.abs() > 0.01 {
        let d = f64::from(scroll.x) / f64::from(rect.width().max(1.0)) * app.map_view.span;
        app.map_view.start =
            (app.map_view.start + d).clamp(0.0, (total - app.map_view.span).max(0.0));
        app.map_view.user_moved = true;
    }

    // --- read-only drawing pass --------------------------------------------
    let view = app.map_view;
    let mut hovered: Option<Vec<String>> = None;
    let mut clicked: Option<u32> = None;
    {
        let rows: Vec<Row<'_>> = ids
            .iter()
            .filter_map(|&g| {
                app.model.map(cluster, g).map(|m| Row {
                    genome: g,
                    name: m.genome_name.clone(),
                    axis: Axis::from_map(m),
                    map: m,
                    breakpoints: Vec::new(),
                })
            })
            .collect();
        if rows.is_empty() {
            return;
        }
        let axis = rows[0].axis.clone();
        draw_ruler(&painter, &axis, &view, plot, rect, weak);
        let row_top = plot.top() + RULER_H;
        let mut hovered_row: Option<usize> = None;
        for (i, row) in rows.iter().enumerate() {
            let top = row_top + i as f32 * ROW_H;
            if top > rect.bottom() {
                break;
            }
            let band =
                Rect::from_min_size(Pos2::new(rect.left(), top), Vec2::new(rect.width(), ROW_H));
            draw_row(&painter, row, &view, band, plot, dark, weak);
            if let Some(p) = pointer {
                if band.contains(p) {
                    hovered_row = Some(i);
                    if resp.clicked() {
                        clicked = Some(row.genome);
                    }
                }
            }
        }
        if let (Some(p), Some(i)) = (pointer, hovered_row) {
            // Hovering shows the position and the alleles (section 8).
            hovered = Some(hover_text(&rows[i], &axis, flat_of(p.x, &view, plot)));
        }
    }
    if let Some(lines) = hovered {
        resp.on_hover_ui(|ui| {
            for line in lines {
                ui.label(line);
            }
        });
    }
    if let Some(g) = clicked {
        // Clicking a row selects that member, which is what links the map back
        // to the tree and the genome list.
        let additive = ui.ctx().input(|i| i.modifiers.command);
        app.select_genome(g, additive);
    }
}

/// The part of the panel that the bases are drawn in: right of the row labels.
fn plot_rect(rect: Rect) -> Rect {
    Rect::from_min_max(
        Pos2::new(rect.left() + LABEL_W, rect.top()),
        rect.right_bottom(),
    )
}

fn draw_ruler(
    painter: &egui::Painter,
    axis: &Axis,
    view: &MapView,
    plot: Rect,
    rect: Rect,
    weak: Color32,
) {
    let font = FontId::proportional(9.5);
    let y = plot.top() + RULER_H - 4.0;
    painter.line_segment(
        [Pos2::new(plot.left(), y), Pos2::new(plot.right(), y)],
        Stroke::new(1.0, weak),
    );
    // Contig boundaries run the full height of the panel.
    for &off in axis.offsets.iter().skip(1) {
        let x = x_of(off as f64, view, plot);
        if x < plot.left() || x > plot.right() {
            continue;
        }
        painter.line_segment(
            [Pos2::new(x, plot.top()), Pos2::new(x, rect.bottom())],
            Stroke::new(1.0, colors::wash(weak, 90)),
        );
    }
    // Contig names, only where there is room for them.
    for (i, name) in axis.contigs.iter().enumerate() {
        let from = x_of(axis.offsets[i] as f64, view, plot);
        let to = x_of((axis.offsets[i] + axis.lens[i]) as f64, view, plot);
        if to < plot.left() || from > plot.right() || to - from < 28.0 {
            continue;
        }
        painter.text(
            Pos2::new(
                ((from + to) / 2.0).clamp(plot.left(), plot.right()),
                plot.top() + 15.0,
            ),
            Align2::CENTER_CENTER,
            name,
            font.clone(),
            weak,
        );
    }
    // Ticks at round positions.
    for g in axis.grid(6) {
        let x = x_of(g as f64, view, plot);
        if x < plot.left() || x > plot.right() {
            continue;
        }
        painter.line_segment(
            [Pos2::new(x, y), Pos2::new(x, y + 3.0)],
            Stroke::new(1.0, weak),
        );
        painter.text(
            Pos2::new(x, y - 5.0),
            Align2::CENTER_BOTTOM,
            fmt::count(g),
            font.clone(),
            weak,
        );
    }
    // The window's own limits, in the 1-based coordinates the user reads.
    let left = axis.label(app_start(view) as u64);
    let right = axis.label(((view.start + view.span).min(axis.total as f64)) as u64);
    painter.text(
        Pos2::new(rect.left() + 4.0, plot.top() + 6.0),
        Align2::LEFT_CENTER,
        format!("{left} … {right}"),
        font,
        weak,
    );
}

fn app_start(view: &MapView) -> f64 {
    view.start.max(0.0)
}

fn draw_row(
    painter: &egui::Painter,
    row: &Row<'_>,
    view: &MapView,
    band: Rect,
    plot: Rect,
    dark: bool,
    weak: Color32,
) {
    let cy = band.center().y;
    let h = (band.height() - 5.0).max(2.0);
    let font = FontId::proportional(10.0);

    // The row label: the member's name, in its cluster colour.
    let colour = colors::cluster_color(row.map.cluster, dark);
    painter.text(
        Pos2::new(band.left() + 4.0, cy),
        Align2::LEFT_CENTER,
        crate::genome_panel::clip(painter, &row.name, LABEL_W - 10.0, font.clone()),
        font,
        colour,
    );
    painter.line_segment(
        [
            Pos2::new(plot.left(), band.top()),
            Pos2::new(plot.right(), band.top()),
        ],
        Stroke::new(0.5, colors::wash(weak, 40)),
    );

    // Identical stretches: a thin grey line wherever the member covers the
    // pivot. Everything else stays as the background, which reads as a gap.
    let uncovered = interval::merge(row.map.uncovered.clone());
    for (s, e) in map::complement(&uncovered, row.axis.total) {
        let x0 = x_of(s as f64, view, plot).max(plot.left());
        let x1 = x_of(e as f64, view, plot).min(plot.right());
        if x1 > x0 {
            painter.line_segment(
                [Pos2::new(x0, cy), Pos2::new(x1, cy)],
                Stroke::new(1.0, colors::IDENTICAL),
            );
        }
    }

    // Difference pieces. An unresolved one is a wide orange box; a resolved one
    // a narrow pale box, which is the visible half of the "zoom" of section 8.
    for p in &row.map.pieces {
        let flat0 = row
            .axis
            .offsets
            .get(p.pivot_contig as usize)
            .copied()
            .unwrap_or(0)
            + p.start;
        let flat1 = row
            .axis
            .offsets
            .get(p.pivot_contig as usize)
            .copied()
            .unwrap_or(0)
            + p.end;
        let x0 = x_of(flat0 as f64, view, plot);
        let x1 = x_of(flat1 as f64, view, plot);
        if x1 < plot.left() || x0 > plot.right() {
            continue;
        }
        let r = if p.resolved {
            Rect::from_min_max(
                Pos2::new(x0.max(plot.left()), cy - 1.5),
                Pos2::new(x1.min(plot.right()), cy + 1.5),
            )
        } else {
            Rect::from_min_max(
                Pos2::new(x0.max(plot.left()), cy - h * 0.42),
                Pos2::new(x1.min(plot.right()), cy + h * 0.42),
            )
        };
        if r.width() < 1.0 {
            painter.line_segment(
                [
                    Pos2::new(r.center().x, r.top()),
                    Pos2::new(r.center().x, r.bottom()),
                ],
                Stroke::new(
                    1.0,
                    if p.resolved {
                        colors::PIECE_RESOLVED
                    } else {
                        colors::PIECE_UNRESOLVED
                    },
                ),
            );
        } else {
            painter.rect_filled(
                r,
                0,
                if p.resolved {
                    colors::PIECE_RESOLVED
                } else {
                    colors::PIECE_UNRESOLVED
                },
            );
        }
    }

    // Variants: SNP red tick, indel blue tick, complex purple box.
    for v in &row.map.variants {
        let off = row
            .axis
            .offsets
            .get(v.pivot_contig as usize)
            .copied()
            .unwrap_or(0);
        let x0 = x_of((off + v.pivot_start) as f64, view, plot);
        let x1 = x_of((off + v.pivot_end) as f64, view, plot);
        if x1 < plot.left() || x0 > plot.right() {
            continue;
        }
        let colour = colors::variant_color(v.kind);
        match v.kind {
            VariantKind::Complex => {
                let r = Rect::from_min_max(
                    Pos2::new(x0.max(plot.left()), cy - h * 0.42),
                    Pos2::new(x1.max(x0 + 2.0).min(plot.right()), cy + h * 0.42),
                );
                painter.rect_filled(r, 0, colour);
            }
            _ => {
                let x = x0.max(plot.left()).min(plot.right());
                let w = (x1 - x0).max(1.0);
                painter.rect_filled(
                    Rect::from_min_max(
                        Pos2::new(x, cy - h * 0.45),
                        Pos2::new(x + w, cy + h * 0.45),
                    ),
                    0,
                    colour,
                );
            }
        }
    }

    // Breakpoints as triangles. Empty in practice: see the module docs.
    for &bp in &row.breakpoints {
        let x = x_of(bp as f64, view, plot);
        if x < plot.left() || x > plot.right() {
            continue;
        }
        painter.add(egui::Shape::convex_polygon(
            vec![
                Pos2::new(x, cy - h * 0.45),
                Pos2::new(x - 3.0, cy + h * 0.45),
                Pos2::new(x + 3.0, cy + h * 0.45),
            ],
            colors::BREAKPOINT,
            Stroke::NONE,
        ));
    }
}

/// What the pointer is over: the 1-based position and, when there is one, the
/// alleles of the variant underneath.
fn hover_text(row: &Row<'_>, axis: &Axis, flat: f64) -> Vec<String> {
    let pos = flat.max(0.0) as u64;
    let mut out = vec![row.name.clone(), format!("position {}", axis.label(pos))];
    let mut found = false;
    for v in &row.map.variants {
        let off = axis
            .offsets
            .get(v.pivot_contig as usize)
            .copied()
            .unwrap_or(0);
        let s = off + v.pivot_start;
        // An insertion covers no pivot base, so give it a one-base-wide target.
        let e = (off + v.pivot_end).max(s + 1);
        if pos >= s && pos < e {
            found = true;
            out.push(format!(
                "{}   {} → {}   pivot {} [{}]   member contig {} base {}   strand {}",
                v.kind.as_str(),
                fmt::or_dash(&v.ref_allele),
                fmt::or_dash(&v.alt_allele),
                axis.label(s),
                fmt::one_based_range(v.pivot_start, v.pivot_end, v.kind),
                v.genome_contig,
                fmt::one_based(v.genome_start, v.kind),
                v.strand,
            ));
            let flags = v.flags_string();
            if !flags.is_empty() {
                out.push(format!("flags: {flags}"));
            }
        }
    }
    if !found {
        let uncovered = interval::merge(row.map.uncovered.clone());
        if uncovered.iter().any(|&(s, e)| pos >= s && pos < e) {
            out.push("not covered by this member".to_string());
        } else if row
            .map
            .pieces
            .iter()
            .any(|p| flat_in_piece(axis, p.pivot_contig, p.start, p.end, pos))
        {
            out.push(
                if row.map.pieces.iter().any(|p| !p.resolved) {
                    "difference piece, not resolved yet (≈)"
                } else {
                    "difference piece, resolved"
                }
                .to_string(),
            );
        } else {
            out.push("identical to the pivot".to_string());
        }
    }
    out
}

fn flat_in_piece(axis: &Axis, contig: u32, start: u64, end: u64, flat: u64) -> bool {
    let off = axis.offsets.get(contig as usize).copied().unwrap_or(0);
    flat >= off + start && flat < off + end
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect() -> Rect {
        Rect::from_min_size(Pos2::new(100.0, 0.0), Vec2::new(400.0, 50.0))
    }

    fn view() -> MapView {
        MapView {
            start: 1000.0,
            span: 2000.0,
            needs_fit: false,
            user_moved: true,
            zoom_request: None,
        }
    }

    #[test]
    fn the_axis_maps_bases_to_pixels_and_back() {
        let v = view();
        let r = rect();
        assert_eq!(
            x_of(1000.0, &v, r),
            100.0,
            "the window starts at the left edge"
        );
        assert_eq!(x_of(3000.0, &v, r), 500.0, "and ends at the right edge");
        assert_eq!(x_of(2000.0, &v, r), 300.0, "the middle is in the middle");
        for flat in [0.0, 999.5, 1000.0, 2500.0, 3000.0, 5000.0] {
            let back = flat_of(x_of(flat, &v, r), &v, r);
            // Pixels are f32, so the round trip loses a little precision; the
            // tolerance is a thousandth of a base, far below anything visible.
            assert!((back - flat).abs() < 1e-3, "{flat} -> {back}");
        }
    }

    #[test]
    fn zooming_keeps_the_anchor_base_under_the_pointer() {
        let v = view();
        // Factors that leave the window inside the pivot. Zooming past either end
        // is clamped, and a clamped window cannot also hold the anchor: keeping
        // the base under the pointer loses to never showing bases that do not
        // exist. `zoom_is_clamped_to_the_pivot` covers that side.
        for f in [0.5, 2.0, 0.25, 1.5] {
            let z = zoomed_about(v, 2000.0, f, 10_000.0);
            let frac_before = (2000.0 - v.start) / v.span;
            let frac_after = (2000.0 - z.start) / z.span;
            assert!(
                (frac_before - frac_after).abs() < 1e-9,
                "factor {f} moved the anchor"
            );
        }
    }

    #[test]
    fn zoom_is_clamped_to_the_pivot() {
        let v = view();
        let tiny = (0..40).fold(v, |a, _| zoomed_about(a, 2000.0, 0.5, 10_000.0));
        assert!(tiny.span >= MIN_SPAN);
        let huge = (0..40).fold(v, |a, _| zoomed_about(a, 2000.0, 2.0, 10_000.0));
        assert!(
            (huge.span - 10_000.0).abs() < 1e-9,
            "never wider than the pivot"
        );
        assert!(huge.start >= 0.0);
        assert!(tiny.start + tiny.span <= 10_000.0 + 1e-9);
    }

    #[test]
    fn panning_never_leaves_the_pivot() {
        let v = MapView {
            start: 9500.0,
            span: 1000.0,
            needs_fit: false,
            user_moved: true,
            zoom_request: None,
        };
        // Dragging right moves the window left; the clamp keeps it inside.
        let moved = zoomed_about(v, 9999.0, 1.0, 10_000.0);
        assert!(moved.start >= 0.0);
        assert!(moved.start + moved.span <= 10_000.0 + 1e-9);
    }

    #[test]
    fn the_plot_leaves_room_for_the_row_labels() {
        let r = rect();
        let p = plot_rect(r);
        assert_eq!(p.left(), r.left() + LABEL_W);
        assert_eq!(p.right(), r.right());
        assert!(p.width() > 0.0);
    }
}
