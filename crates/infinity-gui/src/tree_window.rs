//! The tree window: a second native window with the whole tree as a
//! publication figure (white page, one line width, SNPs at the tips, SNP
//! groups shaded, scale bar and legend), whatever the theme of the main
//! window.
//!
//! It is an egui *deferred* viewport: it repaints on its own, so a big tree in
//! it never slows the main window, and closing it leaves the main window as it
//! was. The two talk through [`Shared`]: the main window puts a [`Figure`]
//! there whenever the model, the selection or the SNP groups change, and picks
//! up what the tree window collected (clicks on leaves select, with the same
//! modifiers as everywhere else; the SNP-group threshold can be changed here
//! too).
//!
//! Saving: SVG is written from the same primitives the window paints
//! ([`crate::pubtree`] and [`crate::svg::from_scene`]); PNG is a screenshot of
//! this window's canvas.

use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use eframe::egui;
use egui::{Color32, Pos2, Rect, Sense, Vec2};

use crate::app::InfinityApp;
use crate::pubtree::{self, Figure, Geometry, Hit, Ink, PubOptions};
use crate::selection::ClickMods;

/// Row height and tree width of the saved SVG.
const SVG_ROW_PX: f32 = 16.0;
const SVG_TREE_W: f32 = 640.0;
const MARGIN: f32 = 16.0;

pub fn viewport_id() -> egui::ViewportId {
    egui::ViewportId::from_hash_of("infinity-tree-window")
}

/// Pan and zoom of the tree window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PubView {
    pub px_per_row: f32,
    pub px_per_unit: f32,
    pub offset: Vec2,
    pub needs_fit: bool,
}

impl Default for PubView {
    fn default() -> Self {
        Self {
            px_per_row: 14.0,
            px_per_unit: 500.0,
            offset: Vec2::ZERO,
            needs_fit: true,
        }
    }
}

/// What the main window and the tree window share.
#[derive(Debug, Default)]
pub struct Shared {
    /// False once the user closed the tree window.
    pub open: bool,
    pub figure: Option<Figure>,
    /// What the figure was built from, to rebuild it only on a change.
    pub built_from: Option<(u64, u64, u64, u64)>,
    pub view: PubView,
    pub opts: PubOptions,
    /// Every cluster opened (the full tree), or the clusters opened in the
    /// main window.
    pub all_open: bool,
    /// Clicks on leaves, for the main window to apply.
    pub clicks: Vec<(u32, ClickMods)>,
    /// Clicks on clusters (brackets, wedges): make it the map's cluster.
    pub cluster_clicks: Vec<u32>,
    /// Clicks on a SNP group's label: select its genomes.
    pub group_clicks: Vec<(usize, bool)>,
    /// The SNP-group switch and threshold as the main window has them; a
    /// change made here is sent back through `groups_request`.
    pub groups: (bool, u32),
    pub groups_request: Option<(bool, u32)>,
    /// Where pictures are saved: the main window's export folder.
    pub export_dir: String,
    pub status: String,
    /// The canvas on screen and the pixels per point, for cropping a PNG.
    pub canvas: Option<(Rect, f32)>,
    /// A PNG requested by a button (or the screenshot harness): the path,
    /// and whether to crop to the canvas (false keeps the whole window).
    pub png_request: Option<(PathBuf, bool)>,
    pub png_result: Arc<Mutex<Option<Result<PathBuf, String>>>>,
}

pub type SharedRef = Arc<Mutex<Shared>>;

pub fn new_shared() -> SharedRef {
    Arc::new(Mutex::new(Shared {
        all_open: true,
        groups: (true, 5),
        ..Default::default()
    }))
}

fn lock(s: &SharedRef) -> std::sync::MutexGuard<'_, Shared> {
    s.lock().unwrap_or_else(|e| e.into_inner())
}

/// A cheap fingerprint of a list of ids.
fn signature(ids: impl IntoIterator<Item = u32>) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for id in ids {
        h ^= u64::from(id) + 1;
        h = h.wrapping_mul(0x100_0000_01b3).rotate_left(7);
    }
    h
}

/// Build the figure the tree window draws from the main window's model.
pub fn build_figure(app: &InfinityApp, all_open: bool) -> Option<Figure> {
    let open: HashSet<u32> = if all_open {
        app.model.clusters.iter().map(|c| c.id).collect()
    } else {
        app.expanded.clone()
    };
    let display = app.model.display_tree(&open)?;
    let layout = crate::layout::layout(&display);
    let marks = app.marks(&layout, false);
    Some(Figure { layout, marks })
}

/// Open the tree window (the header button). Calling it again while it is
/// open brings it to the front.
pub fn open(app: &mut InfinityApp, ctx: &egui::Context) {
    {
        let mut s = lock(&app.tree_window);
        if !s.open {
            s.open = true;
            s.view.needs_fit = true;
            s.built_from = None;
        }
    }
    if app.show_tree_window {
        ctx.send_viewport_cmd_to(viewport_id(), egui::ViewportCommand::Focus);
    }
    app.show_tree_window = true;
}

/// The main window's half: keep the figure current, apply what the tree
/// window collected, and keep the viewport alive while it is open.
pub fn show(ctx: &egui::Context, app: &mut InfinityApp) {
    if !app.show_tree_window {
        return;
    }
    let shared = Arc::clone(&app.tree_window);
    let (clicks, cluster_clicks, group_clicks, request, open, all_open) = {
        let mut s = lock(&shared);
        (
            std::mem::take(&mut s.clicks),
            std::mem::take(&mut s.cluster_clicks),
            std::mem::take(&mut s.group_clicks),
            s.groups_request.take(),
            s.open,
            s.all_open,
        )
    };
    if !open {
        app.show_tree_window = false;
        return;
    }
    let order = {
        let s = lock(&shared);
        s.figure
            .as_ref()
            .map(|f| crate::tree_panel::row_order(&f.layout, &app.model))
            .unwrap_or_default()
    };
    for (id, mods) in clicks {
        app.click_genome(id, mods, &order);
    }
    for c in cluster_clicks {
        app.set_cluster(c, true);
    }
    for (g, additive) in group_clicks {
        app.select_group(g, additive);
    }
    if let Some((on, n)) = request {
        app.tree.groups_on = on;
        app.tree.group_n = n;
        app.refresh_groups();
    }

    // Rebuild the figure only when what it shows changed.
    let groups_sig = if app.tree.groups_on {
        signature(
            app.groups
                .set
                .groups
                .iter()
                .flat_map(|g| g.members.iter().copied().chain([u32::MAX])),
        ) ^ u64::from(app.groups.set.threshold)
    } else {
        0
    };
    let key = (
        app.model.version,
        signature(app.selection.ids().iter().copied()),
        if all_open {
            0
        } else {
            app.expanded_signature()
        },
        groups_sig,
    );
    let stale = lock(&shared).built_from != Some(key);
    {
        let mut s = lock(&shared);
        s.groups = (app.tree.groups_on, app.tree.group_n);
        s.opts.groups = app.tree.groups_on;
    }
    if stale {
        let fig = build_figure(app, all_open);
        let mut s = lock(&shared);
        let rows_changed = match (&s.figure, &fig) {
            (Some(a), Some(b)) => (a.layout.rows - b.layout.rows).abs() > 1e-9,
            _ => true,
        };
        if rows_changed {
            s.view.needs_fit = true;
        }
        s.figure = fig;
        s.built_from = Some(key);
        s.export_dir = app.export_text.trim().to_string();
        ctx.request_repaint_of(viewport_id());
    }

    let builder = egui::ViewportBuilder::default()
        .with_title(format!(
            "Infinity-Win {} — tree ({} genomes)",
            infinity_core::VERSION,
            app.model.genome_count()
        ))
        .with_inner_size([1040.0, 880.0])
        .with_min_inner_size([480.0, 360.0]);
    let child = Arc::clone(&shared);
    ctx.show_viewport_deferred(viewport_id(), builder, move |ui, _class| {
        window_ui(ui, &child);
    });
}

/// The tree window's own pass.
fn window_ui(ui: &mut egui::Ui, shared: &SharedRef) {
    let ctx = ui.ctx().clone();
    if ctx.input(|i| i.viewport().close_requested()) {
        lock(shared).open = false;
        ctx.request_repaint_of(egui::ViewportId::ROOT);
        return;
    }
    // A light page whatever the main window's theme: the program's own light
    // style, for this window only.
    *ui.visuals_mut() = ctx.style_of(egui::Theme::Light).visuals.clone();
    let mut s = lock(shared);
    poll_png(&mut s);
    toolbar(ui, &mut s);
    canvas(ui, &mut s);
}

fn poll_png(s: &mut Shared) {
    let taken = s.png_result.lock().ok().and_then(|mut g| g.take());
    if let Some(r) = taken {
        s.status = match r {
            Ok(p) => format!("wrote {}", p.display()),
            Err(e) => format!("cannot write the picture: {e}"),
        };
    }
}

fn export_path(s: &Shared, file: &str) -> PathBuf {
    let dir = if s.export_dir.is_empty() {
        PathBuf::from(".")
    } else {
        PathBuf::from(&s.export_dir)
    };
    dir.join(file)
}

/// Write the whole figure as SVG, at a readable row height.
pub fn save_svg(fig: &Figure, opts: &PubOptions, path: &std::path::Path) -> Result<(), String> {
    let svg = crate::svg::render(
        fig,
        &crate::svg::SvgOptions {
            opts: *opts,
            tree_width: SVG_TREE_W,
            row_px: SVG_ROW_PX,
        },
    );
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    std::fs::write(path, svg).map_err(|e| format!("{}: {e}", path.display()))
}

fn toolbar(ui: &mut egui::Ui, s: &mut Shared) {
    egui::Panel::top("tree-window-toolbar")
        .frame(
            egui::Frame::new()
                .fill(Color32::from_rgb(246, 247, 249))
                .inner_margin(egui::Margin::symmetric(10, 6)),
        )
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                if ui.button("Fit").on_hover_text("fit the whole tree").clicked() {
                    s.view.needs_fit = true;
                }
                if ui.button("−").on_hover_text("zoom out").clicked() {
                    zoom(&mut s.view, 1.0 / 1.25);
                }
                if ui.button("+").on_hover_text("zoom in").clicked() {
                    zoom(&mut s.view, 1.25);
                }
                ui.separator();
                let before = s.all_open;
                ui.checkbox(&mut s.all_open, "Every cluster opened")
                    .on_hover_text("off: open the same clusters as the main window");
                if s.all_open != before {
                    s.built_from = None;
                    ui.ctx().request_repaint_of(egui::ViewportId::ROOT);
                }
                ui.menu_button("Show ▾", |ui| {
                    ui.checkbox(&mut s.opts.tip_snps, "SNPs at the tips (+n)");
                    ui.checkbox(&mut s.opts.branch_labels, "SNP counts on internal branches");
                    ui.checkbox(&mut s.opts.cluster_brackets, "Cluster brackets");
                    ui.checkbox(&mut s.opts.aligned_tips, "Aligned tip labels");
                    ui.checkbox(&mut s.opts.color_strip, "Cluster colour strip");
                    ui.checkbox(&mut s.opts.legend, "Scale bar legend");
                });
                ui.separator();
                let (mut on, mut n) = s.groups;
                let changed_on = ui
                    .checkbox(&mut on, "SNP groups ≤")
                    .on_hover_text("shade genomes linked by chains of at most this many SNPs")
                    .changed();
                let changed_n = ui
                    .add(egui::DragValue::new(&mut n).range(0..=100).suffix(" SNPs"))
                    .changed();
                if changed_on || changed_n {
                    s.groups = (on, n);
                    s.groups_request = Some((on, n));
                    ui.ctx().request_repaint_of(egui::ViewportId::ROOT);
                }
                ui.separator();
                if ui
                    .button("Save PNG")
                    .on_hover_text("this view, as on screen, to tree-window.png in the export folder")
                    .clicked()
                {
                    let p = export_path(s, "tree-window.png");
                    s.png_request = Some((p, true));
                }
                if ui
                    .button("Save SVG")
                    .on_hover_text("the whole tree as a vector figure, to tree-window.svg in the export folder")
                    .clicked()
                {
                    let p = export_path(s, "tree-window.svg");
                    s.status = match &s.figure {
                        Some(f) => match save_svg(f, &s.opts, &p) {
                            Ok(()) => format!("wrote {}", p.display()),
                            Err(e) => format!("cannot write the SVG: {e}"),
                        },
                        None => "no tree to save yet".to_string(),
                    };
                }
            });
            let hint = if s.status.is_empty() {
                format!(
                    "export folder: {} · drag to pan, scroll to zoom, click a name to select \
                     (Ctrl, Shift as in the main window), click a group label to select the group",
                    if s.export_dir.is_empty() { "." } else { &s.export_dir }
                )
            } else {
                s.status.clone()
            };
            ui.label(egui::RichText::new(hint).small().color(Color32::from_gray(90)));
        });
}

fn zoom(v: &mut PubView, f: f32) {
    v.px_per_row = (v.px_per_row * f).clamp(0.3, 200.0);
    v.px_per_unit = (v.px_per_unit * f).clamp(5.0, 60_000.0);
}

fn canvas(ui: &mut egui::Ui, s: &mut Shared) {
    egui::CentralPanel::default()
        .frame(egui::Frame::new().fill(pubtree::PAPER))
        .show(ui, |ui| {
            let rect = ui.max_rect();
            let resp = ui.allocate_rect(rect, Sense::click_and_drag());
            let ppp = ui.ctx().pixels_per_point();
            s.canvas = Some((rect, ppp));
            let painter = ui.painter().with_clip_rect(rect);
            let ink = Ink::PAPER;
            painter.rect_filled(rect, 0, ink.paper);
            let Some(fig) = s.figure.as_ref() else {
                painter.text(
                    rect.center(),
                    egui::Align2::CENTER_CENTER,
                    "No tree yet: it appears as soon as the genomes are clustered.",
                    egui::FontId::proportional(15.0),
                    ink.text_weak,
                );
                return;
            };
            let l = &fig.layout;

            // --- view: fit, pan, zoom ---------------------------------------
            let base = rect.min + Vec2::new(MARGIN, MARGIN);
            let geometry = |v: &PubView| Geometry {
                origin: base + v.offset + Vec2::new(0.0, v.px_per_row * 0.5),
                px_per_row: v.px_per_row,
                px_per_unit: v.px_per_unit,
                clip: rect,
            };
            let footer_h = |v: &PubView| {
                let g = geometry(v);
                let used = pubtree::Used {
                    exact: true,
                    estimate: true,
                    groups: s.opts.groups,
                    selection: !fig.marks.selection.is_empty(),
                    wedges: true,
                    tip_snps: s.opts.tip_snps,
                };
                pubtree::footer_height(l, &fig.marks, &g, &used, s.opts.legend, rect.width())
            };
            if s.view.needs_fit {
                let h = (rect.height() - 2.0 * MARGIN - footer_h(&s.view)).max(40.0);
                let ppr = (h / l.rows.max(1.0) as f32).clamp(0.3, 22.0);
                let room = pubtree::label_room(l, &fig.marks, &s.opts, ppr);
                let w = (rect.width() - 2.0 * MARGIN - room).max(rect.width() * 0.3);
                s.view.px_per_row = ppr;
                s.view.px_per_unit = (w / l.x_max.max(1e-6) as f32).clamp(5.0, 60_000.0);
                s.view.offset = Vec2::ZERO;
                s.view.needs_fit = false;
            }
            if resp.dragged() {
                s.view.offset += resp.drag_delta();
            }
            let pointer = resp.hover_pos();
            let scroll = ui.ctx().input(|i| i.smooth_scroll_delta);
            if scroll.y.abs() > 0.01 {
                let anchor = pointer.unwrap_or(rect.center());
                let f = (scroll.y * 0.004).exp();
                let o = geometry(&s.view).origin;
                let lx = (anchor.x - o.x) / s.view.px_per_unit;
                let ly = (anchor.y - o.y) / s.view.px_per_row;
                zoom(&mut s.view, f);
                let o2 = Pos2::new(
                    anchor.x - lx * s.view.px_per_unit,
                    anchor.y - ly * s.view.px_per_row,
                );
                s.view.offset = o2 - base - Vec2::new(0.0, s.view.px_per_row * 0.5);
            }
            if scroll.x.abs() > 0.01 {
                s.view.offset.x += scroll.x;
            }
            let g = geometry(&s.view);

            // --- draw --------------------------------------------------------
            let sc = pubtree::scene(l, &fig.marks, &g, &s.opts, &ink, None);
            pubtree::paint(&painter, &sc.prims);
            let fh =
                pubtree::footer_height(l, &fig.marks, &g, &sc.used, s.opts.legend, rect.width());
            let area = Rect::from_min_max(
                Pos2::new(rect.left(), rect.bottom() - fh - 4.0),
                rect.right_bottom(),
            );
            let footer =
                pubtree::footer(l, &fig.marks, &g, &sc.used, s.opts.legend, &ink, area, true);
            pubtree::paint(&painter, &footer);

            // --- clicks ------------------------------------------------------
            let hit = pointer.filter(|p| !area.contains(*p)).and_then(|p| {
                sc.hits
                    .iter()
                    .rev()
                    .find(|(_, r)| r.contains(p))
                    .map(|(h, _)| *h)
            });
            if hit.is_some() {
                ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
            }
            if resp.clicked() {
                let mods = ui.ctx().input(|i| ClickMods::from_modifiers(&i.modifiers));
                match hit {
                    Some(Hit::Genome(id)) if id != u32::MAX => {
                        s.clicks.push((id, mods));
                        ui.ctx().request_repaint_of(egui::ViewportId::ROOT);
                    }
                    Some(Hit::Cluster(id)) | Some(Hit::Badge(id)) => {
                        s.cluster_clicks.push(id);
                        ui.ctx().request_repaint_of(egui::ViewportId::ROOT);
                    }
                    Some(Hit::Group(i)) => {
                        s.group_clicks.push((i, mods.toggle));
                        ui.ctx().request_repaint_of(egui::ViewportId::ROOT);
                    }
                    _ => {}
                }
            }

            // --- a PNG of this window, asked for by the button ---------------
            if let Some((path, crop)) = s.png_request.take() {
                let slot = Arc::clone(&s.png_result);
                let region = crop.then_some((rect, ppp));
                let wake = ui.ctx().clone();
                ui.ctx().request_screenshot(move |img| {
                    let res = crate::export::write_png(&path, &img, region)
                        .map(|()| path.clone())
                        .map_err(|e| e.to_string());
                    if let Ok(mut g) = slot.lock() {
                        *g = Some(res);
                    }
                    // Show "wrote …" without waiting for the next event.
                    wake.request_repaint_of(viewport_id());
                });
                ui.ctx().request_repaint();
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signatures_change_with_the_ids_and_their_order() {
        assert_eq!(signature([1, 2, 3]), signature([1, 2, 3]));
        assert_ne!(signature([1, 2, 3]), signature([1, 2]));
        assert_ne!(signature([1, 2]), signature([2, 1]));
        assert_ne!(signature([]), signature([0]));
    }

    #[test]
    fn the_svg_of_the_window_is_written_from_its_primitives() {
        let dir = std::env::temp_dir().join(format!("infinity-tw-{}", std::process::id()));
        let path = dir.join("tree-window.svg");
        let f = crate::pubtree::tests_support::figure();
        save_svg(&f, &PubOptions::default(), &path).unwrap();
        let svg = std::fs::read_to_string(&path).unwrap();
        assert!(svg.starts_with("<svg xmlns"));
        assert!(svg.contains(r##"fill="#ffffff""##), "white page");
        assert!(svg.contains(">g1<tspan"));
        assert!(
            svg.contains(r#"font-weight="bold">g1<tspan"#),
            "g1 is selected"
        );
        assert!(svg.contains("SNPs between clusters"));
        assert!(svg.contains("G1 · 2 genomes"));
        assert!(svg.trim_end().ends_with("</svg>"));
        std::fs::remove_dir_all(&dir).ok();
    }
}
