//! The genome list (section 8): name, cluster colour, length, contigs, focus.
//! Sortable, filterable, and its selection is what links the list to the tree
//! and to the map.
//!
//! Rows are virtualised with `ScrollArea::show_rows`, so a thousand genomes cost
//! the same per frame as twenty: only the visible ones are built.

use eframe::egui;
use egui::{Align, Layout, RichText, Sense, Vec2};

use crate::app::InfinityApp;
use crate::selection::ClickMods;
use crate::state::{sorted_indices, SortColumn};
use crate::{colors, fmt, theme};

/// Row height in pixels; also the virtualisation unit.
const ROW_H: f32 = 16.0;
/// Column widths of the fixed columns.
const W_ID: f32 = 30.0;
const W_SWATCH: f32 = 14.0;
const W_LEN: f32 = 54.0;
const W_CONTIGS: f32 = 30.0;
const W_FOCUS: f32 = 36.0;

/// `total_w` is the width of the whole window below the toolbar: the list may
/// take at most about a third of it, so the tree always keeps room to draw.
pub fn show(ui: &mut egui::Ui, app: &mut InfinityApp, total_w: f32) {
    let max_w = (total_w * 0.34).max(200.0);
    egui::Panel::left("genomes")
        .default_size(260.0)
        .size_range(190.0..=max_w)
        .resizable(true)
        .show(ui, |ui| {
            header(ui, app);
            list(ui, app);
        });
}

fn header(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let mut found: Option<u32> = None;
    let mut clear_sel = false;
    let mut select_all = false;
    let mut sort_by: Option<SortColumn> = None;
    let p = theme::palette(ui);
    theme::panel_header(ui, "Genomes", |ui| {
        ui.label(RichText::new(fmt::count(app.model.genome_count() as u64)).color(p.text_weak));
        if app.model.genome_count() > 0
            && app.selection.len() < app.model.genome_count()
            && ui
                .small_button("all")
                .on_hover_text("select every genome the list shows (Ctrl+A)")
                .clicked()
        {
            select_all = true;
        }
        if !app.selection.is_empty() {
            ui.label(RichText::new(format!("· {} selected", app.selection.len())).color(p.accent));
            if ui
                .small_button("clear")
                .on_hover_text("clear the selection")
                .clicked()
            {
                clear_sel = true;
            }
        }
    });
    ui.horizontal(|ui| {
        // N50 and GC are in the model but not worth a column of their own;
        // they stay reachable as sort keys.
        ui.menu_button("⇅", |ui| {
            if ui.button("sort by N50").clicked() {
                sort_by = Some(SortColumn::N50);
                ui.close();
            }
            if ui.button("sort by GC").clicked() {
                sort_by = Some(SortColumn::Gc);
                ui.close();
            }
        })
        .response
        .on_hover_text("more sort orders: N50, GC");
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.filter)
                .desired_width(f32::INFINITY)
                .hint_text("🔍 filter by name"),
        );
        if resp.changed() {
            // Keep the tree's search box and the list filter in step: both are
            // "find a genome by name".
            app.tree.search = app.filter.clone();
            app.tree.search_hit = app.model.search(&app.filter);
        }
        if !app.filter.is_empty() && ui.small_button("✕").clicked() {
            app.filter.clear();
            app.tree.search.clear();
            app.tree.search_hit = None;
        }
    });
    if let Some(id) = app.tree.search_hit {
        if let Some(g) = app.model.genome(id) {
            let name = g.name.clone();
            ui.horizontal(|ui| {
                theme::hint(ui, "found:");
                if ui.small_button(name.as_str()).clicked() {
                    found = Some(id);
                }
            });
        }
    }
    column_headers(ui, app);
    if clear_sel {
        app.clear_selection();
    }
    if select_all {
        app.select_all_visible();
    }
    if let Some(id) = found {
        app.select_genome(id, false);
    }
    if let Some(c) = sort_by {
        app.sort = app.sort.toggled(c);
    }
}

/// The x positions and widths of the columns for a row `width` wide. The rows
/// and the header share this, so they always line up, and nothing in either is
/// ever wider than the panel.
fn columns(left: f32, width: f32) -> [(SortColumn, &'static str, f32, f32); 6] {
    let name_w = (width - W_ID - W_SWATCH - W_LEN - W_CONTIGS - W_FOCUS - 12.0).max(40.0);
    let mut x = left + 2.0;
    let mut out = [(SortColumn::Id, "", 0.0, 0.0); 6];
    for (i, (col, label, w)) in [
        (SortColumn::Id, "#", W_ID),
        (SortColumn::Cluster, "c", W_SWATCH),
        (SortColumn::Name, "name", name_w),
        (SortColumn::Length, "length", W_LEN),
        (SortColumn::Contigs, "ctg", W_CONTIGS),
        (SortColumn::Focus, "focus", W_FOCUS),
    ]
    .into_iter()
    .enumerate()
    {
        out[i] = (col, label, x, w);
        x += w;
    }
    out
}

/// One painted header row. Each cell is a click target that sorts by its
/// column; clicking the sorted column again flips the direction.
fn column_headers(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let p = theme::palette(ui);
    let (rect, _) = ui.allocate_exact_size(Vec2::new(ui.available_width(), 18.0), Sense::hover());
    let painter = ui.painter_at(rect);
    painter.line_segment(
        [rect.left_bottom(), rect.right_bottom()],
        egui::Stroke::new(1.0, p.border),
    );
    let mut clicked: Option<SortColumn> = None;
    for (col, label, x, w) in columns(rect.left(), rect.width()) {
        let cell =
            egui::Rect::from_min_size(egui::pos2(x, rect.top()), Vec2::new(w, rect.height()));
        let resp = ui
            .interact(cell, ui.scope_id().with(("sort", label)), Sense::click())
            .on_hover_text(format!("sort by {label}"));
        let sorted = app.sort.column == col;
        let arrow = match (sorted, app.sort.ascending) {
            (false, _) => "",
            (true, true) => " ↑",
            (true, false) => " ↓",
        };
        let colour = if sorted {
            p.accent
        } else if resp.hovered() {
            p.text_strong
        } else {
            p.text_weak
        };
        // Numbers are right-aligned in the rows, so their headers are too.
        let (anchor, at) = match col {
            SortColumn::Length | SortColumn::Contigs | SortColumn::Focus => (
                egui::Align2::RIGHT_CENTER,
                egui::pos2(x + w - 4.0, cell.center().y),
            ),
            _ => (egui::Align2::LEFT_CENTER, egui::pos2(x, cell.center().y)),
        };
        painter.text(
            at,
            anchor,
            format!("{label}{arrow}"),
            egui::FontId::proportional(10.5),
            colour,
        );
        if resp.clicked() {
            clicked = Some(col);
        }
    }
    if let Some(c) = clicked {
        app.sort = app.sort.toggled(c);
    }
}

fn list(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let order = sorted_indices(&app.model.genomes, &app.filter, &app.sort);
    let selected = &app.selection;
    let dark = ui.visuals().dark_mode;
    let total = order.len();

    if total == 0 {
        if app.model.genome_count() == 0 {
            empty_state(ui, app);
        } else {
            ui.vertical_centered(|ui| {
                ui.add_space(24.0);
                ui.weak("nothing matches the filter");
            });
        }
        return;
    }

    // Filled in by the row closure, applied afterwards: the closure only needs
    // to read the model, never to mutate the app.
    let mut click: Option<(u32, ClickMods, bool)> = None;

    let model = &app.model;
    let accent = theme::palette(ui).accent;
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show_rows(ui, ROW_H, total, |ui, range| {
            ui.set_min_width(ui.available_width());
            for i in range {
                let g = &model.genomes[order[i]];
                let is_selected = selected.contains(g.id);
                let colour = g
                    .cluster
                    .map(|c| colors::cluster_color(c, dark))
                    .unwrap_or_else(|| ui.visuals().weak_text_color());
                let (rect, resp) =
                    ui.allocate_exact_size(Vec2::new(ui.available_width(), ROW_H), Sense::click());
                let painter = ui.painter();
                if is_selected {
                    painter.rect_filled(rect, 3, theme::with_alpha(accent, 48));
                    painter.rect_filled(
                        egui::Rect::from_min_size(rect.min, Vec2::new(2.0, rect.height())),
                        0,
                        accent,
                    );
                } else if resp.hovered() {
                    painter.rect_filled(rect, 3, ui.visuals().faint_bg_color);
                } else if i % 2 == 1 {
                    // A faint zebra stripe keeps a long list readable.
                    painter.rect_filled(
                        rect,
                        0,
                        theme::with_alpha(ui.visuals().faint_bg_color, 120),
                    );
                }
                let mut x = rect.left() + 2.0;
                painter.text(
                    pos(x, rect.center().y),
                    egui::Align2::LEFT_CENTER,
                    g.id.to_string(),
                    egui::FontId::proportional(10.0),
                    ui.visuals().weak_text_color(),
                );
                x += W_ID;
                painter.circle_filled(pos(x + 5.0, rect.center().y), 4.5, colour);
                x += W_SWATCH;

                let name_w =
                    (rect.width() - W_ID - W_SWATCH - W_LEN - W_CONTIGS - W_FOCUS - 12.0).max(40.0);
                let name_colour = if g.is_pivot {
                    ui.visuals().strong_text_color()
                } else {
                    ui.visuals().text_color()
                };
                let name = if g.is_pivot {
                    format!("◆ {}", g.name)
                } else {
                    g.name.clone()
                };
                painter.text(
                    pos(x, rect.center().y),
                    egui::Align2::LEFT_CENTER,
                    clip(painter, &name, name_w, egui::FontId::proportional(11.5)),
                    egui::FontId::proportional(11.5),
                    name_colour,
                );
                x += name_w;
                painter.text(
                    pos(x + W_LEN - 4.0, rect.center().y),
                    egui::Align2::RIGHT_CENTER,
                    fmt::bp(g.length),
                    egui::FontId::proportional(10.0),
                    ui.visuals().weak_text_color(),
                );
                x += W_LEN;
                painter.text(
                    pos(x + W_CONTIGS - 4.0, rect.center().y),
                    egui::Align2::RIGHT_CENTER,
                    g.contigs.to_string(),
                    egui::FontId::proportional(10.0),
                    ui.visuals().weak_text_color(),
                );
                x += W_CONTIGS;
                // The focus column carries the "≈" marker whenever the numbers
                // behind it are estimates (section 0 rule 8).
                let focus_text = format!("F{}{}", g.focus, if g.approximate { "≈" } else { "" });
                painter.text(
                    pos(x + W_FOCUS - 4.0, rect.center().y),
                    egui::Align2::RIGHT_CENTER,
                    focus_text,
                    egui::FontId::proportional(10.0),
                    colors::focus_color(g.focus, dark),
                );

                if resp.clicked() {
                    let mods = ui.ctx().input(|i| ClickMods::from_modifiers(&i.modifiers));
                    click = Some((g.id, mods, false));
                }
                if resp.double_clicked() {
                    click = Some((g.id, ClickMods::PLAIN, true));
                }
                resp.on_hover_ui(|ui| {
                    ui.with_layout(Layout::top_down(Align::LEFT), |ui| {
                        ui.label(RichText::new(&g.name).strong());
                        ui.label(format!(
                            "length {} | contigs {} | N50 {} | GC {}",
                            fmt::bp(g.length),
                            fmt::count(u64::from(g.contigs)),
                            fmt::bp(g.n50),
                            fmt::pct(Some(g.gc))
                        ));
                        ui.label(format!(
                            "cluster {} | {}",
                            g.cluster
                                .map(|c| c.to_string())
                                .unwrap_or_else(|| "-".to_string()),
                            if g.is_pivot { "pivot" } else { "member" }
                        ));
                        ui.label(fmt::focus_label(g.focus, g.approximate));
                        if let Some(n) = &g.numbers {
                            ui.separator();
                            ui.label(format!(
                                "{}diff pieces {}",
                                fmt::approx_mark(g.approximate),
                                fmt::count(n.diff_pieces)
                            ));
                            ui.label(format!(
                                "{}SNPs {}  indels {}",
                                fmt::approx_mark(g.approximate),
                                fmt::count(n.snps),
                                fmt::count(n.indels)
                            ));
                            ui.label(format!(
                                "uncovered {} | novel {}",
                                fmt::bp(n.uncovered_bp),
                                fmt::bp(n.novel_bp)
                            ));
                        } else {
                            ui.weak("no numbers yet");
                        }
                        ui.separator();
                        ui.weak(
                            "click to select · Ctrl-click to add or remove · Shift-click for a \
                             range · Ctrl+A for all · double-click to open its cluster",
                        );
                    });
                });
            }
        });

    if let Some((id, mods, open_cluster)) = click {
        let ids: Vec<u32> = order.iter().map(|&i| app.model.genomes[i].id).collect();
        app.click_genome(id, mods, &ids);
        if open_cluster {
            if let Some(c) = app.model.cluster_of(id) {
                if !app.expanded.contains(&c) {
                    app.toggle_cluster(c);
                } else {
                    app.set_cluster(c, true);
                }
            }
        }
    }
}

/// What the list shows before any run: the two ways in, as big buttons, and
/// what is queued for the next run.
fn empty_state(ui: &mut egui::Ui, app: &mut InfinityApp) {
    use crate::dialogs::PickKind;
    let p = theme::palette(ui);
    let busy = app.pending_pick.is_some();
    ui.vertical_centered(|ui| {
        ui.add_space(28.0);
        ui.label(
            RichText::new("Add your genomes")
                .size(15.0)
                .strong()
                .color(p.text_strong),
        );
        ui.add_space(4.0);
        theme::hint(ui, "assembled genomes, one FASTA file each");
        ui.add_space(10.0);
        let w = (ui.available_width() - 24.0).clamp(120.0, 200.0);
        if ui
            .add_enabled(
                !busy,
                theme::primary_button("＋ Add files…", &p).min_size(Vec2::new(w, 28.0)),
            )
            .on_hover_text("choose FASTA files (.fasta, .fa, .fna, .fas, also .gz)")
            .clicked()
        {
            app.request_pick(PickKind::Files);
        }
        ui.add_space(4.0);
        if ui
            .add_enabled(
                !busy,
                egui::Button::new("＋ Add folder…").min_size(Vec2::new(w, 28.0)),
            )
            .on_hover_text("every FASTA file directly inside the folder is one genome")
            .clicked()
        {
            app.request_pick(PickKind::Folders);
        }
        ui.add_space(10.0);
        theme::hint(ui, "or drop files and folders anywhere on the window");
        if !app.inputs.is_empty() {
            ui.add_space(14.0);
            ui.label(
                RichText::new(format!(
                    "Ready to run: {} in {} item(s)",
                    crate::dialogs::count_text(Some(app.input_genome_total())),
                    app.inputs.len()
                ))
                .color(p.accent),
            );
            theme::hint(ui, "press ▶ Run in the toolbar");
        }
    });
}

fn pos(x: f32, y: f32) -> egui::Pos2 {
    egui::pos2(x, y)
}

/// Truncate a label to a width, with an ellipsis.
///
/// Measured once, cut proportionally, then walked back only if the guess was too
/// generous: the list is redrawn every frame, so a character-by-character scan
/// would cost more than the rest of the row together.
pub(crate) fn clip(painter: &egui::Painter, text: &str, width: f32, font: egui::FontId) -> String {
    // `width` comes from a panel the user can resize, so NaN and negatives are
    // both possible and both mean "no room at all".
    if !width.is_finite() || width <= 0.0 {
        return String::new();
    }
    let measure = |s: &str| {
        painter
            .layout_no_wrap(s.to_string(), font.clone(), egui::Color32::BLACK)
            .size()
            .x
    };
    let full = measure(text);
    if full <= width {
        return text.to_string();
    }
    let chars: Vec<char> = text.chars().collect();
    let mut keep = ((chars.len() as f32 * width / full).floor() as usize)
        .saturating_sub(1)
        .clamp(1, chars.len());
    loop {
        let mut out: String = chars.iter().take(keep).collect();
        out.push('…');
        if measure(&out) <= width || keep <= 1 {
            return out;
        }
        keep -= 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Run `f` against a real painter with no window. Fonts only exist inside a
    /// pass, so the callback runs inside `Context::run_ui`.
    fn with_painter<R>(f: impl FnOnce(&egui::Painter) -> R) -> R {
        let ctx = egui::Context::default();
        ctx.set_fonts(egui::FontDefinitions::default());
        // `run_ui` takes an `FnMut`, so the one-shot callback is handed over
        // through an `Option`.
        let mut f = Some(f);
        let mut out: Option<R> = None;
        let full = ctx.run_ui(egui::RawInput::default(), |ui| {
            if let Some(f) = f.take() {
                out = Some(f(ui.painter()));
            }
        });
        full.drop_without_applying_deltas();
        out.expect("the callback ran")
    }

    #[test]
    fn a_short_name_is_left_alone() {
        let out = with_painter(|p| {
            let font = egui::FontId::proportional(11.5);
            clip(p, "strainA01", 200.0, font)
        });
        assert_eq!(out, "strainA01");
    }

    #[test]
    fn a_long_name_is_cut_to_fit_and_ends_in_an_ellipsis() {
        let text = "Salmonella_enterica_serovar_Typhimurium_strain_LT2_assembly_0042";
        let width = 80.0;
        let out = with_painter(|p| {
            let font = egui::FontId::proportional(11.5);
            let cut = clip(p, text, width, font.clone());
            let w = p
                .layout_no_wrap(cut.clone(), font, egui::Color32::BLACK)
                .size()
                .x;
            (cut, w)
        });
        assert!(out.0.ends_with('…'), "{:?}", out.0);
        assert!(out.0.len() < text.len());
        assert!(
            out.1 <= width + 0.5,
            "the clipped label is {} px, over {width}",
            out.1
        );
        // And it kept a real prefix of the name, not a mangled middle.
        assert!(text.starts_with(out.0.trim_end_matches('…')));
    }

    #[test]
    fn a_width_that_cannot_hold_anything_yields_nothing_rather_than_an_overflow() {
        with_painter(|p| {
            let font = egui::FontId::proportional(11.5);
            assert_eq!(clip(p, "anything", 0.0, font.clone()), "");
            assert_eq!(clip(p, "anything", -5.0, font.clone()), "");
            assert_eq!(clip(p, "anything", f32::NAN, font.clone()), "");
            // One pixel is too narrow for even a single character plus ellipsis,
            // so the shortest possible answer comes back rather than a loop.
            assert!(clip(p, "anything", 1.0, font).chars().count() <= 2);
        });
    }

    #[test]
    fn a_wide_column_keeps_the_whole_name_even_when_it_is_long() {
        let text = "a_really_quite_long_genome_name_that_goes_on";
        let out = with_painter(|p| clip(p, text, 1000.0, egui::FontId::proportional(11.5)));
        assert_eq!(out, text);
    }
}
