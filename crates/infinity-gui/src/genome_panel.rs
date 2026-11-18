//! The genome list (section 8): name, cluster colour, length, contigs, focus.
//! Sortable, filterable, and its selection is what links the list to the tree
//! and to the map.
//!
//! Rows are virtualised with `ScrollArea::show_rows`, so a thousand genomes cost
//! the same per frame as twenty: only the visible ones are built.

use eframe::egui;
use egui::{Align, Layout, RichText, Sense, Vec2};

use crate::app::InfinityApp;
use crate::state::{sorted_indices, SortColumn};
use crate::{colors, fmt};

/// Row height in pixels; also the virtualisation unit.
const ROW_H: f32 = 16.0;
/// Column widths of the fixed columns.
const W_ID: f32 = 34.0;
const W_SWATCH: f32 = 14.0;
const W_LEN: f32 = 58.0;
const W_CONTIGS: f32 = 34.0;
const W_FOCUS: f32 = 40.0;

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
    ui.horizontal(|ui| {
        ui.label(
            RichText::new(format!(
                "GENOMES ({})",
                fmt::count(app.model.genome_count() as u64)
            ))
            .strong(),
        );
        if !app.selection.is_empty() {
            ui.weak(format!("{} selected", app.selection.len()));
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
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.filter)
                .desired_width(f32::INFINITY)
                .hint_text("filter by name"),
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
                ui.weak("found:");
                if ui.button(name.as_str()).clicked() {
                    found = Some(id);
                }
            });
        }
    }
    column_headers(ui, app);
    ui.separator();
    if clear_sel {
        app.clear_selection();
    }
    if let Some(id) = found {
        app.select_genome(id, false);
    }
}

fn column_headers(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let mut clicked: Option<SortColumn> = None;
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        header_button(ui, &mut clicked, app.sort, SortColumn::Id, "#", W_ID);
        // The cluster swatch column doubles as the "group by cluster" sort.
        header_button(
            ui,
            &mut clicked,
            app.sort,
            SortColumn::Cluster,
            "c",
            W_SWATCH,
        );
        header_button(
            ui,
            &mut clicked,
            app.sort,
            SortColumn::Name,
            "name",
            f32::INFINITY,
        );
        header_button(
            ui,
            &mut clicked,
            app.sort,
            SortColumn::Length,
            "length",
            W_LEN,
        );
        header_button(
            ui,
            &mut clicked,
            app.sort,
            SortColumn::Contigs,
            "ctg",
            W_CONTIGS,
        );
        header_button(
            ui,
            &mut clicked,
            app.sort,
            SortColumn::Focus,
            "focus",
            W_FOCUS,
        );
        // N50 and GC are in the model but not worth a column of their own; they
        // stay reachable as sort keys.
        ui.menu_button("⋮", |ui| {
            if ui.button("sort by N50").clicked() {
                clicked = Some(SortColumn::N50);
                ui.close();
            }
            if ui.button("sort by GC").clicked() {
                clicked = Some(SortColumn::Gc);
                ui.close();
            }
        });
    });
    if let Some(c) = clicked {
        app.sort = app.sort.toggled(c);
    }
}

fn header_button(
    ui: &mut egui::Ui,
    clicked: &mut Option<SortColumn>,
    sort: crate::state::SortSpec,
    column: SortColumn,
    label: &str,
    width: f32,
) {
    let arrow = if sort.column != column {
        ""
    } else if sort.ascending {
        " ↑"
    } else {
        " ↓"
    };
    let mut text = RichText::new(format!("{label}{arrow}")).small();
    if sort.column == column {
        text = text.strong();
    }
    let resp = if width.is_finite() {
        ui.add_sized(
            Vec2::new(width, ROW_H),
            egui::Button::new(text).frame(false),
        )
    } else {
        ui.add(egui::Button::new(text).frame(false))
    };
    if resp.clicked() {
        *clicked = Some(column);
    }
    resp.on_hover_text(format!("sort by {label}"));
}

fn list(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let order = sorted_indices(&app.model.genomes, &app.filter, &app.sort);
    let selected = app.selection.clone();
    let dark = ui.visuals().dark_mode;
    let total = order.len();

    if total == 0 {
        ui.vertical_centered(|ui| {
            ui.add_space(24.0);
            ui.weak(if app.model.genome_count() == 0 {
                "no genomes yet: add a file or a folder above"
            } else {
                "nothing matches the filter"
            });
        });
        return;
    }

    // Filled in by the row closure, applied afterwards: the closure only needs
    // to read the model, never to mutate the app.
    let mut click: Option<(u32, bool, bool)> = None;

    let model = &app.model;
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show_rows(ui, ROW_H, total, |ui, range| {
            ui.set_min_width(ui.available_width());
            for i in range {
                let g = &model.genomes[order[i]];
                let is_selected = selected.contains(&g.id);
                let colour = g
                    .cluster
                    .map(|c| colors::cluster_color(c, dark))
                    .unwrap_or_else(|| ui.visuals().weak_text_color());
                let (rect, resp) =
                    ui.allocate_exact_size(Vec2::new(ui.available_width(), ROW_H), Sense::click());
                let painter = ui.painter();
                if is_selected {
                    painter.rect_filled(rect, 0, colors::wash(colors::SELECTED, 60));
                } else if resp.hovered() {
                    painter.rect_filled(rect, 0, ui.visuals().faint_bg_color);
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

                let name_w = (rect.width() - W_ID - W_SWATCH - W_LEN - W_CONTIGS - W_FOCUS - 12.0)
                    .max(40.0);
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
                let focus_text = format!(
                    "F{}{}",
                    g.focus,
                    if g.approximate { "≈" } else { "" }
                );
                painter.text(
                    pos(x + W_FOCUS - 4.0, rect.center().y),
                    egui::Align2::RIGHT_CENTER,
                    focus_text,
                    egui::FontId::proportional(10.0),
                    colors::focus_color(g.focus, dark),
                );

                if resp.clicked() {
                    let additive = ui.ctx().input(|i| i.modifiers.command);
                    click = Some((g.id, additive, false));
                }
                if resp.double_clicked() {
                    click = Some((g.id, false, true));
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
                        ui.weak("click to select, ctrl-click for a pair, double-click to open its cluster");
                    });
                });
            }
        });

    if let Some((id, additive, open_cluster)) = click {
        app.select_genome(id, additive);
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
