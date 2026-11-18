//! The genome list (section 8): name, cluster colour, length, contigs, focus.
//! Sortable, filterable, and its selection is what links the list to the tree
//! and to the map.
//!
//! Rows are virtualised with `ScrollArea::show_rows`, so a thousand genomes cost
//! the same per frame as twenty: only the visible ones are built.

use eframe::egui;
use egui::{Align, Layout, RichText, Sense, Vec2};

use crate::app::MtatoApp;
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

pub fn show(ui: &mut egui::Ui, app: &mut MtatoApp) {
    egui::Panel::left("genomes")
        .default_size(300.0)
        .min_size(210.0)
        .resizable(true)
        .show(ui, |ui| {
            header(ui, app);
            list(ui, app);
        });
}

fn header(ui: &mut egui::Ui, app: &mut MtatoApp) {
    let mut found: Option<u32> = None;
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
    if let Some(id) = found {
        app.select_genome(id, false);
    }
}

fn column_headers(ui: &mut egui::Ui, app: &mut MtatoApp) {
    let mut clicked: Option<SortColumn> = None;
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        header_button(ui, &mut clicked, app.sort, SortColumn::Id, "#", W_ID);
        ui.add_sized(
            Vec2::new(W_SWATCH, ROW_H),
            egui::Label::new(RichText::new("c").weak()),
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
        ui.add_sized(Vec2::new(width, ROW_H), egui::Button::new(text).frame(false))
    } else {
        ui.add(egui::Button::new(text).frame(false))
    };
    if resp.clicked() {
        *clicked = Some(column);
    }
    resp.on_hover_text(format!("sort by {label}"));
}

fn list(ui: &mut egui::Ui, app: &mut MtatoApp) {
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
    if !(width > 0.0) {
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

    /// A painter that needs no window, so text metrics can be tested headless.
    fn painter() -> egui::Painter {
        egui::Context::default().debug_painter()
    }

    #[test]
    fn a_short_name_is_left_alone() {
        let p = painter();
        let font = egui::FontId::proportional(11.5);
        assert_eq!(clip(&p, "strainA01", 200.0, font.clone()), "strainA01");
    }

    #[test]
    fn a_long_name_is_cut_to_fit_and_ends_in_an_ellipsis() {
        let p = painter();
        let font = egui::FontId::proportional(11.5);
        let text = "Salmonella_enterica_serovar_Typhimurium_strain_LT2_assembly_0042";
        let width = 80.0;
        let out = clip(&p, text, width, font.clone());
        assert!(out.ends_with('…'), "{out}");
        assert!(out.len() < text.len());
        let w = p
            .layout_no_wrap(out, font, egui::Color32::BLACK)
            .size()
            .x;
        assert!(w <= width + 0.5, "the clipped label is {w} px, over {width}");
    }

    #[test]
    fn a_width_that_cannot_hold_anything_yields_nothing_rather_than_a_overflow() {
        let p = painter();
        let font = egui::FontId::proportional(11.5);
        assert_eq!(clip(&p, "anything", 0.0, font.clone()), "");
        assert_eq!(clip(&p, "anything", -5.0, font), "");
        // One pixel is too narrow for even a single character plus ellipsis, so
        // the shortest possible answer is returned rather than looping.
        assert!(clip(&p, "anything", 1.0, font).chars().count() <= 2);
    }
}


