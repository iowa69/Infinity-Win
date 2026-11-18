//! The right-hand "Run" panel: the live counters of section 7, the run's result,
//! the warnings, and the exports.
//!
//! The parameters that used to live here are in the Settings window now (see
//! `settings_panel`), so this panel is about the run in front of the user, not
//! about the next one. It can be closed from the toolbar (or dragged shut) and
//! its contents scroll rather than overflow, so at the smallest window size it
//! costs the tree nothing but its own width.

use eframe::egui;
use egui::RichText;

use crate::app::InfinityApp;
use crate::chrome::ExportAction;
use crate::export::Item;
use crate::state::Note;
use crate::{fmt, theme};

pub fn show(ui: &mut egui::Ui, app: &mut InfinityApp, total_w: f32) {
    let max_w = (total_w * 0.30).max(220.0);
    let mut open = app.show_side;
    egui::Panel::right("run-panel")
        .default_size(250.0)
        .size_range(210.0..=max_w)
        .resizable(true)
        .show_collapsible(ui, &mut open, |ui| {
            theme::panel_header(ui, "Run", |ui| {
                let p = theme::palette(ui);
                let state = if app.running() {
                    "running"
                } else if app.model.finished {
                    "finished"
                } else if app.model.genome_count() > 0 {
                    "stopped"
                } else {
                    "idle"
                };
                ui.label(RichText::new(state).small().color(p.text_weak));
            });
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    counters(ui, app);
                    summary(ui, app);
                    notes(ui, app);
                    export_area(ui, app);
                });
        });
    app.show_side = open;
}

fn counters(ui: &mut egui::Ui, app: &InfinityApp) {
    theme::section(ui, "Live counters");
    let c = &app.model.counters;
    // "blocks identical 99.7 %" is the number the whole idea is judged on, so it
    // goes first and gets a bar.
    let identity = c.block_identity();
    theme::stat_row(ui, "blocks identical", &fmt::pct(identity));
    if let Some(f) = identity {
        ui.add(egui::ProgressBar::new(f as f32).desired_height(6.0));
    }
    theme::stat_row(ui, "blocks looked up", &fmt::count(c.blocks_looked_up));
    theme::stat_row(
        ui,
        "matched L1/L2/L3",
        &format!(
            "{} / {} / {}",
            fmt::count(c.blocks_matched[0]),
            fmt::count(c.blocks_matched[1]),
            fmt::count(c.blocks_matched[2])
        ),
    );
    theme::stat_row(ui, "pieces compared", &fmt::count(c.pieces_compared));
    theme::stat_row(ui, "cache hits", &fmt::count(c.cache_hits));
    theme::stat_row(ui, "alignments run", &fmt::count(c.alignments_run));
    theme::stat_row(ui, "bases aligned", &fmt::count(c.bases_aligned));
    theme::stat_row(ui, "k-mers scanned", &fmt::count(c.kmers_scanned));
    theme::stat_row(ui, "genomes scanned", &fmt::count(c.genomes_scanned));
    theme::stat_row(ui, "end extensions", &fmt::count(c.end_extensions));
    theme::stat_row(
        ui,
        "stray points dropped",
        &fmt::count(c.stray_points_dropped),
    );
    theme::stat_row(ui, "breakpoints", &fmt::count(c.breakpoints));
}

fn summary(ui: &mut egui::Ui, app: &InfinityApp) {
    let Some(s) = &app.model.summary else {
        return;
    };
    theme::section(ui, if s.cancelled { "Stopped" } else { "Result" });
    theme::stat_row(ui, "genomes", &fmt::count(u64::from(s.genomes)));
    theme::stat_row(ui, "clusters", &fmt::count(u64::from(s.clusters)));
    theme::stat_row(
        ui,
        "focus 0/1/2",
        &format!(
            "{}/{}/{}",
            s.focus_counts[0], s.focus_counts[1], s.focus_counts[2]
        ),
    );
    theme::stat_row(
        ui,
        "first tree",
        &s.timings
            .first_tree_ms
            .map(fmt::secs)
            .unwrap_or_else(|| "not reached".to_string()),
    );
    theme::stat_row(
        ui,
        "exact tree",
        &s.timings
            .exact_tree_ms
            .map(fmt::secs)
            .unwrap_or_else(|| "not reached".to_string()),
    );
    if let Some(t) = s.timings.total_ms {
        theme::stat_row(ui, "total", &fmt::secs(t));
    }
    match s.peak_memory_bytes {
        // An unknown peak is reported as unknown, never as zero.
        Some(b) => theme::stat_row(ui, "peak RSS", &fmt::bytes(b)),
        None => theme::stat_row(ui, "peak RSS", "unknown"),
    }
    for n in &s.notes {
        theme::hint(ui, n);
    }
}

fn notes(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let count = app.model.notes.len();
    theme::section(ui, &format!("Warnings and notes ({count})"));
    if count == 0 {
        theme::hint(ui, "none");
        return;
    }
    ui.checkbox(&mut app.show_notes, "show");
    if !app.show_notes {
        return;
    }
    let notes: Vec<Note> = app.model.notes.clone();
    egui::ScrollArea::vertical()
        .id_salt("notes-scroll")
        .max_height(160.0)
        .auto_shrink([false, true])
        .stick_to_bottom(true)
        .show(ui, |ui| {
            for n in notes.iter().rev() {
                let colour = if n.local {
                    ui.visuals().warn_fg_color
                } else {
                    ui.visuals().weak_text_color()
                };
                ui.label(
                    RichText::new(format!("[{:7.2}s] {}", n.t_ms / 1000.0, n.text))
                        .small()
                        .color(colour),
                );
            }
        });
}

fn export_area(ui: &mut egui::Ui, app: &mut InfinityApp) {
    theme::section(ui, "Export");
    ui.add(
        egui::TextEdit::singleline(&mut app.export_text)
            .desired_width(f32::INFINITY)
            .hint_text("folder to write into"),
    );
    let has_data = app.model.genome_count() > 0;
    let mut action: Option<ExportAction> = None;
    ui.add_enabled_ui(has_data, |ui| {
        for (group, items) in Item::groups() {
            ui.label(
                RichText::new(*group)
                    .small()
                    .color(theme::palette(ui).text_weak),
            );
            ui.horizontal_wrapped(|ui| {
                for item in items.iter() {
                    if ui
                        .small_button(item.short_label())
                        .on_hover_text(item.label())
                        .clicked()
                    {
                        action = Some(ExportAction::Item(*item));
                    }
                }
            });
        }
        ui.add_space(2.0);
        let p = theme::palette(ui);
        if ui
            .add(theme::primary_button("Export everything", &p))
            .on_hover_text("every file above, into the export folder")
            .clicked()
        {
            action = Some(ExportAction::All);
        }
    });
    if !has_data {
        theme::hint(ui, "nothing to export yet");
    }
    if let Some(a) = action {
        let ctx = ui.ctx().clone();
        app.run_export(a, &ctx);
    }
}
