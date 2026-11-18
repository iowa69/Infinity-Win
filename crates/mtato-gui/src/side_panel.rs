//! The right-hand panel: every parameter of section 6.1, the live counters of
//! section 7, the warnings, and the export folder.
//!
//! Only the cluster distance acts while a run is going, and that slider lives in
//! the toolbar where it cannot be missed. Everything here applies at the next
//! run, which is stated on the panel rather than left for the user to discover.

use eframe::egui;
use egui::{Context, RichText};

use crate::app::MtatoApp;
use crate::fmt;
use crate::state::Note;

pub fn show(ctx: &Context, app: &mut MtatoApp) {
    egui::SidePanel::right("parameters")
        .default_width(290.0)
        .min_width(230.0)
        .resizable(true)
        .show(ctx, |ui| {
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    parameters(ui, app);
                    ui.separator();
                    counters(ui, app);
                    ui.separator();
                    summary(ui, app);
                    ui.separator();
                    notes(ui, app);
                    ui.separator();
                    export_area(ui, app);
                });
        });
}

fn parameters(ui: &mut egui::Ui, app: &mut MtatoApp) {
    ui.label(RichText::new("PARAMETERS").strong());
    ui.weak("applied at the next run, except the cluster distance");
    let mut changed = false;

    ui.horizontal(|ui| {
        ui.label("k");
        changed |= ui
            .add(egui::DragValue::new(&mut app.prefs.params.k).range(1..=31))
            .changed();
        ui.weak("odd, ≤ 31");
    });

    ui.label("densities (levels 1/2/3)");
    ui.horizontal(|ui| {
        for (i, d) in app.prefs.params.densities.iter_mut().enumerate() {
            changed |= ui
                .add(egui::DragValue::new(d).range(1..=1_000_000).speed(16))
                .on_hover_text(format!("one anchor per this many k-mers at level {}", i + 1))
                .changed();
        }
    });

    ui.horizontal(|ui| {
        ui.label("sketch level");
        changed |= ui
            .add(egui::DragValue::new(&mut app.prefs.params.sketch_level).range(1..=3))
            .changed();
    });

    ui.horizontal(|ui| {
        ui.label("cluster dist");
        changed |= ui
            .add(
                egui::DragValue::new(&mut app.prefs.params.cluster_dist)
                    .range(0.000_01..=0.99)
                    .speed(0.0001)
                    .fixed_decimals(5),
            )
            .changed();
    });

    ui.horizontal(|ui| {
        ui.label("max span ratio");
        changed |= ui
            .add(
                egui::DragValue::new(&mut app.prefs.params.max_span_ratio)
                    .range(1.0..=1000.0)
                    .speed(0.1),
            )
            .changed();
    });
    ui.horizontal(|ui| {
        ui.label("max span extra");
        changed |= ui
            .add(
                egui::DragValue::new(&mut app.prefs.params.max_span_extra)
                    .range(0..=10_000_000)
                    .speed(1000),
            )
            .changed();
    });
    ui.horizontal(|ui| {
        ui.label("dp max len");
        changed |= ui
            .add(egui::DragValue::new(&mut app.prefs.params.dp_max_len).range(1..=1_000_000))
            .changed();
    });
    ui.horizontal(|ui| {
        ui.label("min identity");
        changed |= ui
            .add(
                egui::DragValue::new(&mut app.prefs.params.min_identity)
                    .range(0.01..=1.0)
                    .speed(0.01)
                    .fixed_decimals(2),
            )
            .changed();
    });
    ui.horizontal(|ui| {
        ui.label("xdrop");
        changed |= ui
            .add(egui::DragValue::new(&mut app.prefs.params.xdrop).range(1..=1000))
            .changed();
    });
    ui.horizontal(|ui| {
        ui.label("threads");
        changed |= ui
            .add(egui::DragValue::new(&mut app.prefs.params.threads).range(0..=512))
            .changed();
        if app.prefs.params.threads == 0 {
            ui.weak("all cores");
        }
    });

    ui.add_space(4.0);
    ui.label(RichText::new("Experiment switches").underline());
    for (flag, label, hint) in [
        (
            &mut app.prefs.params.block_skip,
            "block skip",
            "hash blocks and skip the identical ones (6.6 A)",
        ),
        (
            &mut app.prefs.params.hierarchical,
            "hierarchical",
            "use the nested levels top-down instead of leaf blocks only (6.10)",
        ),
        (
            &mut app.prefs.params.piece_cache,
            "piece cache",
            "resolve an identical difference piece once (6.10)",
        ),
        (
            &mut app.prefs.params.profile_grouping,
            "profile grouping",
            "handle members with identical profiles once (6.9)",
        ),
        (
            &mut app.prefs.params.single_cluster,
            "single cluster",
            "compare everything to one pivot instead of clustering first",
        ),
    ] {
        changed |= ui.checkbox(flag, label).on_hover_text(hint).changed();
    }
    ui.horizontal(|ui| {
        ui.label("memory budget");
        let mut mib = (app.prefs.params.memory_budget_bytes / (1 << 20)) as i64;
        if ui
            .add(egui::DragValue::new(&mut mib).range(0..=1_000_000).suffix(" MiB"))
            .changed()
        {
            app.prefs.params.memory_budget_bytes = (mib.max(0) as u64) << 20;
            changed = true;
        }
        ui.weak(if app.prefs.params.memory_budget_bytes == 0 {
            "a quarter of RAM"
        } else {
            ""
        });
    });

    if changed {
        app.validate_params();
        app.save_prefs();
    }
    if let Some(e) = &app.params_error {
        ui.colored_label(ui.visuals().error_fg_color, e);
    }
}

fn counters(ui: &mut egui::Ui, app: &MtatoApp) {
    ui.label(RichText::new("LIVE COUNTERS").strong());
    let c = &app.model.counters;
    // "blocks identical 99.7 %" is the number the whole idea is judged on, so it
    // goes first and gets a bar.
    let identity = c.block_identity();
    ui.horizontal(|ui| {
        ui.label("blocks identical");
        ui.strong(fmt::pct(identity));
    });
    if let Some(f) = identity {
        ui.add(egui::ProgressBar::new(f as f32).show_percentage());
    }
    row(ui, "blocks looked up", fmt::count(c.blocks_looked_up));
    row(
        ui,
        "matched L1/L2/L3",
        format!(
            "{} / {} / {}",
            fmt::count(c.blocks_matched[0]),
            fmt::count(c.blocks_matched[1]),
            fmt::count(c.blocks_matched[2])
        ),
    );
    row(ui, "pieces compared", fmt::count(c.pieces_compared));
    row(ui, "cache hits", fmt::count(c.cache_hits));
    row(ui, "alignments run", fmt::count(c.alignments_run));
    row(ui, "bases aligned", fmt::count(c.bases_aligned));
    row(ui, "k-mers scanned", fmt::count(c.kmers_scanned));
    row(ui, "genomes scanned", fmt::count(c.genomes_scanned));
    row(ui, "end extensions", fmt::count(c.end_extensions));
    row(ui, "stray points dropped", fmt::count(c.stray_points_dropped));
    row(ui, "breakpoints", fmt::count(c.breakpoints));
}

fn row(ui: &mut egui::Ui, label: &str, value: String) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).weak());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(value);
        });
    });
}

fn summary(ui: &mut egui::Ui, app: &MtatoApp) {
    let Some(s) = &app.model.summary else {
        return;
    };
    ui.label(RichText::new(if s.cancelled {
        "STOPPED"
    } else {
        "FINISHED"
    })
    .strong());
    row(ui, "genomes", fmt::count(u64::from(s.genomes)));
    row(ui, "clusters", fmt::count(u64::from(s.clusters)));
    row(
        ui,
        "focus 0/1/2",
        format!("{}/{}/{}", s.focus_counts[0], s.focus_counts[1], s.focus_counts[2]),
    );
    row(
        ui,
        "first tree",
        s.timings
            .first_tree_ms
            .map(fmt::secs)
            .unwrap_or_else(|| "not reached".to_string()),
    );
    row(
        ui,
        "exact tree",
        s.timings
            .exact_tree_ms
            .map(fmt::secs)
            .unwrap_or_else(|| "not reached".to_string()),
    );
    if let Some(t) = s.timings.total_ms {
        row(ui, "total", fmt::secs(t));
    }
    match s.peak_memory_bytes {
        // An unknown peak is reported as unknown, never as zero.
        Some(b) => row(ui, "peak RSS", fmt::bytes(b)),
        None => row(ui, "peak RSS", "unknown".to_string()),
    }
    for n in &s.notes {
        ui.weak(n);
    }
}

fn notes(ui: &mut egui::Ui, app: &mut MtatoApp) {
    let count = app.model.notes.len();
    let header = RichText::new(format!("WARNINGS AND NOTES ({count})")).strong();
    ui.checkbox(&mut app.show_notes, header);
    if !app.show_notes || count == 0 {
        return;
    }
    let notes: Vec<Note> = app.model.notes.clone();
    egui::ScrollArea::vertical()
        .max_height(160.0)
        .auto_shrink([false, false])
        .stick_to_bottom(true)
        .show(ui, |ui| {
            for n in notes.iter().rev() {
                let colour = if n.local {
                    ui.visuals().warn_fg_color
                } else {
                    ui.visuals().weak_text_color()
                };
                ui.colored_label(colour, format!("[{:7.2}s] {}", n.t_ms / 1000.0, n.text));
            }
        });
}

fn export_area(ui: &mut egui::Ui, app: &mut MtatoApp) {
    ui.label(RichText::new("EXPORT").strong());
    ui.add(
        egui::TextEdit::singleline(&mut app.export_text)
            .desired_width(f32::INFINITY)
            .hint_text("folder to write into"),
    );
    let has_data = app.model.genome_count() > 0;
    ui.add_enabled_ui(has_data, |ui| {
        ui.horizontal_wrapped(|ui| {
            for item in crate::export::Item::files() {
                if ui.button(short_label(item)).clicked() {
                    app.export_files(Some(item));
                }
            }
            if ui.button("PNG").clicked() {
                app.export_png(ui.ctx());
            }
            if ui.button("all").clicked() {
                app.export_files(None);
            }
        });
    });
    if !has_data {
        ui.weak("nothing to export yet");
    }
    if let Some(cid) = app.selected_cluster {
        ui.weak(format!("cluster {cid} is the one the map shows"));
    }
}

/// A button-sized name for an export item.
fn short_label(item: crate::export::Item) -> &'static str {
    match item {
        crate::export::Item::Tree => "tree.nwk",
        crate::export::Item::Clusters => "clusters",
        crate::export::Item::Genomes => "genomes",
        crate::export::Item::Distances => "distances",
        crate::export::Item::Variants => "variants",
        crate::export::Item::Png => "PNG",
    }
}
