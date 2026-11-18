//! The window's chrome: the toolbar along the top and the status bar along the
//! bottom (section 8).
//!
//! The toolbar carries the three things that must be reachable at any moment —
//! adding genomes, running and stopping — plus the cluster distance slider,
//! which is the one parameter section 8 requires to take effect immediately.
//! Everything else lives behind the File / Export / View menus so that the bar
//! stays one line high.

use eframe::egui;
use egui::{Context, RichText};

use crate::app::MtatoApp;
use crate::export::Item;
use crate::fmt;
use crate::prefs::MAX_RECENT;

/// The top bar: menus, transport controls, cluster distance and threads.
pub fn toolbar(ctx: &Context, app: &mut MtatoApp) {
    let mut show_inputs = app.show_inputs;
    let mut run_now = false;
    let mut stop_now = false;
    let mut open_now: Option<String> = None;
    let mut export_now: Option<Option<Item>> = None;
    let mut dist_changed = false;
    let mut params_changed = false;

    egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
        menu_bar(ui, app, &mut open_now, &mut export_now);

        ui.horizontal(|ui| {
            if ui
                .selectable_label(show_inputs, "Add genomes")
                .on_hover_text(
                    "Type or drop a path. No native file dialog is used, because the crates \
                     that provide one need system libraries this project must not require.",
                )
                .clicked()
            {
                show_inputs = !show_inputs;
            }
            let can_run = !app.running() && app.params_error.is_none() && !app.inputs.is_empty();
            ui.add_enabled_ui(can_run, |ui| {
                if ui
                    .button(RichText::new("▶ Run").strong())
                    .on_hover_text("Compare the listed genomes")
                    .clicked()
                {
                    run_now = true;
                }
            });
            ui.add_enabled_ui(app.running(), |ui| {
                if ui
                    .button("■ Stop")
                    .on_hover_text("Ask the engine to stop after the current genome")
                    .clicked()
                {
                    stop_now = true;
                }
            });
            if app.running() {
                ui.spinner();
            }
            ui.separator();

            // The one parameter that acts while a run is going.
            let before = app.prefs.params.cluster_dist;
            ui.label("cluster distance");
            ui.add(
                egui::Slider::new(&mut app.prefs.params.cluster_dist, 0.000_1..=0.05)
                    .logarithmic(true)
                    .show_value(true)
                    .fixed_decimals(4),
            );
            if (app.prefs.params.cluster_dist - before).abs() > f64::EPSILON {
                dist_changed = true;
            }

            ui.separator();
            ui.label("threads");
            if ui
                .add(egui::DragValue::new(&mut app.prefs.params.threads).range(0..=512))
                .changed()
            {
                params_changed = true;
            }
            if app.prefs.params.threads == 0 {
                ui.weak("(all cores)");
            }
            ui.separator();
            ui.label("output");
            let out = ui
                .add(egui::TextEdit::singleline(&mut app.outdir_text).desired_width(160.0))
                .on_hover_text("Where the run writes its files, including events.jsonl");
            if out.changed() {
                params_changed = true;
            }
        });

        if show_inputs {
            input_area(ui, app, &mut params_changed);
        }
    });

    app.show_inputs = show_inputs;
    if dist_changed {
        app.validate_params();
        app.push_cluster_dist();
    }
    if params_changed {
        app.validate_params();
        app.save_prefs();
    }
    if run_now {
        app.start_run(ctx);
    }
    if stop_now {
        app.stop_run();
    }
    if let Some(dir) = open_now {
        app.open_folder(&dir, ctx);
    }
    if let Some(item) = export_now {
        match item {
            Some(Item::Png) => app.export_png(ctx),
            other => app.export_files(other),
        }
    }
}

fn menu_bar(
    ui: &mut egui::Ui,
    app: &mut MtatoApp,
    open_now: &mut Option<String>,
    export_now: &mut Option<Option<Item>>,
) {
    egui::menu::MenuBar::default().ui(ui, |ui| {
        ui.menu_button("File", |ui| {
            ui.label("Open a finished output folder:");
            ui.horizontal(|ui| {
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut app.open_text)
                        .desired_width(240.0)
                        .hint_text("/path/to/mtato-out"),
                );
                let go = ui.button("Open").clicked()
                    || (resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)));
                if go {
                    *open_now = Some(app.open_text.clone());
                    ui.close();
                }
            });
            if !app.prefs.recent_inputs.is_empty() {
                ui.separator();
                ui.label("Recent:");
                let recent = app.prefs.recent_inputs.clone();
                for p in recent.iter().take(MAX_RECENT) {
                    let clicked = ui.button(p.as_str()).clicked();
                    if clicked {
                        app.open_text = p.clone();
                        *open_now = Some(p.clone());
                        ui.close();
                        break;
                    }
                }
            }
            ui.separator();
            if ui.button("Quit").clicked() {
                ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
            }
        });

        ui.menu_button("Export", |ui| {
            let has_data = app.model.genome_count() > 0;
            ui.add_enabled_ui(has_data, |ui| {
                for item in Item::files() {
                    if ui.button(item.label()).clicked() {
                        *export_now = Some(Some(item));
                        ui.close();
                    }
                }
                ui.separator();
                if ui.button("Tree picture (PNG)").clicked() {
                    *export_now = Some(Some(Item::Png));
                    ui.close();
                }
                if ui.button("Everything").clicked() {
                    *export_now = Some(None);
                    ui.close();
                }
            });
            if !has_data {
                ui.weak("nothing to export yet");
            }
            ui.separator();
            ui.label("Export folder:");
            ui.add(
                egui::TextEdit::singleline(&mut app.export_text)
                    .desired_width(240.0)
                    .hint_text("/path/to/export"),
            );
        });

        ui.menu_button("View", |ui| {
            ui.checkbox(&mut app.tree.show_legend, "Tree legend");
            ui.checkbox(&mut app.show_notes, "Warnings and notes");
            if ui.button("Fit tree").clicked() {
                app.tree.needs_fit = true;
                app.tree.user_moved = false;
            }
            if ui.button("Fit map").clicked() {
                app.map_view.needs_fit = true;
            }
        });
    });
}

/// The input area: a path field, the remembered folders, and the current list.
fn input_area(ui: &mut egui::Ui, app: &mut MtatoApp, params_changed: &mut bool) {
    ui.separator();
    ui.horizontal(|ui| {
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.input_text)
                .desired_width(320.0)
                .hint_text("a .fa / .fasta / .fna file, or a folder of them"),
        );
        let mut add = ui.button("Add").clicked();
        add |= resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        if add {
            let text = app.input_text.trim().to_string();
            if !text.is_empty() {
                let path = std::path::PathBuf::from(&text);
                if path.exists() {
                    if app.add_input_path(&path) {
                        app.status = format!("added {}", path.display());
                    } else {
                        app.status = format!("{} is already in the list", path.display());
                    }
                    app.input_text.clear();
                    *params_changed = true;
                } else {
                    app.status = format!("{} does not exist", path.display());
                }
            }
        }
        if ui.button("Clear list").clicked() {
            app.inputs.clear();
            *params_changed = true;
        }
        ui.weak("or drop files and folders anywhere on this window");
    });

    if !app.prefs.recent_inputs.is_empty() {
        let recent = app.prefs.recent_inputs.clone();
        ui.horizontal_wrapped(|ui| {
            ui.weak("recent:");
            for p in recent {
                let clicked = ui.selectable_label(false, p.as_str()).clicked();
                if clicked {
                    let path = std::path::PathBuf::from(&p);
                    if app.add_input_path(&path) {
                        app.status = format!("added {}", path.display());
                        *params_changed = true;
                    }
                    break;
                }
            }
        });
    }

    if !app.inputs.is_empty() {
        let mut remove: Option<usize> = None;
        egui::ScrollArea::vertical()
            .max_height(72.0)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                for (i, p) in app.inputs.iter().enumerate() {
                    ui.horizontal(|ui| {
                        if ui.small_button("✕").clicked() {
                            remove = Some(i);
                        }
                        ui.label(p.display().to_string());
                    });
                }
            });
        if let Some(i) = remove {
            app.remove_input(i);
            *params_changed = true;
        }
    }
}

/// The bottom bar: one segment per stage, with its progress and elapsed time.
pub fn status_bar(ctx: &Context, app: &MtatoApp) {
    egui::TopBottomPanel::bottom("statusbar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            for stage in &app.model.stages {
                let text = match stage.fraction() {
                    Some(f) => format!(
                        "{} {}/{} {}",
                        stage.name,
                        fmt::count(stage.done),
                        fmt::count(stage.total),
                        fmt::secs(stage.elapsed_ms())
                    ),
                    // A stage with no total has no fraction to show; only its
                    // clock is honest.
                    None => format!("{} {}", stage.name, fmt::secs(stage.elapsed_ms())),
                };
                let colour = if stage.finished {
                    ui.visuals().weak_text_color()
                } else {
                    ui.visuals().strong_text_color()
                };
                let resp = ui.colored_label(colour, text);
                if let Some(f) = stage.fraction() {
                    resp.on_hover_ui(|ui| {
                        ui.add(egui::ProgressBar::new(f as f32).show_percentage());
                    });
                }
                ui.separator();
            }
            if app.model.stages.is_empty() {
                ui.weak("no stage has started");
                ui.separator();
            }
            ui.label(RichText::new(&app.status).strong());
        });
        ui.horizontal(|ui| {
            ui.weak(format!(
                "{} genomes | {} clusters | focus 0/1/2 = {}/{}/{} | events {}",
                fmt::count(app.model.genome_count() as u64),
                fmt::count(app.model.cluster_count() as u64),
                app.model
                    .genomes
                    .iter()
                    .filter(|g| g.focus == fmt::FOCUS_SKETCH)
                    .count(),
                app.model
                    .genomes
                    .iter()
                    .filter(|g| g.focus == fmt::FOCUS_BLOCKS)
                    .count(),
                app.model
                    .genomes
                    .iter()
                    .filter(|g| g.focus == fmt::FOCUS_EXACT)
                    .count(),
                fmt::count(app.model.events_seen),
            ));
            match &app.model.source {
                crate::state::Source::Live => ui.weak("(live run)"),
                crate::state::Source::Folder(p) => {
                    ui.weak(format!("(opened {})", p.display()))
                }
                crate::state::Source::Empty => {}
            }
        });
    });
}
