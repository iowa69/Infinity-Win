//! The window's chrome: the toolbar along the top and the status bar along the
//! bottom (section 8).
//!
//! The toolbar is one line: the menus, the three things that must be reachable
//! at any moment (genomes, run, stop), the cluster distance (the one parameter
//! section 8 requires to act immediately, shown with its % identity), and on
//! the right the Settings window and the side-panel toggle. Every other
//! parameter lives in Settings. Below it, the "Genomes" drawer holds the input
//! paths and the output folder: everything a run needs, in one place.

use eframe::egui;
use egui::RichText;

use crate::app::InfinityApp;
use crate::export::Item;
use crate::prefs::MAX_RECENT;
use crate::{fmt, settings, theme};

/// The top bar: menus, transport controls, cluster distance, settings.
pub fn toolbar(ui: &mut egui::Ui, app: &mut InfinityApp) {
    let mut run_now = false;
    let mut stop_now = false;
    let mut open_now: Option<String> = None;
    let mut export_now: Option<ExportAction> = None;
    let mut dist_changed = false;
    let mut params_changed = false;
    let p = theme::palette(ui);

    egui::Panel::top("toolbar")
        .frame(
            egui::Frame::side_top_panel(ui.style())
                .inner_margin(egui::Margin::symmetric(10, 6))
                .fill(p.panel),
        )
        .show(ui, |ui| {
            egui::menu::MenuBar::default().ui(ui, |ui| {
                ui.label(
                    RichText::new("∞ Infinity-Win")
                        .strong()
                        .size(14.0)
                        .color(p.accent),
                );
                ui.add_space(6.0);
                menus(ui, app, &mut open_now, &mut export_now);
                ui.separator();

                add_buttons(ui, app);
                let n = app.inputs.len();
                let label = if n == 0 {
                    "Inputs ▾".to_string()
                } else {
                    format!("Inputs ({}) ▾", fmt::count(app.input_genome_total() as u64))
                };
                if ui
                    .selectable_label(app.show_inputs, label)
                    .on_hover_text(
                        "The genome files and folders of the next run, how many genomes each \
                         holds, and where the results go.",
                    )
                    .clicked()
                {
                    app.show_inputs = !app.show_inputs;
                }
                let can_run =
                    !app.running() && app.params_error.is_none() && !app.inputs.is_empty();
                let run = ui
                    .add_enabled(can_run, theme::primary_button("▶ Run", &p))
                    .on_hover_text("Compare the listed genomes")
                    .on_disabled_hover_text(if app.inputs.is_empty() {
                        "add genomes first"
                    } else if app.running() {
                        "a run is going"
                    } else {
                        "fix the settings first"
                    });
                if run.clicked() {
                    run_now = true;
                }
                if ui
                    .add_enabled(app.running(), egui::Button::new("■ Stop"))
                    .on_hover_text("Stop after the current genome")
                    .clicked()
                {
                    stop_now = true;
                }
                if app.running() {
                    ui.spinner();
                }
                ui.separator();

                // The one parameter that acts while a run is going.
                let before = app.prefs.params.cluster_dist;
                ui.label(RichText::new("Cluster distance").color(p.text_weak))
                    .on_hover_text(
                        "How different genomes may be and still share a cluster. Acts \
                         immediately while a run is going: clusters are recomputed from \
                         the sketches.",
                    );
                ui.add(
                    egui::Slider::new(&mut app.prefs.params.cluster_dist, 0.000_1..=0.05)
                        .logarithmic(true)
                        .show_value(true)
                        .fixed_decimals(4),
                );
                ui.label(
                    RichText::new(format!(
                        "{:.2} % id",
                        settings::identity_pct(app.prefs.params.cluster_dist)
                    ))
                    .small()
                    .color(p.text_weak),
                );
                if (app.prefs.params.cluster_dist - before).abs() > f64::EPSILON {
                    dist_changed = true;
                }

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui
                        .selectable_label(app.show_side, "◨")
                        .on_hover_text("show or hide the run panel on the right")
                        .clicked()
                    {
                        app.show_side = !app.show_side;
                    }
                    let warn =
                        app.params_error.is_some() || app.prefs.resolution.validate().is_err();
                    let gear = if warn {
                        "⚙ Settings ⚠"
                    } else {
                        "⚙ Settings"
                    };
                    if ui
                        .selectable_label(app.show_settings, gear)
                        .on_hover_text("resolution, SNP rules, threads and the algorithm")
                        .clicked()
                    {
                        app.show_settings = !app.show_settings;
                    }
                });
            });

            if app.show_inputs {
                input_area(ui, app, &mut params_changed);
            }
        });

    let ctx = ui.ctx().clone();
    if dist_changed {
        app.validate_params();
        app.push_cluster_dist();
        app.save_prefs();
    }
    if params_changed {
        app.validate_params();
        app.save_prefs();
    }
    if run_now {
        app.start_run(&ctx);
    }
    if stop_now {
        app.stop_run();
    }
    if let Some(dir) = open_now {
        app.open_folder(&dir, &ctx);
    }
    if let Some(a) = export_now {
        app.run_export(a, &ctx);
    }
}

/// The "Add files…" and "Add folder…" buttons: native dialogs, opened off the
/// UI thread (see [`crate::dialogs`]).
pub(crate) fn add_buttons(ui: &mut egui::Ui, app: &mut InfinityApp) {
    use crate::dialogs::PickKind;
    let busy = app.pending_pick.is_some();
    if ui
        .add_enabled(!busy, egui::Button::new("＋ Add files…"))
        .on_hover_text("choose FASTA files (.fasta, .fa, .fna, .fas, also .gz)")
        .on_disabled_hover_text("a file dialog is open")
        .clicked()
    {
        app.request_pick(PickKind::Files);
    }
    if ui
        .add_enabled(!busy, egui::Button::new("＋ Add folder…"))
        .on_hover_text("choose folders: every FASTA file directly inside each is a genome")
        .on_disabled_hover_text("a file dialog is open")
        .clicked()
    {
        app.request_pick(PickKind::Folders);
    }
}

/// What an export menu entry asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportAction {
    /// One export item.
    Item(Item),
    /// Everything the window can export.
    All,
}

fn menus(
    ui: &mut egui::Ui,
    app: &mut InfinityApp,
    open_now: &mut Option<String>,
    export_now: &mut Option<ExportAction>,
) {
    ui.menu_button("File", |ui| {
        if ui.button("Add genome files…").clicked() {
            app.request_pick(crate::dialogs::PickKind::Files);
            ui.close();
        }
        if ui.button("Add a folder of genomes…").clicked() {
            app.request_pick(crate::dialogs::PickKind::Folders);
            ui.close();
        }
        ui.separator();
        ui.label("Open a finished output folder:");
        ui.horizontal(|ui| {
            let resp = ui.add(
                egui::TextEdit::singleline(&mut app.open_text)
                    .desired_width(240.0)
                    .hint_text("/path/to/infinity-out"),
            );
            let go = ui.button("Open").clicked()
                || (resp.lost_focus() && ui.ctx().input(|i| i.key_pressed(egui::Key::Enter)));
            if go {
                *open_now = Some(app.open_text.clone());
                ui.close();
            }
        });
        if !app.prefs.recent_inputs.is_empty() {
            ui.separator();
            ui.label(RichText::new("Recent").small());
            let recent = app.prefs.recent_inputs.clone();
            for p in recent.iter().take(MAX_RECENT) {
                if ui.button(p.as_str()).clicked() {
                    app.open_text = p.clone();
                    *open_now = Some(p.clone());
                    ui.close();
                    break;
                }
            }
        }
        ui.separator();
        if ui.button("Settings…").clicked() {
            app.show_settings = true;
            ui.close();
        }
        if ui.button("Quit").clicked() {
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
        }
    });

    ui.menu_button("Export", |ui| {
        let has_data = app.model.genome_count() > 0;
        ui.add_enabled_ui(has_data, |ui| {
            for item in Item::menu() {
                if ui.button(item.label()).clicked() {
                    *export_now = Some(ExportAction::Item(item));
                    ui.close();
                }
            }
            ui.separator();
            if ui.button("Everything").clicked() {
                *export_now = Some(ExportAction::All);
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
        ui.checkbox(&mut app.show_side, "Run panel");
        ui.checkbox(&mut app.show_notes, "Warnings and notes");
        ui.checkbox(&mut app.tree.opts.legend, "Tree legend");
        ui.separator();
        ui.label(RichText::new("Theme").small());
        let mut theme_name = if app.prefs.theme.is_empty() {
            "dark".to_string()
        } else {
            app.prefs.theme.clone()
        };
        let before = theme_name.clone();
        ui.horizontal(|ui| {
            ui.selectable_value(&mut theme_name, "dark".to_string(), "Dark");
            ui.selectable_value(&mut theme_name, "light".to_string(), "Light");
            ui.selectable_value(&mut theme_name, "system".to_string(), "System");
        });
        if theme_name != before {
            app.prefs.theme = theme_name.clone();
            ui.ctx().set_theme(theme::preference(&theme_name));
            app.save_prefs();
        }
        ui.separator();
        if ui.button("Fit tree").clicked() {
            app.tree.needs_fit = true;
            app.tree.user_moved = false;
        }
        if ui.button("Fit map").clicked() {
            app.map_view.needs_fit = true;
        }
    });
}

/// The genomes drawer: a path field, the remembered folders, the current list,
/// and the output folder.
fn input_area(ui: &mut egui::Ui, app: &mut InfinityApp, params_changed: &mut bool) {
    let p = theme::palette(ui);
    ui.add_space(2.0);
    theme::hairline(ui);
    ui.horizontal_wrapped(|ui| {
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.input_text)
                .desired_width(260.0)
                .hint_text("or type a file or folder path"),
        );
        let mut add = ui.button("Add path").clicked();
        add |= resp.lost_focus() && ui.ctx().input(|i| i.key_pressed(egui::Key::Enter));
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
        ui.label(
            RichText::new("or drop files and folders on the window")
                .small()
                .color(p.text_weak),
        );
        ui.separator();
        ui.label(RichText::new("Results to").color(p.text_weak))
            .on_hover_text("Where the run writes its files, including events.jsonl");
        if ui
            .add(egui::TextEdit::singleline(&mut app.outdir_text).desired_width(200.0))
            .changed()
        {
            *params_changed = true;
        }
    });

    if !app.prefs.recent_inputs.is_empty() {
        let recent = app.prefs.recent_inputs.clone();
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new("recent").small().color(p.text_weak));
            for r in recent {
                let clicked = ui
                    .add(egui::Button::new(RichText::new(r.as_str()).small()).small())
                    .clicked();
                if clicked {
                    let path = std::path::PathBuf::from(&r);
                    if app.add_input_path(&path) {
                        app.status = format!("added {}", path.display());
                        *params_changed = true;
                    }
                    break;
                }
            }
        });
    }

    if app.inputs.is_empty() {
        ui.horizontal_wrapped(|ui| {
            ui.label(
                RichText::new("No genomes in the list yet.")
                    .small()
                    .color(p.text_weak),
            );
            add_buttons(ui, app);
        });
    } else {
        let mut remove: Option<usize> = None;
        let mut clear = false;
        // One line per file or folder, with what it holds. Shrinks to its rows:
        // a fixed height left a band of empty space when one path was listed.
        let entries: Vec<(std::path::PathBuf, Option<usize>)> = app
            .inputs
            .clone()
            .into_iter()
            .map(|path| {
                let n = app.input_count(&path);
                (path, n)
            })
            .collect();
        egui::ScrollArea::vertical()
            .max_height(96.0)
            .auto_shrink([false, true])
            .show(ui, |ui| {
                for (i, (path, n)) in entries.iter().enumerate() {
                    ui.horizontal(|ui| {
                        let what = if path.is_dir() { "folder" } else { "file" };
                        if ui
                            .small_button("✕")
                            .on_hover_text(format!("remove this {what} from the list"))
                            .clicked()
                        {
                            remove = Some(i);
                        }
                        let icon = if path.is_dir() { "🗀" } else { "🗋" };
                        ui.label(RichText::new(icon).color(p.text_weak));
                        ui.label(RichText::new(path.display().to_string()).color(p.text));
                        let colour = match n {
                            Some(0) | None => p.warn,
                            _ => p.text_weak,
                        };
                        ui.label(
                            RichText::new(crate::dialogs::count_text(*n))
                                .small()
                                .color(colour),
                        );
                    });
                }
            });
        ui.horizontal(|ui| {
            ui.label(
                RichText::new(format!(
                    "{} item(s) · {} in total",
                    entries.len(),
                    crate::dialogs::count_text(Some(app.input_genome_total()))
                ))
                .small()
                .color(p.text_weak),
            );
            if ui
                .small_button("Clear all")
                .on_hover_text("empty the list")
                .clicked()
            {
                clear = true;
            }
        });
        if let Some(i) = remove {
            app.remove_input(i);
            *params_changed = true;
        }
        if clear {
            app.inputs.clear();
            app.input_counts.clear();
            *params_changed = true;
        }
    }
}

/// The bottom bar: one chip per stage, with its progress and elapsed time, then
/// the status line and the totals.
pub fn status_bar(ui: &mut egui::Ui, app: &InfinityApp) {
    let p = theme::palette(ui);
    egui::Panel::bottom("statusbar")
        .frame(
            egui::Frame::side_top_panel(ui.style())
                .inner_margin(egui::Margin::symmetric(10, 4))
                .fill(p.panel),
        )
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                for stage in &app.model.stages {
                    stage_chip(ui, stage, &p);
                }
                if app.model.stages.is_empty() {
                    ui.label(RichText::new("idle").small().color(p.text_weak));
                }
                ui.add_space(6.0);
                ui.label(RichText::new(&app.status).color(p.text_strong));
            });
            ui.horizontal_wrapped(|ui| {
                let focus_count = |f: u8| app.model.genomes.iter().filter(|g| g.focus == f).count();
                ui.label(
                    RichText::new(format!(
                        "{} genomes ({}) · {} clusters · focus 0/1/2 = {}/{}/{} · {} events",
                        fmt::count(app.model.genome_count() as u64),
                        fmt::bp(app.model.total_length()),
                        fmt::count(app.model.cluster_count() as u64),
                        focus_count(fmt::FOCUS_SKETCH),
                        focus_count(fmt::FOCUS_BLOCKS),
                        focus_count(fmt::FOCUS_EXACT),
                        fmt::count(app.model.events_seen),
                    ))
                    .small()
                    .color(p.text_weak),
                );
                let source = match &app.model.source {
                    crate::state::Source::Live => "live run".to_string(),
                    crate::state::Source::Folder(dir) => format!("opened {}", dir.display()),
                    crate::state::Source::Empty => String::new(),
                };
                if !source.is_empty() {
                    ui.label(RichText::new(source).small().color(p.text_weak));
                }
            });
        });
}

/// One stage of the status bar: name, count and clock, with a thin progress
/// line under it.
fn stage_chip(ui: &mut egui::Ui, stage: &crate::state::StageState, p: &theme::Palette) {
    let text = match stage.fraction() {
        Some(_) => format!(
            "{} {}/{} · {}",
            stage.name,
            fmt::count(stage.done),
            fmt::count(stage.total),
            fmt::secs(stage.elapsed_ms())
        ),
        // A stage with no total has no fraction to show; only its clock is
        // honest.
        None => format!("{} · {}", stage.name, fmt::secs(stage.elapsed_ms())),
    };
    let colour = if stage.finished {
        p.text_weak
    } else {
        p.text_strong
    };
    let frame = egui::Frame::new()
        .fill(p.surface)
        .corner_radius(10)
        .inner_margin(egui::Margin::symmetric(8, 2));
    let resp = frame
        .show(ui, |ui| {
            ui.label(RichText::new(text).small().color(colour));
        })
        .response;
    if let Some(f) = stage.fraction() {
        let r = resp.rect;
        let y = r.bottom() - 1.5;
        let x1 = r.left() + 6.0 + (r.width() - 12.0) * f as f32;
        ui.painter().line_segment(
            [egui::pos2(r.left() + 6.0, y), egui::pos2(x1, y)],
            egui::Stroke::new(
                1.5,
                if stage.finished {
                    p.text_weak
                } else {
                    p.accent
                },
            ),
        );
        resp.on_hover_ui(|ui| {
            ui.add(egui::ProgressBar::new(f as f32).show_percentage());
        });
    }
}
