//! The Settings window: every run parameter, grouped by what a user is trying to
//! decide, with the common choices first and the algorithm's knobs folded away.
//!
//! * **Run** — where results go and how many cores to use.
//! * **Resolution** — how far refinement goes, the cluster distance (also shown
//!   as % identity), whether the distances between clusters are sharpened too,
//!   and the SNP cut-offs for outbreak-style SNP clusters.
//! * **Accuracy (SNP rules)** — which SNPs are trusted, with an on/off switch
//!   and an "advanced" fold for the numbers.
//! * **Advanced (algorithm)** — section 6.1's parameters and the experiment
//!   switches of section 11, folded by default.
//!
//! Every label carries a tooltip in plain words. Everything here applies at the
//! next run, except the cluster distance, which acts immediately while a run is
//! going (section 8) and is also on the toolbar.

use eframe::egui;
use egui::RichText;

use crate::app::InfinityApp;
use crate::settings::{self, Resolution, ResolutionSettings};
use crate::theme;

/// Show the window when it is open. Changes are validated and saved as they
/// are made, so there is no "Apply" button to forget.
pub fn show(ctx: &egui::Context, app: &mut InfinityApp) {
    if !app.show_settings {
        return;
    }
    let mut open = true;
    let mut changed = false;
    let mut dist_changed = false;
    let mut close = false;
    egui::Window::new(RichText::new("Settings").strong())
        .id(egui::Id::unique("settings-window"))
        .open(&mut open)
        .collapsible(false)
        .resizable(true)
        .default_size([460.0, 600.0])
        .min_width(380.0)
        .max_height(ctx.content_rect().height() - 40.0)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            egui::ScrollArea::vertical()
                .auto_shrink([false, true])
                .show(ui, |ui| {
                    theme::hint(
                        ui,
                        "Applied at the next run. The cluster distance also acts at once \
                         while a run is going.",
                    );
                    run_section(ui, app, &mut changed);
                    resolution_section(ui, app, &mut changed, &mut dist_changed);
                    accuracy_section(ui, app, &mut changed);
                    advanced_section(ui, app, &mut changed);
                });
            ui.add_space(4.0);
            theme::hairline(ui);
            if let Some(e) = &app.params_error {
                ui.colored_label(ui.visuals().error_fg_color, e);
            }
            ui.horizontal(|ui| {
                if ui
                    .button("Restore defaults")
                    .on_hover_text("every setting back to its default; folders are kept")
                    .clicked()
                {
                    let threads = app.prefs.params.threads;
                    app.prefs.params = infinity_core::params::Parameters::default();
                    app.prefs.params.threads = threads;
                    app.prefs.resolution = ResolutionSettings::default();
                    changed = true;
                    dist_changed = true;
                }
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let p = theme::palette(ui);
                    if ui.add(theme::primary_button("Done", &p)).clicked() {
                        close = true;
                    }
                });
            });
        });
    if !open || close {
        app.show_settings = false;
    }
    if changed || dist_changed {
        app.prefs.resolution.normalise();
        app.validate_params();
        app.save_prefs();
    }
    if dist_changed {
        app.push_cluster_dist();
    }
}

fn grid(ui: &mut egui::Ui, id: &str, add: impl FnOnce(&mut egui::Ui)) {
    egui::Grid::new(id)
        .num_columns(2)
        .spacing([14.0, 6.0])
        .min_col_width(130.0)
        .striped(false)
        .show(ui, add);
}

fn run_section(ui: &mut egui::Ui, app: &mut InfinityApp, changed: &mut bool) {
    theme::section(ui, "Run");
    grid(ui, "settings-run", |ui| {
        theme::grid_row(
            ui,
            "Output folder",
            "Where a run writes its results: trees, clusters, distance matrices, variants, \
             and events.jsonl, the recording that File → Open replays.",
            |ui| {
                *changed |= ui
                    .add(egui::TextEdit::singleline(&mut app.outdir_text).desired_width(230.0))
                    .changed();
            },
        );
        theme::grid_row(
            ui,
            "Threads",
            "How many processor cores to use. 0 uses all of them.",
            |ui| {
                *changed |= ui
                    .add(egui::DragValue::new(&mut app.prefs.params.threads).range(0..=512))
                    .changed();
                if app.prefs.params.threads == 0 {
                    theme::hint(ui, "all cores");
                }
            },
        );
    });
}

fn resolution_section(
    ui: &mut egui::Ui,
    app: &mut InfinityApp,
    changed: &mut bool,
    dist_changed: &mut bool,
) {
    theme::section(ui, "Resolution");
    let res = &mut app.prefs.resolution;
    ui.label("How far to refine")
        .on_hover_text("Every run starts blurry and sharpens; this is where it stops.");
    ui.horizontal_wrapped(|ui| {
        for r in Resolution::ALL {
            *changed |= ui
                .radio_value(&mut res.resolution, r, r.label())
                .on_hover_text(r.explain())
                .changed();
        }
    });
    theme::hint(
        ui,
        &format!(
            "{} The run stops at Focus {}.",
            res.resolution.explain(),
            res.resolution.focus()
        ),
    );
    ui.add_space(4.0);

    grid(ui, "settings-resolution", |ui| {
        theme::grid_row(
            ui,
            "Cluster distance",
            "How different two genomes may be and still be compared in detail, as a \
             sketch distance. 0.005 means about 99.5 % identical. Smaller makes more, \
             tighter clusters. This one acts immediately while a run is going.",
            |ui| {
                let before = app.prefs.params.cluster_dist;
                ui.add(
                    egui::Slider::new(&mut app.prefs.params.cluster_dist, 0.000_1..=0.05)
                        .logarithmic(true)
                        .fixed_decimals(4),
                );
                if (app.prefs.params.cluster_dist - before).abs() > f64::EPSILON {
                    *dist_changed = true;
                }
            },
        );
        theme::grid_row(
            ui,
            "",
            "The same threshold as a percent identity (1 − distance).",
            |ui| {
                let mut pct = settings::identity_pct(app.prefs.params.cluster_dist);
                let r = ui.add(
                    egui::DragValue::new(&mut pct)
                        .range(95.0..=99.999)
                        .speed(0.01)
                        .fixed_decimals(3)
                        .suffix(" % identity"),
                );
                if r.changed() {
                    app.prefs.params.cluster_dist = settings::dist_from_identity(pct);
                    *dist_changed = true;
                }
            },
        );
        theme::grid_row(
            ui,
            "Sharpen between clusters",
            "Also refine the distances between clusters, so the branches joining \
             clusters become SNP counts too, not only estimates from the sketches.",
            |ui| {
                *changed |= ui
                    .checkbox(&mut app.prefs.resolution.refine_backbone, "")
                    .changed();
            },
        );
    });

    ui.add_space(2.0);
    ui.label("SNP thresholds").on_hover_text(
        "Cut-offs for outbreak-style SNP clusters: genomes within this many SNPs of \
         each other are grouped, one grouping per threshold.",
    );
    let mut remove: Option<usize> = None;
    ui.horizontal_wrapped(|ui| {
        for (i, t) in app.prefs.resolution.snp_thresholds.iter().enumerate() {
            let chip = ui
                .button(format!("{t} SNPs  ✕"))
                .on_hover_text("remove this threshold");
            if chip.clicked() {
                remove = Some(i);
            }
        }
        let resp = ui.add(
            egui::TextEdit::singleline(&mut app.threshold_text)
                .desired_width(56.0)
                .hint_text("add"),
        );
        let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        if (ui.small_button("+").clicked() || enter) && !app.threshold_text.trim().is_empty() {
            match app.prefs.resolution.add_threshold(&app.threshold_text) {
                Ok(()) => {
                    app.threshold_text.clear();
                    app.threshold_error = None;
                    *changed = true;
                }
                Err(e) => app.threshold_error = Some(e),
            }
        }
    });
    if let Some(i) = remove {
        app.prefs.resolution.snp_thresholds.remove(i);
        *changed = true;
    }
    if let Some(e) = &app.threshold_error {
        ui.colored_label(ui.visuals().warn_fg_color, e);
    }
}

fn accuracy_section(ui: &mut egui::Ui, app: &mut InfinityApp, changed: &mut bool) {
    theme::section(ui, "Accuracy (SNP rules)");
    let f = &mut app.prefs.resolution.snp_filter;
    *changed |= ui
        .checkbox(&mut f.enabled, "Filter unreliable SNPs")
        .on_hover_text(
            "Leave out SNPs that are more likely assembly or repeat artefacts than real \
             differences: on short contigs, near contig ends, in high-copy contigs, or \
             in suspicious dense clusters.",
        )
        .changed();
    ui.add_enabled_ui(f.enabled, |ui| {
        *changed |= ui
            .checkbox(&mut f.core_only, "Core genome only")
            .on_hover_text(
                "Count only positions that every genome of the cluster covers, so a \
                 genome is not 'closer' just because part of it is missing.",
            )
            .changed();
        egui::CollapsingHeader::new("advanced")
            .id_salt("snp-rules-advanced")
            .default_open(false)
            .show(ui, |ui| {
                grid(ui, "settings-snp-rules", |ui| {
                    theme::grid_row(
                        ui,
                        "Min contig length",
                        "SNPs on contigs shorter than this are ignored: short contigs are \
                         often fragments of repeats or contamination.",
                        |ui| {
                            *changed |= ui
                                .add(
                                    egui::DragValue::new(&mut f.min_contig_len)
                                        .range(0..=1_000_000)
                                        .speed(100)
                                        .suffix(" bp"),
                                )
                                .changed();
                        },
                    );
                    theme::grid_row(
                        ui,
                        "Contig-end margin",
                        "SNPs this close to the end of a contig are ignored: assembly ends \
                         are where base calls are least reliable.",
                        |ui| {
                            *changed |= ui
                                .add(
                                    egui::DragValue::new(&mut f.contig_end_margin)
                                        .range(0..=100_000)
                                        .speed(5)
                                        .suffix(" bp"),
                                )
                                .changed();
                        },
                    );
                    theme::grid_row(
                        ui,
                        "High-copy contig ratio",
                        "A contig that looks this many times more abundant than the genome \
                         average is treated as a repeat (plasmid, phage, rRNA) and its SNPs \
                         are not counted.",
                        |ui| {
                            *changed |= ui
                                .add(
                                    egui::DragValue::new(&mut f.high_copy_ratio)
                                        .range(1.0..=100.0)
                                        .speed(0.05)
                                        .fixed_decimals(2)
                                        .suffix(" ×"),
                                )
                                .changed();
                        },
                    );
                    theme::grid_row(
                        ui,
                        "SNP-cluster window",
                        "Several SNPs packed into this window, more than chance allows, are \
                         treated as one recombination or mis-assembly event.",
                        |ui| {
                            *changed |= ui
                                .add(
                                    egui::DragValue::new(&mut f.cluster_window)
                                        .range(1..=1_000_000)
                                        .speed(10)
                                        .suffix(" bp"),
                                )
                                .changed();
                        },
                    );
                    theme::grid_row(
                        ui,
                        "Alpha",
                        "How surprising a dense SNP cluster must be before it is filtered: \
                         smaller filters less.",
                        |ui| {
                            *changed |= ui
                                .add(
                                    egui::DragValue::new(&mut f.alpha)
                                        .range(0.000_1..=0.5)
                                        .speed(0.001)
                                        .fixed_decimals(4),
                                )
                                .changed();
                        },
                    );
                });
            });
    });
    if let Err(e) = app.prefs.resolution.validate() {
        ui.colored_label(ui.visuals().error_fg_color, e);
    }
}

fn advanced_section(ui: &mut egui::Ui, app: &mut InfinityApp, changed: &mut bool) {
    ui.add_space(4.0);
    egui::CollapsingHeader::new(RichText::new("Advanced (algorithm)").strong())
        .id_salt("settings-advanced")
        .default_open(false)
        .show(ui, |ui| {
            theme::hint(
                ui,
                "The defaults suit bacterial genomes. Change these to experiment, not to \
                 get a better answer.",
            );
            let p = &mut app.prefs.params;
            grid(ui, "settings-advanced-grid", |ui| {
                theme::grid_row(
                    ui,
                    "k-mer length",
                    "Length of the words used as anchors. Must be odd and at most 31; \
                     longer is more specific, shorter tolerates more differences.",
                    |ui| {
                        *changed |= ui
                            .add(egui::DragValue::new(&mut p.k).range(1..=31))
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "Block densities",
                    "One anchor per this many k-mers at the coarse, middle and fine \
                     levels (about 16 kb, 2 kb and 256 bp blocks). Each must divide the \
                     previous one.",
                    |ui| {
                        for (i, d) in p.densities.iter_mut().enumerate() {
                            *changed |= ui
                                .add(egui::DragValue::new(d).range(1..=1_000_000).speed(16))
                                .on_hover_text(format!(
                                    "one anchor per this many k-mers at level {}",
                                    i + 1
                                ))
                                .changed();
                        }
                    },
                );
                theme::grid_row(
                    ui,
                    "Sketch level",
                    "Which level's anchors make up the sketch used for clustering (1 is \
                     coarsest).",
                    |ui| {
                        *changed |= ui
                            .add(egui::DragValue::new(&mut p.sketch_level).range(1..=3))
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "Max span ratio",
                    "How much longer one side of a matched stretch may be than the other \
                     before the anchors are judged out of order.",
                    |ui| {
                        *changed |= ui
                            .add(
                                egui::DragValue::new(&mut p.max_span_ratio)
                                    .range(1.0..=1000.0)
                                    .speed(0.1),
                            )
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "Max span slack",
                    "Extra bases allowed on top of the span ratio.",
                    |ui| {
                        *changed |= ui
                            .add(
                                egui::DragValue::new(&mut p.max_span_extra)
                                    .range(0..=10_000_000)
                                    .speed(1000)
                                    .suffix(" bp"),
                            )
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "Max aligned piece",
                    "Longest differing piece that is aligned base by base; longer pieces \
                     are reported as one complex difference.",
                    |ui| {
                        *changed |= ui
                            .add(
                                egui::DragValue::new(&mut p.dp_max_len)
                                    .range(1..=1_000_000)
                                    .suffix(" bp"),
                            )
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "Min identity",
                    "Below this identity an aligned piece is reported as one complex \
                     difference rather than as SNPs and indels.",
                    |ui| {
                        *changed |= ui
                            .add(
                                egui::DragValue::new(&mut p.min_identity)
                                    .range(0.01..=1.0)
                                    .speed(0.01)
                                    .fixed_decimals(2),
                            )
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "X-drop",
                    "How far the score may fall before extension at a contig end stops.",
                    |ui| {
                        *changed |= ui
                            .add(egui::DragValue::new(&mut p.xdrop).range(1..=1000))
                            .changed();
                    },
                );
                theme::grid_row(
                    ui,
                    "Memory budget",
                    "Member sequences are kept in memory while they fit in this budget. \
                     0 means a quarter of RAM.",
                    |ui| {
                        let mut mib = (p.memory_budget_bytes / (1 << 20)) as i64;
                        if ui
                            .add(
                                egui::DragValue::new(&mut mib)
                                    .range(0..=1_000_000)
                                    .suffix(" MiB"),
                            )
                            .changed()
                        {
                            p.memory_budget_bytes = (mib.max(0) as u64) << 20;
                            *changed = true;
                        }
                    },
                );
            });
            ui.add_space(4.0);
            ui.label(RichText::new("Experiment switches").strong())
                .on_hover_text(
                    "Each removes one part of the method while keeping the results \
                     identical; they exist to measure what each part is worth.",
                );
            for (flag, label, hint) in [
                (
                    &mut p.block_skip,
                    "Skip identical blocks",
                    "hash blocks and skip the ones the pivot already has",
                ),
                (
                    &mut p.hierarchical,
                    "Coarse-to-fine levels",
                    "use the nested levels top-down instead of fine blocks only",
                ),
                (
                    &mut p.piece_cache,
                    "Piece cache",
                    "resolve a difference shared by several genomes only once",
                ),
                (
                    &mut p.profile_grouping,
                    "Profile grouping",
                    "handle genomes with identical difference profiles once",
                ),
                (
                    &mut p.single_cluster,
                    "Single cluster",
                    "compare everything to one pivot instead of clustering first",
                ),
            ] {
                *changed |= ui.checkbox(flag, label).on_hover_text(hint).changed();
            }
        });
}
