//! The divergent-assembly warning: a modal dialog that lists the genomes the
//! engine flagged (another species, contamination, an incomplete assembly) and
//! offers to leave them out and start the same run again at once.
//!
//! The engine reports them right after the scan, before any refinement, so the
//! choice costs seconds. "Continue anyway" keeps the run going and leaves a
//! warning in the Run panel for as long as these results are on screen.

use eframe::egui;
use egui::RichText;

use infinity_core::events::OutlierInfo;
use infinity_core::outlier::OutlierThresholds;

use crate::app::InfinityApp;
use crate::state::UiModel;
use crate::{fmt, theme};

/// The dialog's own state: which flagged genomes are ticked.
#[derive(Debug, Clone, PartialEq)]
pub struct OutlierDialog {
    pub list: Vec<OutlierInfo>,
    pub ticked: Vec<bool>,
}

impl OutlierDialog {
    pub fn new(list: Vec<OutlierInfo>) -> Self {
        let ticked = vec![true; list.len()];
        Self { list, ticked }
    }

    /// Names of the ticked genomes, in the engine's order.
    pub fn ticked_names(&self) -> Vec<String> {
        self.list
            .iter()
            .zip(&self.ticked)
            .filter(|(_, &t)| t)
            .map(|(o, _)| o.name.clone())
            .collect()
    }
}

/// Add names to an exclusion list, keeping its order and skipping repeats.
/// Returns how many were new.
pub fn add_exclusions(list: &mut Vec<String>, names: &[String]) -> usize {
    let mut added = 0;
    for n in names {
        let n = n.trim();
        if !n.is_empty() && !list.iter().any(|x| x == n) {
            list.push(n.to_string());
            added += 1;
        }
    }
    added
}

/// Short tags for the table: why, in a few words. The full sentences are in
/// the tooltip and in `outliers.tsv`.
pub fn tags(o: &OutlierInfo, model: &UiModel) -> Vec<String> {
    let t = OutlierThresholds::default();
    let mut out = Vec::new();
    if o.nearest_dist >= 1.0 {
        out.push("shares nothing with any other genome".to_string());
    } else if o.nearest_dist > t.max_nearest_dist {
        out.push(format!(
            "no close relative (≈{:.0}% identity)",
            (1.0 - o.nearest_dist) * 100.0
        ));
    }
    let med_len = model.median_length() as f64;
    if med_len > 0.0 {
        let r = o.length as f64 / med_len;
        if r > t.max_length_ratio || r < t.min_length_ratio {
            out.push(format!("length {r:.2}× median"));
        }
    }
    let gcs: Vec<f64> = model.genomes.iter().map(|g| g.gc).collect();
    if let Some(m) = infinity_core::outlier::median(&gcs) {
        let d = (o.gc - m) * 100.0;
        if d.abs() > t.max_gc_diff * 100.0 {
            out.push(format!("GC {d:+.1} points"));
        }
    }
    if out.is_empty() {
        // The model has not seen every genome (an opened folder written by
        // an older version, say): fall back to the engine's own words.
        out = o
            .reasons
            .iter()
            .map(|r| r.split(':').next().unwrap_or(r).to_string())
            .collect();
    }
    out
}

/// What the user chose this frame.
#[derive(Debug, Clone, PartialEq)]
enum Choice {
    Restart(Vec<String>),
    AddOnly(Vec<String>),
    Continue,
}

/// Show the dialog while it is open.
pub fn show(ctx: &egui::Context, app: &mut InfinityApp) {
    let Some(dialog) = app.outlier_dialog.as_mut() else {
        return;
    };
    let model = &app.model;
    // A finished folder that was opened has no genomes to re-run with.
    let can_restart = !app.inputs.is_empty() && model.source == crate::state::Source::Live;
    let total = model.genome_count();
    let mut choice: Option<Choice> = None;

    egui::Modal::new(egui::Id::unique("outlier-dialog")).show(ctx, |ui| {
        let p = theme::palette(ui);
        ui.set_max_width((ctx.content_rect().width() - 80.0).clamp(360.0, 860.0));
        ui.horizontal(|ui| {
            ui.label(RichText::new("⚠").size(22.0).color(p.warn));
            ui.label(
                RichText::new("Divergent assemblies")
                    .size(17.0)
                    .strong()
                    .color(p.text_strong),
            );
        });
        ui.add_space(2.0);
        ui.label(format!(
            "{} of {} genomes do not look like the rest: probably another species, a \
             mislabelled genome, a contaminated or an incomplete assembly. Compared with the \
             others they can distort the clusters, shrink the core genome and change the \
             SNP counts.",
            dialog.list.len(),
            total
        ));
        ui.add_space(6.0);

        egui::ScrollArea::vertical()
            .max_height((ctx.content_rect().height() * 0.5).max(160.0))
            .auto_shrink([false, true])
            .show(ui, |ui| {
                egui::Grid::new("outlier-table")
                    .striped(true)
                    .num_columns(6)
                    .spacing([12.0, 6.0])
                    .show(ui, |ui| {
                        for h in ["exclude", "genome", "length", "GC", "nearest genome", "why"] {
                            ui.label(RichText::new(h).small().strong().color(p.text_weak));
                        }
                        ui.end_row();
                        for (i, o) in dialog.list.iter().enumerate() {
                            ui.checkbox(&mut dialog.ticked[i], "");
                            ui.label(RichText::new(&o.name).strong())
                                .on_hover_text(o.reasons.join("\n"));
                            ui.label(fmt::bp(o.length));
                            ui.label(format!("{:.1} %", o.gc * 100.0));
                            if o.nearest_name.is_empty() {
                                ui.label(RichText::new("none").color(p.text_weak));
                            } else {
                                ui.label(format!("{} ({:.3})", o.nearest_name, o.nearest_dist))
                                    .on_hover_text(
                                        "sketch distance to the closest other genome; \
                                     about 1 − identity",
                                    );
                            }
                            ui.label(RichText::new(tags(o, model).join(" · ")).color(p.warn))
                                .on_hover_text(o.reasons.join("\n"));
                            ui.end_row();
                        }
                    });
            });
        ui.horizontal(|ui| {
            if ui.small_button("tick all").clicked() {
                dialog.ticked.iter_mut().for_each(|t| *t = true);
            }
            if ui.small_button("tick none").clicked() {
                dialog.ticked.iter_mut().for_each(|t| *t = false);
            }
            theme::hint(
                ui,
                "Excluded genomes are listed in Settings → Run, where the list can be cleared.",
            );
        });
        ui.add_space(6.0);
        theme::hairline(ui);
        let names = dialog.ticked_names();
        ui.horizontal(|ui| {
            let label = format!(
                "Exclude {} ticked and restart",
                if names.len() == 1 {
                    "the one".to_string()
                } else {
                    names.len().to_string()
                }
            );
            if can_restart {
                let r = ui
                    .add_enabled(!names.is_empty(), theme::primary_button(&label, &p))
                    .on_hover_text(
                        "stop this run, leave the ticked genomes out and start the same run \
                         again at once",
                    );
                if r.clicked() {
                    choice = Some(Choice::Restart(names.clone()));
                }
            } else if ui
                .add_enabled(
                    !names.is_empty(),
                    theme::primary_button("Add ticked to the exclusion list", &p),
                )
                .on_hover_text("the next run leaves them out")
                .clicked()
            {
                choice = Some(Choice::AddOnly(names.clone()));
            }
            if ui
                .button("Continue anyway")
                .on_hover_text(
                    "keep every genome; the Run panel keeps a warning that the results may be \
                     compromised",
                )
                .clicked()
            {
                choice = Some(Choice::Continue);
            }
        });
    });

    match choice {
        Some(Choice::Restart(names)) => {
            app.outlier_dialog = None;
            app.exclude_and_restart(&names, ctx);
        }
        Some(Choice::AddOnly(names)) => {
            app.outlier_dialog = None;
            let n = add_exclusions(&mut app.prefs.params.exclude, &names);
            app.save_prefs();
            app.status = format!("{n} genome(s) added to the exclusion list for the next run");
        }
        Some(Choice::Continue) => {
            app.outlier_dialog = None;
            app.outliers_acknowledged = true;
        }
        None => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn info(name: &str, length: u64, gc: f64, d: f64) -> OutlierInfo {
        OutlierInfo {
            id: 0,
            name: name.into(),
            length,
            gc,
            nearest_dist: d,
            nearest_name: if d < 1.0 { "x".into() } else { String::new() },
            reasons: vec!["engine words: detail".into()],
        }
    }

    #[test]
    fn every_flagged_genome_starts_ticked() {
        let mut d = OutlierDialog::new(vec![info("a", 1, 0.3, 1.0), info("b", 1, 0.3, 1.0)]);
        assert_eq!(d.ticked_names(), vec!["a", "b"]);
        d.ticked[0] = false;
        assert_eq!(d.ticked_names(), vec!["b"]);
    }

    #[test]
    fn exclusions_are_added_once() {
        let mut list = vec!["a".to_string()];
        let n = add_exclusions(&mut list, &["a".into(), "b".into(), " ".into(), "b".into()]);
        assert_eq!(n, 1);
        assert_eq!(list, vec!["a", "b"]);
    }

    #[test]
    fn tags_say_why_in_a_few_words() {
        let mut m = UiModel::default();
        for i in 0..4u32 {
            m.apply(
                0.0,
                crate::state::scan_event(i, &format!("g{i}"), 2_800_000, 10),
            );
        }
        let t = tags(&info("odd", 5_600_000, 0.65, 1.0), &m);
        assert_eq!(t[0], "shares nothing with any other genome");
        assert_eq!(t[1], "length 2.00× median");
        assert_eq!(t[2], "GC +15.0 points");
        let t = tags(&info("odd", 2_800_000, 0.5, 0.14), &m);
        assert_eq!(t, vec!["no close relative (≈86% identity)"]);
        // Nothing the model can explain: the engine's words, shortened.
        let t = tags(&info("odd", 2_800_000, 0.5, 0.01), &m);
        assert_eq!(t, vec!["engine words"]);
    }
}
