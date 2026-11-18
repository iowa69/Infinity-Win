//! `mtato-gui`, the desktop program of section 8.
//!
//! One window that shows the analysis sharpening while it runs: a fast blurry
//! tree first, then block-level differences, then exact base-level variants.
//!
//! # Structure
//!
//! Everything that is not drawing lives in a module with no `egui::Ui` in it, so
//! it can be unit-tested without a window:
//!
//! | Module | Role |
//! |---|---|
//! | [`state`] | folds the engine's event stream (section 7) into the UI model |
//! | [`layout`] | turns a tree into node positions, plus the move animation |
//! | [`map`] | the pivot axis: contig offsets, covered spans, row planning |
//! | [`pairwise`] | recomputes 6.9 from two members' variant lists |
//! | [`fmt`] | number formatting and the 0-based to 1-based conversion |
//! | [`colors`] | the visual vocabulary of focus levels, clusters and variants |
//! | [`prefs`] | remembered folders and parameters, as JSON in the config dir |
//! | [`export`] | the files of section 9, written from the UI model |
//! | [`run`] | starting the engine and waking the interface from its threads |
//! | [`open`] | reading a finished output folder back through the same fold |
//!
//! The remaining modules are renderers over that model: [`chrome`] (toolbar and
//! status bar), [`genome_panel`], [`tree_panel`], [`map_panel`] and
//! [`side_panel`]. [`app`] holds them together and owns the repaint policy.
//!
//! # Responsiveness
//!
//! The window must never stall for more than about 50 ms (section 8). Three
//! things guarantee that: the engine runs on its own threads; a relay thread
//! wakes the interface once per *burst* of events rather than once per event;
//! and the update loop drains at most [`run::BATCH`] events per frame, so a
//! thousand `GenomeScanned` events arriving at once cannot starve the drawing.

#![forbid(unsafe_code)]

mod app;
mod chrome;
mod colors;
mod export;
mod fmt;
mod genome_panel;
mod layout;
mod map;
mod map_panel;
mod open;
mod pairwise;
mod prefs;
mod run;
mod side_panel;
mod state;
mod tree_panel;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1480.0, 920.0])
            .with_min_inner_size([940.0, 560.0])
            .with_title(format!("M-Tato {}", mtato_core::VERSION))
            // Genome folders are dropped onto the window: there is no native
            // file dialog, and that is deliberate (see `chrome::toolbar`).
            .with_drag_and_drop(true),
        ..Default::default()
    };
    eframe::run_native(
        "M-Tato",
        options,
        Box::new(|_cc| Ok(Box::new(app::MtatoApp::new()))),
    )
}
