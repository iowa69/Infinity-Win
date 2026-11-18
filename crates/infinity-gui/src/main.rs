//! `Infinity-Win`, the desktop program of section 8.
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
//! | [`settings`] | the Resolution and SNP-rule settings, mirrored for the engine |
//! | [`theme`] | palette, spacing, typography and the shared building blocks |
//!
//! The remaining modules are renderers over that model: [`chrome`] (toolbar and
//! status bar), [`genome_panel`], [`tree_panel`], [`map_panel`],
//! [`side_panel`] (the Run panel) and [`settings_panel`] (the Settings window).
//! [`app`] holds them together and owns the repaint policy.
//!
//! # Responsiveness
//!
//! The window must never stall for more than about 50 ms (section 8). Three
//! things guarantee that: the engine runs on its own threads; a relay thread
//! wakes the interface once per *burst* of events rather than once per event;
//! and the update loop drains at most [`run::BATCH`] events per frame, so a
//! thousand `GenomeScanned` events arriving at once cannot starve the drawing.
//!
//! # Headless check
//!
//! `Infinity-Win --check <output-folder>` folds a finished run and prints what the
//! window would show — genomes, clusters, focus counts, map rows, tree rows
//! collapsed and opened — then verifies the layout is finite and complete. It
//! needs no display, so the reading half of the program can be exercised in a
//! terminal or in CI.

// On Windows, a release build is a GUI program: no console window opens behind it
// when Infinity-Win.exe is double-clicked. Debug builds keep the console for logs.
#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]
#![forbid(unsafe_code)]

use eframe::egui;

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
mod settings;
mod settings_panel;
mod side_panel;
mod state;
mod theme;
mod tree_panel;

fn main() -> eframe::Result {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.as_slice() {
        [] => {}
        // A headless smoke test of the whole reading path: fold a finished
        // output folder and report what the window would have shown. No display
        // is needed, so it also runs in CI and over ssh.
        [flag, path] if flag == "--check" => {
            std::process::exit(check_folder(std::path::PathBuf::from(path)));
        }
        _ => {
            eprintln!("usage: Infinity-Win [--check <output-folder>]");
            std::process::exit(2);
        }
    }

    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1480.0, 920.0])
            .with_min_inner_size([940.0, 560.0])
            .with_title(format!("Infinity-Win {}", infinity_core::VERSION))
            // Genome folders are dropped onto the window: there is no native
            // file dialog, and that is deliberate (see `chrome::toolbar`).
            .with_drag_and_drop(true),
        ..Default::default()
    };
    eframe::run_native(
        "Infinity-Win",
        options,
        Box::new(|_cc| Ok(Box::new(app::InfinityApp::new()))),
    )
}

/// Read `dir` as the window would and print what it found. Returns a process
/// exit code, so a regression shows up as a non-zero status.
fn check_folder(dir: std::path::PathBuf) -> i32 {
    use std::collections::HashSet;
    use std::sync::atomic::AtomicU64;

    let read_bytes = AtomicU64::new(0);
    let model = match open::read_folder(&dir, &read_bytes, &egui::Context::default()) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Infinity-Win --check: {e}");
            return 1;
        }
    };
    println!("folder            {}", dir.display());
    println!("events            {}", model.events_seen);
    println!("genomes           {}", model.genome_count());
    println!("clusters          {}", model.cluster_count());
    println!(
        "focus 0/1/2       {}/{}/{}",
        model
            .genomes
            .iter()
            .filter(|g| g.focus == fmt::FOCUS_SKETCH)
            .count(),
        model
            .genomes
            .iter()
            .filter(|g| g.focus == fmt::FOCUS_BLOCKS)
            .count(),
        model
            .genomes
            .iter()
            .filter(|g| g.focus == fmt::FOCUS_EXACT)
            .count()
    );
    println!("stages            {}", model.stages.len());
    println!("maps              {}", count_maps(&model));
    println!("warnings          {}", model.notes.len());

    // Everything collapsed, then everything opened: both must lay out, and the
    // opened tree must have one row per genome.
    let collapsed = match model.display_tree(&HashSet::new()) {
        Some(d) => layout::layout(&d),
        None => {
            eprintln!("Infinity-Win --check: no backbone tree in the recording");
            return 1;
        }
    };
    let all: HashSet<u32> = model.clusters.iter().map(|c| c.id).collect();
    let opened = match model.display_tree(&all) {
        Some(d) => layout::layout(&d),
        None => {
            eprintln!("Infinity-Win --check: no backbone tree in the recording");
            return 1;
        }
    };
    println!(
        "tree rows         {} collapsed / {} opened",
        collapsed.rows as usize, opened.rows as usize
    );
    println!(
        "tree nodes        {} / {}   regions {}",
        collapsed.nodes.len(),
        opened.nodes.len(),
        opened.regions.len()
    );
    let genomes_with_maps = model
        .clusters
        .iter()
        .filter(|c| !model.mapped_members(c.id).is_empty())
        .count();
    println!("clusters with map {}", genomes_with_maps);

    let mut bad = 0;
    if opened.rows as usize != model.genome_count() {
        eprintln!(
            "FAIL: the fully opened tree has {} leaf rows but there are {} genomes",
            opened.rows as usize,
            model.genome_count()
        );
        bad += 1;
    }
    if !opened
        .nodes
        .iter()
        .all(|n| n.x.is_finite() && n.y.is_finite())
    {
        eprintln!("FAIL: a node was laid out at a non-finite position");
        bad += 1;
    }
    if collapsed.rows as usize > opened.rows as usize {
        eprintln!("FAIL: collapsing clusters increased the row count");
        bad += 1;
    }
    if bad == 0 {
        println!("OK");
        0
    } else {
        1
    }
}

fn count_maps(model: &state::UiModel) -> usize {
    model.maps.values().map(|m| m.len()).sum()
}
