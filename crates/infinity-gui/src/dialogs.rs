//! The native "Add files…" and "Add folder…" dialogs, and what the genome
//! drawer says about each input (how many genomes a folder holds).
//!
//! The dialogs come from `rfd`: the Win32 common dialog on Windows, the
//! desktop's own file chooser through xdg-desktop-portal on Linux (no GTK, see
//! `DECISIONS.md` D3). A file dialog is modal and can stay open for minutes, so
//! it never runs on the UI thread: [`start`] creates `rfd`'s asynchronous
//! dialog (which needs the window handle, so it is built on the UI thread) and
//! hands the future to a worker thread, which drives it to completion with a
//! ten-line executor and wakes the window. The update loop polls
//! [`PendingPick::take`] once per frame.

use std::future::Future;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::task::{Context, Poll, Wake};

/// What the dialog picks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PickKind {
    /// One or more genome files.
    Files,
    /// One or more folders of genome files.
    Folders,
}

impl PickKind {
    pub fn title(self) -> &'static str {
        match self {
            PickKind::Files => "Add genome files",
            PickKind::Folders => "Add folders of genomes",
        }
    }
}

/// The extensions the file dialog offers by default. An explicitly chosen file
/// is always read, whatever its name; a folder contributes only the files the
/// engine accepts (see [`count_genomes`]).
pub const GENOME_EXTENSIONS: &[&str] = &[
    "fasta", "fa", "fna", "fas", "fasta.gz", "fa.gz", "fna.gz", "fas.gz",
];

/// A dialog that is open, or has just closed.
#[derive(Debug)]
pub struct PendingPick {
    pub kind: PickKind,
    slot: Arc<Mutex<Option<Vec<PathBuf>>>>,
}

impl PendingPick {
    /// `Some(paths)` once the dialog has closed (empty when it was cancelled),
    /// `None` while it is still open. Collected once.
    pub fn take(&self) -> Option<Vec<PathBuf>> {
        self.slot.lock().ok()?.take()
    }
}

/// Open a dialog without blocking the caller. `parent` makes it modal to the
/// main window where the platform supports it; `dir` is where it starts.
pub fn start(
    kind: PickKind,
    parent: Option<&eframe::Frame>,
    dir: Option<&Path>,
    ctx: egui::Context,
) -> Result<PendingPick, String> {
    let mut d = rfd::AsyncFileDialog::new().set_title(kind.title());
    if let Some(p) = parent {
        d = d.set_parent(p);
    }
    if let Some(dir) = dir.filter(|d| d.is_dir()) {
        d = d.set_directory(dir);
    }
    let slot: Arc<Mutex<Option<Vec<PathBuf>>>> = Arc::new(Mutex::new(None));
    let out = Arc::clone(&slot);
    let finish = move |paths: Vec<PathBuf>| {
        if let Ok(mut g) = out.lock() {
            *g = Some(paths);
        }
        ctx.request_repaint();
    };
    let spawned = match kind {
        PickKind::Files => {
            let fut = d
                .add_filter("Genome assemblies (FASTA)", GENOME_EXTENSIONS)
                .add_filter("All files", &["*"])
                .pick_files();
            std::thread::Builder::new()
                .name("infinity-gui-dialog".into())
                .spawn(move || finish(handles_to_paths(block_on(fut))))
        }
        PickKind::Folders => {
            let fut = d.pick_folders();
            std::thread::Builder::new()
                .name("infinity-gui-dialog".into())
                .spawn(move || finish(handles_to_paths(block_on(fut))))
        }
    };
    spawned.map_err(|e| format!("cannot open the file dialog: {e}"))?;
    Ok(PendingPick { kind, slot })
}

fn handles_to_paths(picked: Option<Vec<rfd::FileHandle>>) -> Vec<PathBuf> {
    picked
        .unwrap_or_default()
        .iter()
        .map(|h| h.path().to_path_buf())
        .collect()
}

/// Wakes the worker thread that is waiting on the dialog.
struct ThreadWaker(std::thread::Thread);

impl Wake for ThreadWaker {
    fn wake(self: Arc<Self>) {
        self.0.unpark();
    }
}

/// Drive a future to completion on the current thread, parking between polls.
/// `rfd` resolves its dialog future from its own thread, so this is all the
/// executor the dialog needs.
fn block_on<F: Future>(fut: F) -> F::Output {
    let mut fut = std::pin::pin!(fut);
    let waker = Arc::new(ThreadWaker(std::thread::current())).into();
    let mut cx = Context::from_waker(&waker);
    loop {
        match fut.as_mut().poll(&mut cx) {
            Poll::Ready(v) => return v,
            Poll::Pending => std::thread::park(),
        }
    }
}

/// How many genomes an input contributes, by the engine's own rule: a folder
/// is every accepted file directly inside it, a file is one genome. `None`
/// when the path cannot be read (it was moved, or a drive was unplugged).
pub fn count_genomes(path: &Path) -> Option<usize> {
    if path.is_file() {
        return Some(1);
    }
    infinity_core::fasta::collect_inputs(&[path.to_path_buf()])
        .ok()
        .map(|v| v.len())
}

/// "182 genomes", "1 genome", "no genomes", or "cannot read".
pub fn count_text(n: Option<usize>) -> String {
    match n {
        None => "cannot read".to_string(),
        Some(0) => "no genomes".to_string(),
        Some(1) => "1 genome".to_string(),
        Some(n) => format!("{} genomes", crate::fmt::count(n as u64)),
    }
}

/// Where the next dialog should start: the folder of the most recent input.
pub fn start_dir(recent: &[String]) -> Option<PathBuf> {
    let p = PathBuf::from(recent.first()?);
    if p.is_dir() {
        Some(p)
    } else {
        p.parent().map(Path::to_path_buf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_executor_drives_a_ready_and_a_late_future() {
        assert_eq!(block_on(async { 7 }), 7);
        // A future completed from another thread, as rfd's is.
        let slot = Arc::new(Mutex::new(None::<u32>));
        let s2 = Arc::clone(&slot);
        let waker_slot: Arc<Mutex<Option<std::task::Waker>>> = Arc::new(Mutex::new(None));
        let w2 = Arc::clone(&waker_slot);
        std::thread::spawn(move || {
            std::thread::sleep(std::time::Duration::from_millis(20));
            *s2.lock().unwrap() = Some(42);
            if let Some(w) = w2.lock().unwrap().take() {
                w.wake();
            }
        });
        let got = block_on(std::future::poll_fn(|cx| {
            if let Some(v) = *slot.lock().unwrap() {
                return Poll::Ready(v);
            }
            *waker_slot.lock().unwrap() = Some(cx.waker().clone());
            // The value may have landed between the check and the waker.
            match *slot.lock().unwrap() {
                Some(v) => Poll::Ready(v),
                None => Poll::Pending,
            }
        }));
        assert_eq!(got, 42);
    }

    #[test]
    fn a_folder_counts_the_genome_files_the_engine_would_read() {
        let dir = std::env::temp_dir().join(format!("infinity-dlg-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for f in ["a.fasta", "b.fna", "c.fa.gz", "notes.txt", "d.fa"] {
            std::fs::write(dir.join(f), ">x\nACGT\n").unwrap();
        }
        assert_eq!(count_genomes(&dir), Some(4));
        assert_eq!(count_genomes(&dir.join("a.fasta")), Some(1));
        assert_eq!(count_genomes(&dir.join("missing")), None);
        assert_eq!(count_text(Some(4)), "4 genomes");
        assert_eq!(count_text(Some(1)), "1 genome");
        assert_eq!(count_text(Some(0)), "no genomes");
        assert_eq!(count_text(None), "cannot read");
        assert_eq!(
            start_dir(&[dir.join("a.fasta").display().to_string()]),
            Some(dir.clone())
        );
        assert_eq!(start_dir(&[dir.display().to_string()]), Some(dir.clone()));
        assert_eq!(start_dir(&[]), None);
        std::fs::remove_dir_all(&dir).ok();
    }
}
