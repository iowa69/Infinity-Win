//! Opening a finished output folder and viewing it without recomputing
//! (section 8).
//!
//! A run appends every event, with its timestamp, to `events.jsonl` (section 7),
//! so the folder *is* a recording of the run. Reading it back through exactly the
//! same fold the live window uses is what guarantees that an opened result and a
//! watched one look identical.
//!
//! The lines are streamed rather than collected with `read_event_log`, because a
//! thousand-genome run logs one `MapUpdate` per member with every variant in it:
//! holding all of those envelopes in memory at once, on top of the model, is
//! avoidable. Streaming also gives the progress bar something to show.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use egui::Context;

use mtato_core::events::Envelope;
use mtato_core::summary::RunSummary;

use crate::state::{Source, UiModel};

/// Name of the recording inside an output folder.
pub const LOG: &str = "events.jsonl";
/// Name of the summary inside an output folder.
pub const RUN_JSON: &str = "run.json";
/// How often the reader wakes the interface to report progress.
const WAKE_EVERY: u64 = 256;

/// A folder being read on a background thread.
pub struct OpenHandle {
    model: Arc<Mutex<Option<UiModel>>>,
    error: Arc<Mutex<Option<String>>>,
    read_bytes: Arc<AtomicU64>,
    total_bytes: u64,
    dir: PathBuf,
}

impl std::fmt::Debug for OpenHandle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("OpenHandle")
            .field("dir", &self.dir)
            .field("total_bytes", &self.total_bytes)
            .field("read_bytes", &self.read_bytes.load(Ordering::Relaxed))
            .finish()
    }
}

impl OpenHandle {
    /// Fraction read, in `[0, 1]`, or `None` when the file's size is unknown.
    pub fn progress(&self) -> Option<f32> {
        if self.total_bytes == 0 {
            return None;
        }
        let r = self.read_bytes.load(Ordering::Relaxed);
        Some((r as f32 / self.total_bytes as f32).clamp(0.0, 1.0))
    }

    pub fn dir(&self) -> &Path {
        &self.dir
    }

    pub fn is_done(&self) -> bool {
        self.done().is_some()
    }

    /// `Some(Ok(model))` or `Some(Err(message))` once the reader has finished,
    /// `None` while it is still going. Collected once.
    pub fn done(&self) -> Option<Result<UiModel, String>> {
        if let Ok(mut g) = self.model.lock() {
            if g.is_some() {
                return Some(Ok(g.take().expect("checked above")));
            }
        }
        if let Ok(g) = self.error.lock() {
            if let Some(e) = g.as_ref() {
                return Some(Err(e.clone()));
            }
        }
        None
    }
}

/// Start reading `dir` on a background thread.
pub fn start(dir: PathBuf, ctx: Context) -> Result<OpenHandle, String> {
    let log = dir.join(LOG);
    if !log.is_file() {
        return Err(format!(
            "{} has no {}: not an M-Tato output folder",
            dir.display(),
            LOG
        ));
    }
    let total_bytes = std::fs::metadata(&log).map(|m| m.len()).unwrap_or(0);
    let handle = OpenHandle {
        model: Arc::new(Mutex::new(None)),
        error: Arc::new(Mutex::new(None)),
        read_bytes: Arc::new(AtomicU64::new(0)),
        total_bytes,
        dir: dir.clone(),
    };
    let m2 = Arc::clone(&handle.model);
    let e2 = Arc::clone(&handle.error);
    let b2 = Arc::clone(&handle.read_bytes);
    std::thread::Builder::new()
        .name("mtato-gui-open".into())
        .spawn(move || {
            match read_folder(&dir, &b2, &ctx) {
                Ok(model) => {
                    if let Ok(mut g) = m2.lock() {
                        *g = Some(model);
                    }
                }
                Err(err) => {
                    if let Ok(mut g) = e2.lock() {
                        *g = Some(err);
                    }
                }
            }
            ctx.request_repaint();
        })
        .map_err(|e| format!("cannot start reading {}: {e}", dir.display()))?;
    Ok(handle)
}

/// Stream `events.jsonl` through the fold. Exposed for tests.
pub fn read_folder(dir: &Path, read_bytes: &AtomicU64, ctx: &Context) -> Result<UiModel, String> {
    let path = dir.join(LOG);
    let file = File::open(&path).map_err(|e| format!("opening {}: {e}", path.display()))?;
    let mut model = UiModel::new(Source::Folder(dir.to_path_buf()));
    let mut read = 0u64;
    let mut lines = 0u64;
    let mut bad = 0u64;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| format!("reading {}: {e}", path.display()))?;
        read += line.len() as u64 + 1;
        lines += 1;
        if line.trim().is_empty() {
            continue;
        }
        match serde_json::from_str::<Envelope>(&line) {
            Ok(env) => model.apply_envelope(env),
            Err(e) => {
                bad += 1;
                if bad <= 3 {
                    model.push_local_note(format!("skipping an unreadable event: {e}"));
                }
            }
        }
        if lines % WAKE_EVERY == 0 {
            read_bytes.store(read, Ordering::Relaxed);
            ctx.request_repaint();
        }
    }
    read_bytes.store(read, Ordering::Relaxed);
    if model.events_seen == 0 {
        return Err(format!("{} is empty", path.display()));
    }
    if bad > 3 {
        model.push_local_note(format!("skipped {bad} unreadable events in total"));
    }
    // `run.json` carries the parameters when the recording was truncated by a
    // crash before the `Finished` event.
    if model.summary.is_none() {
        if let Ok(text) = std::fs::read_to_string(dir.join(RUN_JSON)) {
            if let Ok(s) = serde_json::from_str::<RunSummary>(&text) {
                model.apply(0.0, mtato_core::events::Event::Finished(Box::new(s)));
            }
        }
    }
    Ok(model)
}

#[cfg(test)]
mod tests {
    use super::*;
    use mtato_core::events::{event_bus, Event, GenomeScannedInfo};

    fn tmp_dir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("mtato-gui-open-{}-{tag}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// Write a small but real recording using the engine's own sink, so the
    /// reader is tested against the format the engine actually produces.
    fn write_recording(dir: &Path, genomes: u32) {
        let (sink, _rx) = event_bus();
        sink.attach_log(dir).unwrap();
        for i in 0..genomes {
            sink.emit(Event::GenomeScanned(GenomeScannedInfo {
                id: i,
                name: format!("g{i}"),
                length: 1000 + u64::from(i),
                contigs: 1,
                n50: 1000,
                gc: 0.5,
                sketch_size: 10,
            }));
        }
        sink.warn("a recorded warning");
        sink.progress(mtato_core::events::stage::SCAN, genomes as u64, genomes as u64);
        sink.flush_log();
    }

    #[test]
    fn an_opened_folder_folds_back_into_the_same_model() {
        let dir = tmp_dir("basic");
        write_recording(&dir, 3);
        let counter = AtomicU64::new(0);
        let ctx = Context::default();
        let model = read_folder(&dir, &counter, &ctx).expect("reads");
        assert_eq!(model.genome_count(), 3);
        assert_eq!(model.genome(2).unwrap().length, 1002);
        assert!(model.notes.iter().any(|n| n.text == "a recorded warning"));
        assert_eq!(model.stages.len(), 1);
        assert!(model.stages[0].finished);
        assert!(counter.load(Ordering::Relaxed) > 0);
        assert_eq!(model.source, Source::Folder(dir.clone()));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_folder_without_a_recording_is_reported_clearly() {
        let dir = tmp_dir("empty-dir");
        let err = start(dir.clone(), Context::default()).expect_err("no events.jsonl");
        assert!(err.contains("not an M-Tato output folder"), "{err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_empty_recording_is_an_error_and_a_broken_line_is_skipped() {
        let dir = tmp_dir("empty-log");
        std::fs::write(dir.join(LOG), "").unwrap();
        let counter = AtomicU64::new(0);
        let ctx = Context::default();
        assert!(read_folder(&dir, &counter, &ctx).is_err());

        let dir2 = tmp_dir("broken-log");
        write_recording(&dir2, 2);
        // Append garbage, as a run killed mid-write would leave behind.
        let mut f = std::fs::OpenOptions::new()
            .append(true)
            .open(dir2.join(LOG))
            .unwrap();
        std::io::Write::write_all(&mut f, b"{\"event\":").unwrap();
        drop(f);
        let model = read_folder(&dir2, &counter, &ctx).expect("the good lines still load");
        assert_eq!(model.genome_count(), 2);
        assert!(model.notes.iter().any(|n| n.text.contains("unreadable")));

        for d in [dir, dir2] {
            std::fs::remove_dir_all(d).ok();
        }
    }

    #[test]
    fn run_json_supplies_a_summary_the_recording_lost() {
        let dir = tmp_dir("run-json");
        write_recording(&dir, 1);
        let s = RunSummary::new("test", mtato_core::params::Parameters::default());
        std::fs::write(dir.join(RUN_JSON), serde_json::to_string_pretty(&s).unwrap()).unwrap();
        let counter = AtomicU64::new(0);
        let model = read_folder(&dir, &counter, &Context::default()).unwrap();
        assert!(model.summary.is_some());
        assert!(model.finished);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn the_reader_finishes_on_its_own_and_hands_the_model_over_once() {
        let dir = tmp_dir("progress");
        write_recording(&dir, 1);
        let h = start(dir.clone(), Context::default()).expect("starts");
        assert!(h.total_bytes > 0);
        if let Some(p) = h.progress() {
            assert!((0.0..=1.0).contains(&p));
        }
        // The folder is tiny, so the reader is done almost at once.
        for _ in 0..200 {
            if h.is_done() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        let done = h.done().expect("finished reading");
        assert_eq!(done.unwrap().genome_count(), 1);
        assert!(h.done().is_none(), "the model is collected only once");
        std::fs::remove_dir_all(&dir).ok();
    }
}
