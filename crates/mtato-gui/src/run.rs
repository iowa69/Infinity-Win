//! Starting a run and opening a finished one, both entirely off the interface
//! thread (section 8: "All computation happens off the interface thread").
//!
//! The engine talks through a crossbeam channel, and egui only redraws when
//! something asks it to. [`Pump`] is the join between the two: a relay thread
//! blocks on the engine's channel, pushes whatever arrives into a queue and
//! wakes the interface *once per burst* rather than once per event. The update
//! loop then drains at most [`BATCH`] events per frame, so a thousand
//! `GenomeScanned` events cannot starve the drawing and the window cannot stall.

use std::collections::VecDeque;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

use crossbeam_channel::{unbounded, Sender};
use egui::Context;

use mtato_core::engine::{self, RunConfig, RunResults};
use mtato_core::events::{event_bus, Command, Envelope};
use mtato_core::params::Parameters;

/// How many events one frame consumes. At 60 fps this keeps up with about
/// 30 000 events per second while still leaving the frame budget to draw.
pub const BATCH: usize = 512;

/// A thread-safe queue between the relay thread and the update loop.
#[derive(Debug, Default)]
pub struct Pump {
    inner: Mutex<VecDeque<Envelope>>,
}

impl Pump {
    pub fn new() -> Arc<Self> {
        Arc::new(Self::default())
    }

    pub fn push(&self, e: Envelope) {
        if let Ok(mut g) = self.inner.lock() {
            g.push_back(e);
        }
    }

    /// Take at most `max` events, oldest first.
    pub fn drain(&self, max: usize) -> Vec<Envelope> {
        let Ok(mut g) = self.inner.lock() else {
            return Vec::new();
        };
        let n = max.min(g.len());
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            match g.pop_front() {
                Some(e) => out.push(e),
                None => break,
            }
        }
        out
    }

    pub fn len(&self) -> usize {
        self.inner.lock().map(|g| g.len()).unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A run in progress, or one that has just finished and whose results the
/// interface has not collected yet.
pub struct RunHandle {
    pump: Arc<Pump>,
    commands: Sender<Command>,
    results: Arc<Mutex<Option<Result<RunResults, String>>>>,
    running: Arc<AtomicBool>,
    outdir: PathBuf,
    /// A problem worth mentioning that did not stop the run.
    pub startup_note: Option<String>,
}

impl std::fmt::Debug for RunHandle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RunHandle")
            .field("outdir", &self.outdir)
            .field("queued", &self.pump.len())
            .field("running", &self.is_running())
            .finish()
    }
}

impl RunHandle {
    pub fn pump(&self) -> &Arc<Pump> {
        &self.pump
    }

    pub fn outdir(&self) -> &Path {
        &self.outdir
    }

    pub fn is_running(&self) -> bool {
        self.running.load(Ordering::SeqCst)
    }

    fn send(&self, cmd: Command) -> bool {
        self.commands.send(cmd).is_ok()
    }

    /// Move a cluster to the front of the refinement queue (6.11.3). Sent when
    /// the user clicks a cluster in the tree, as section 8 requires.
    pub fn prioritize(&self, cluster: u32) -> bool {
        self.send(Command::Prioritize { cluster })
    }

    /// Re-cluster from the stored sketches. Takes effect immediately, which is
    /// why the cluster distance is a slider and not a text field.
    pub fn set_cluster_dist(&self, value: f64) -> bool {
        self.send(Command::SetClusterDist { value })
    }

    pub fn cancel(&self) -> bool {
        self.send(Command::Cancel)
    }

    /// Collect the results once, after the engine has finished. `None` means
    /// "still running".
    pub fn take_results(&self) -> Option<Result<RunResults, String>> {
        self.results.lock().ok()?.take()
    }
}

/// Validate the parameters, then start the engine and the relay thread.
pub fn start(
    inputs: Vec<PathBuf>,
    params: Parameters,
    outdir: PathBuf,
    ctx: Context,
) -> Result<RunHandle, String> {
    params
        .validate()
        .map_err(|e| format!("invalid parameters: {e}"))?;
    if inputs.is_empty() {
        return Err("nothing to compare: add a genome file or a folder first".to_string());
    }
    let cfg = RunConfig {
        inputs,
        params,
        outdir: Some(outdir.clone()),
    };

    // `spawn` would hide the sink, and the sink is what writes `events.jsonl`;
    // without that file a finished run cannot be opened again later.
    let (sink, rx) = event_bus();
    let startup_note = match sink.attach_log(&outdir) {
        Ok(()) => None,
        Err(e) => Some(format!(
            "cannot write events.jsonl in {}: {e}; the run continues but the folder cannot be \
             reopened later",
            outdir.display()
        )),
    };
    let (commands, cmd_rx) = unbounded::<Command>();
    let pump = Pump::new();
    let results: Arc<Mutex<Option<Result<RunResults, String>>>> = Arc::new(Mutex::new(None));
    let running = Arc::new(AtomicBool::new(true));

    let p2 = Arc::clone(&pump);
    let r2 = Arc::clone(&results);
    let n2 = Arc::clone(&running);
    std::thread::Builder::new()
        .name("mtato-gui-relay".into())
        .spawn(move || {
            let engine = std::thread::Builder::new()
                .name("mtato-gui-engine".into())
                .spawn(move || engine::run(cfg, sink, cmd_rx))
                .expect("failed to spawn the engine thread");
            // The receiver errors once the sink is dropped, which happens when
            // `run` returns: that is the end-of-run signal.
            while let Ok(env) = rx.recv() {
                p2.push(env);
                while let Ok(more) = rx.try_recv() {
                    p2.push(more);
                }
                ctx.request_repaint();
            }
            let res = match engine.join() {
                Ok(r) => r,
                Err(_) => Err("the engine thread panicked".to_string()),
            };
            n2.store(false, Ordering::SeqCst);
            if let Ok(mut g) = r2.lock() {
                *g = Some(res);
            }
            ctx.request_repaint();
        })
        .map_err(|e| format!("cannot start the run: {e}"))?;

    Ok(RunHandle {
        pump,
        commands,
        results,
        running,
        outdir,
        startup_note,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use mtato_core::events::Event;

    fn env(seq: u64) -> Envelope {
        Envelope {
            seq,
            t_ms: seq as f64,
            wall_ms: 0,
            event: Event::Warning {
                text: format!("w{seq}"),
            },
        }
    }

    #[test]
    fn the_pump_keeps_the_order_and_drains_at_most_a_batch() {
        let p = Pump::new();
        for i in 0..1000 {
            p.push(env(i));
        }
        assert_eq!(p.len(), 1000);
        let first = p.drain(BATCH);
        assert_eq!(first.len(), BATCH);
        assert_eq!(first[0].seq, 0);
        assert_eq!(first[BATCH - 1].seq, (BATCH - 1) as u64);
        assert_eq!(p.len(), 1000 - BATCH);
        // The window can never be asked to swallow the whole burst at once.
        assert!(first.len() <= BATCH);
        let rest = p.drain(usize::MAX);
        assert_eq!(rest.len(), 1000 - BATCH);
        assert!(p.is_empty());
        assert!(p.drain(BATCH).is_empty());
    }

    #[test]
    fn draining_more_than_is_queued_is_not_an_error() {
        let p = Pump::new();
        p.push(env(1));
        assert_eq!(p.drain(10).len(), 1);
        assert_eq!(p.drain(0).len(), 0);
    }

    #[test]
    fn a_parameter_set_the_engine_would_reject_is_caught_before_a_thread_starts() {
        let ctx = Context::default();
        let mut params = Parameters::default();
        params.k = 30;
        let err = start(vec![PathBuf::from("/tmp")], params, PathBuf::from("/tmp/out"), ctx)
            .expect_err("k must be odd");
        assert!(err.contains("k must be odd"), "{err}");
    }

    #[test]
    fn running_with_nothing_to_compare_is_refused() {
        let ctx = Context::default();
        let err = start(
            Vec::new(),
            Parameters::default(),
            PathBuf::from("/tmp/out"),
            ctx,
        )
        .expect_err("no inputs");
        assert!(err.contains("nothing to compare"), "{err}");
    }
}
