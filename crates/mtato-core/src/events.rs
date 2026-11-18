//! The event stream of section 7.
//!
//! The engine is a library that runs in background threads and talks to its
//! consumers only through these messages (section 0 rule 2). The command-line
//! tool and the GUI are two consumers of exactly the same stream, and every
//! event is also appended to `events.jsonl` with a timestamp so that the
//! benchmark can measure how accuracy grows over time.
//!
//! Coordinates in events are 0-based and half-open, like everywhere inside the
//! program. Files in the output folder use 1-based positions (section 5); that
//! conversion happens in `output.rs` only.

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crossbeam_channel::{unbounded, Receiver, Sender};
use serde::{Deserialize, Serialize};

use crate::counters::Counters;
use crate::summary::RunSummary;
use crate::variant::Variant;

/// What a genome looks like once it has been scanned (6.2).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenomeScannedInfo {
    pub id: u32,
    pub name: String,
    pub length: u64,
    pub contigs: u32,
    pub n50: u64,
    pub gc: f64,
    pub sketch_size: usize,
}

/// One cluster as the consumer sees it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClusterInfo {
    pub id: u32,
    pub pivot: u32,
    pub pivot_name: String,
    pub members: Vec<u32>,
    pub member_names: Vec<String>,
    /// Best focus this cluster has reached so far.
    pub focus: u8,
    /// Sketch distance threshold in force when the cluster was formed.
    pub cluster_dist: f64,
}

/// The per-genome numbers shown next to a leaf in the GUI.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MemberNumbers {
    pub id: u32,
    pub name: String,
    /// Difference pieces at Focus 1, or resolved variants at Focus 2.
    pub diff_pieces: u64,
    pub snps: u64,
    pub indels: u64,
    pub uncovered_bp: u64,
    pub novel_bp: u64,
}

/// A cluster tree at some focus level (6.9).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClusterFocusInfo {
    pub cluster: u32,
    /// 0 = sketch only, 1 = blocks compared, 2 = exact base level (section 5).
    pub focus: u8,
    /// True when the distances behind this tree are estimates and must be
    /// shown with the "≈" marker rather than as exact numbers (section 0
    /// rule 8).
    pub approximate: bool,
    pub newick: String,
    pub members: Vec<MemberNumbers>,
}

/// One difference piece in pivot coordinates, before it is resolved.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PieceView {
    pub pivot_contig: u32,
    pub pivot_contig_name: String,
    pub start: u64,
    pub end: u64,
    /// True once the piece has been resolved into variants (Focus 2).
    pub resolved: bool,
}

/// Everything the map panel needs for one member of one cluster.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MapUpdateInfo {
    pub cluster: u32,
    pub genome: u32,
    pub genome_name: String,
    pub pivot_name: String,
    pub pieces: Vec<PieceView>,
    pub variants: Vec<Variant>,
    /// Pivot intervals this member does not cover (6.8).
    pub uncovered: Vec<(u64, u64)>,
    /// Pivot contig names, in order, so the map can lay out its axis.
    pub pivot_contigs: Vec<String>,
    pub pivot_contig_lens: Vec<u64>,
}

/// The stages a run passes through, named for the status bar.
pub mod stage {
    pub const DISCOVER: &str = "discover";
    pub const SCAN: &str = "scan";
    pub const CLUSTER: &str = "cluster";
    pub const BACKBONE: &str = "backbone";
    pub const REFINE: &str = "refine";
    pub const OUTPUT: &str = "output";
}

/// Engine to consumer (section 7).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "event", content = "data", rename_all = "snake_case")]
pub enum Event {
    GenomeScanned(GenomeScannedInfo),
    Clustered {
        threshold: f64,
        clusters: Vec<ClusterInfo>,
    },
    BackboneTree {
        newick: String,
        /// Names of the pivot genomes, one per leaf, in Newick order.
        pivots: Vec<String>,
    },
    ClusterFocus(ClusterFocusInfo),
    MapUpdate(Box<MapUpdateInfo>),
    Counters(Counters),
    Progress {
        stage: String,
        done: u64,
        total: u64,
        elapsed_ms: u64,
    },
    Warning {
        text: String,
    },
    Finished(RunSummary),
}

/// Consumer to engine (section 7).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "command", rename_all = "snake_case")]
pub enum Command {
    /// Move this cluster to the front of the refinement queue.
    Prioritize { cluster: u32 },
    /// Re-cluster from the stored sketches with a new threshold.
    SetClusterDist { value: f64 },
    /// Stop as soon as the current unit of work finishes.
    Cancel,
}

/// An event with the timing the log file needs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Envelope {
    /// Monotonic sequence number, so a consumer can tell a dropped message
    /// from a late one.
    pub seq: u64,
    /// Milliseconds since the run started. This is the x axis of the
    /// accuracy-over-time figure of section 11.
    pub t_ms: f64,
    /// Unix milliseconds, for correlating with other logs.
    pub wall_ms: u64,
    pub event: Event,
}

impl Envelope {
    pub fn event_name(&self) -> &'static str {
        match &self.event {
            Event::GenomeScanned(_) => "genome_scanned",
            Event::Clustered { .. } => "clustered",
            Event::BackboneTree { .. } => "backbone_tree",
            Event::ClusterFocus(_) => "cluster_focus",
            Event::MapUpdate(_) => "map_update",
            Event::Counters(_) => "counters",
            Event::Progress { .. } => "progress",
            Event::Warning { .. } => "warning",
            Event::Finished(_) => "finished",
        }
    }
}

struct SinkInner {
    tx: Sender<Envelope>,
    start: Instant,
    wall_start: SystemTime,
    seq: AtomicU64,
    /// Optional `events.jsonl` writer. A mutex, because events are emitted
    /// from rayon workers.
    log: Mutex<Option<BufWriter<File>>>,
}

/// Cloneable, thread-safe handle the engine uses to talk to its consumers.
#[derive(Clone)]
pub struct EventSink {
    inner: Arc<SinkInner>,
}

impl std::fmt::Debug for EventSink {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EventSink")
            .field("emitted", &self.emitted())
            .finish()
    }
}

impl EventSink {
    /// A sink with no log file. Pair it with [`event_bus`] to receive events.
    pub fn new(tx: Sender<Envelope>) -> Self {
        Self {
            inner: Arc::new(SinkInner {
                tx,
                start: Instant::now(),
                wall_start: SystemTime::now(),
                seq: AtomicU64::new(0),
                log: Mutex::new(None),
            }),
        }
    }

    /// A sink that drops everything. Useful in tests and for the simulator.
    pub fn null() -> Self {
        let (tx, rx) = unbounded();
        drop(rx);
        Self::new(tx)
    }

    /// Start appending every event to `events.jsonl` inside `dir`.
    ///
    /// Failing to open the log is not fatal: the run continues and the caller
    /// gets the error back to report as a warning, because a missing log file
    /// should not throw away an analysis.
    pub fn attach_log(&self, dir: &Path) -> std::io::Result<()> {
        std::fs::create_dir_all(dir)?;
        let path = dir.join("events.jsonl");
        let f = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        let mut guard = self.inner.log.lock().expect("event log mutex poisoned");
        *guard = Some(BufWriter::new(f));
        Ok(())
    }

    /// Flush the log file. Called at the end of a run and before any point
    /// where a consumer might read `events.jsonl`.
    pub fn flush_log(&self) {
        if let Ok(mut guard) = self.inner.log.lock() {
            if let Some(w) = guard.as_mut() {
                let _ = w.flush();
            }
        }
    }

    /// Emit one event: timestamp it, number it, log it and send it.
    pub fn emit(&self, event: Event) {
        let seq = self.inner.seq.fetch_add(1, Ordering::Relaxed);
        let t_ms = self.inner.start.elapsed().as_secs_f64() * 1000.0;
        let wall_ms = self
            .inner
            .wall_start
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_millis() as u64)
            .unwrap_or(0);
        let env = Envelope {
            seq,
            t_ms,
            wall_ms,
            event,
        };
        if let Ok(mut guard) = self.inner.log.lock() {
            if let Some(w) = guard.as_mut() {
                if let Ok(line) = serde_json::to_string(&env) {
                    let _ = writeln!(w, "{line}");
                }
            }
        }
        // A disconnected receiver is normal: the CLI may have stopped
        // listening after `Finished`. Emitting must never panic.
        let _ = self.inner.tx.send(env);
    }

    /// Milliseconds since the run started, for `Progress` events.
    pub fn elapsed_ms(&self) -> u64 {
        self.inner.start.elapsed().as_millis() as u64
    }

    /// Number of events emitted so far.
    pub fn emitted(&self) -> u64 {
        self.inner.seq.load(Ordering::Relaxed)
    }

    // Convenience emitters, so call sites read like the specification.

    pub fn warn(&self, text: impl Into<String>) {
        self.emit(Event::Warning { text: text.into() });
    }

    pub fn progress(&self, stage: &str, done: u64, total: u64) {
        self.emit(Event::Progress {
            stage: stage.to_string(),
            done,
            total,
            elapsed_ms: self.elapsed_ms(),
        });
    }
}

/// Build a connected `(sink, receiver)` pair.
pub fn event_bus() -> (EventSink, Receiver<Envelope>) {
    let (tx, rx) = unbounded();
    (EventSink::new(tx), rx)
}

/// Collect every event of a finished run, for tests and the CLI report.
pub fn drain(rx: &Receiver<Envelope>) -> Vec<Envelope> {
    let mut out = Vec::new();
    while let Ok(e) = rx.try_recv() {
        out.push(e);
    }
    out
}

/// Read an `events.jsonl` file back. This is what "open results" uses to show
/// a finished run without recomputing it (section 8).
pub fn read_event_log(path: &Path) -> std::io::Result<Vec<Envelope>> {
    use std::io::BufRead;
    let f = File::open(path)?;
    let mut out = Vec::new();
    for (i, line) in std::io::BufReader::new(f).lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        match serde_json::from_str::<Envelope>(&line) {
            Ok(e) => out.push(e),
            Err(e) => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("{} line {}: {e}", path.display(), i + 1),
                ))
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn events_are_numbered_and_timed_in_order() {
        let (sink, rx) = event_bus();
        for i in 0..5 {
            sink.warn(format!("w{i}"));
        }
        let got = drain(&rx);
        assert_eq!(got.len(), 5);
        for (i, e) in got.iter().enumerate() {
            assert_eq!(e.seq, i as u64);
        }
        assert!(got.windows(2).all(|w| w[1].t_ms >= w[0].t_ms));
    }

    #[test]
    fn a_disconnected_receiver_does_not_panic() {
        let (sink, rx) = event_bus();
        drop(rx);
        sink.warn("nobody listening");
        assert_eq!(sink.emitted(), 1);
    }

    #[test]
    fn every_event_serialises_with_its_name() {
        let events = vec![
            Event::GenomeScanned(GenomeScannedInfo {
                id: 0,
                name: "g".into(),
                length: 10,
                contigs: 1,
                n50: 10,
                gc: 0.5,
                sketch_size: 3,
            }),
            Event::Clustered {
                threshold: 0.005,
                clusters: vec![],
            },
            Event::BackboneTree {
                newick: "(a,b);".into(),
                pivots: vec!["a".into(), "b".into()],
            },
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(a,b);".into(),
                members: vec![],
            }),
            Event::Counters(Counters::default()),
            Event::Progress {
                stage: "scan".into(),
                done: 1,
                total: 2,
                elapsed_ms: 3,
            },
            Event::Warning { text: "w".into() },
            Event::Finished(RunSummary::new("t", crate::params::Parameters::default())),
        ];
        for e in &events {
            let env = Envelope {
                seq: 0,
                t_ms: 0.0,
                wall_ms: 0,
                event: e.clone(),
            };
            let j = serde_json::to_string(&env).expect("serialises");
            let back: Envelope = serde_json::from_str(&j).expect("deserialises");
            assert_eq!(back, env, "round trip must be exact for {j}");
            assert!(!env.event_name().is_empty());
        }
    }

    #[test]
    fn event_log_round_trips_through_disk() {
        let dir = std::env::temp_dir().join(format!("mtato-events-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (sink, _rx) = event_bus();
        sink.attach_log(&dir).unwrap();
        sink.warn("hello");
        sink.progress("scan", 1, 1);
        sink.flush_log();
        let back = read_event_log(&dir.join("events.jsonl")).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[0].event_name(), "warning");
        assert_eq!(back[1].event_name(), "progress");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn map_update_is_boxed_but_behaves_like_any_other_event() {
        let e = Event::MapUpdate(Box::new(MapUpdateInfo {
            cluster: 0,
            genome: 1,
            genome_name: "m".into(),
            pivot_name: "p".into(),
            pieces: vec![PieceView {
                pivot_contig: 0,
                pivot_contig_name: "c0".into(),
                start: 10,
                end: 20,
                resolved: false,
            }],
            variants: vec![],
            uncovered: vec![(0, 5)],
            pivot_contigs: vec!["c0".into()],
            pivot_contig_lens: vec![100],
        }));
        let j = serde_json::to_string(&e).unwrap();
        assert_eq!(serde_json::from_str::<Event>(&j).unwrap(), e);
    }
}
