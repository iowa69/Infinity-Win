//! Shared helpers for the integration tests of sections 10.2 and 10.3.
//!
//! These are integration tests rather than unit tests on purpose: they drive the
//! engine through its public entry point, [`mtato_core::engine::run`], reading
//! and writing real files, so they exercise exactly the path the CLI and the GUI
//! use.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use mtato_core::engine::{run, RunConfig, RunResults};
use mtato_core::events::{event_bus, Command, Envelope};
use mtato_core::params::Parameters;

/// A scratch directory that removes itself when dropped.
///
/// Tests must not leave gigabytes of simulated genomes behind, and they must
/// not collide with each other when run in parallel.
pub struct TempDir {
    pub path: PathBuf,
}

impl TempDir {
    pub fn new(label: &str) -> Self {
        let mut path = std::env::temp_dir();
        let unique = format!(
            "mtato-it-{}-{}-{:?}",
            label,
            std::process::id(),
            std::thread::current().name().unwrap_or("main")
        );
        // Thread names can contain characters that are awkward in a path.
        let unique: String = unique
            .chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '-' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        path.push(unique);
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).expect("create temp dir");
        Self { path }
    }

    pub fn join(&self, name: &str) -> PathBuf {
        self.path.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// Write one FASTA file per genome and return their paths, sorted by name so
/// that the run's input order is deterministic.
pub fn write_genomes(dir: &Path, genomes: &[(&str, Vec<String>)]) -> Vec<PathBuf> {
    fs::create_dir_all(dir).expect("create genome dir");
    let mut paths = Vec::new();
    for (name, contigs) in genomes {
        let p = dir.join(format!("{name}.fasta"));
        let mut out = String::new();
        for (i, seq) in contigs.iter().enumerate() {
            out.push_str(&format!(">{name}_ctg{i}\n"));
            for chunk in seq.as_bytes().chunks(70) {
                out.push_str(std::str::from_utf8(chunk).unwrap());
                out.push('\n');
            }
        }
        fs::write(&p, out).expect("write genome");
        paths.push(p);
    }
    paths.sort();
    paths
}

/// Run the engine over a directory of genomes and return both the results and
/// the event stream, so a test can assert on either.
pub fn run_engine(
    inputs: &[PathBuf],
    params: &Parameters,
    outdir: Option<&Path>,
) -> (RunResults, Vec<Envelope>) {
    let (sink, rx) = event_bus();
    if let Some(dir) = outdir {
        sink.attach_log(dir).expect("attach event log");
    }
    let cfg = RunConfig {
        inputs: inputs.to_vec(),
        params: params.clone(),
        outdir: outdir.map(|p| p.to_path_buf()),
    };
    // The engine always needs a command receiver, even when no consumer will
    // ever send a command.
    let (_tx, cmd_rx) = crossbeam_channel::unbounded::<Command>();
    let results = run(cfg, sink, cmd_rx).expect("engine run failed");
    let mut events = Vec::new();
    while let Ok(e) = rx.try_recv() {
        events.push(e);
    }
    (results, events)
}

/// A deterministic pseudo-random DNA string. Not cryptographic, not biology:
/// just something with no repeats and no structure, so that every difference a
/// test sees is one the test put there.
pub fn random_dna(len: usize, seed: u64) -> String {
    let mut x = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let mut out = String::with_capacity(len);
    for _ in 0..len {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        out.push(b"ACGT"[(x >> 33) as usize & 3] as char);
    }
    out
}

/// Introduce `n` SNPs at evenly spaced, deterministic positions.
pub fn with_snps(seq: &str, n: usize, offset: usize) -> (String, Vec<(usize, char, char)>) {
    let bytes: Vec<u8> = seq.bytes().collect();
    let mut out = bytes.clone();
    let mut changes = Vec::new();
    let step = bytes.len() / (n + 1);
    for i in 0..n {
        let pos = offset + step * (i + 1);
        let original = out[pos];
        let replacement = *b"ACGT"
            .iter()
            .find(|&&b| b != original)
            .expect("there is always another base");
        out[pos] = replacement;
        changes.push((pos, original as char, replacement as char));
    }
    (String::from_utf8(out).unwrap(), changes)
}

/// Reverse complement a string.
pub fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect::<Vec<u8>>()
        .into_iter()
        .map(|b| b as char)
        .collect()
}

/// Read every output file of a run into a map, for byte-identity comparisons.
pub fn read_output_dir(dir: &Path) -> HashMap<String, Vec<u8>> {
    let mut out = HashMap::new();
    for entry in fs::read_dir(dir).expect("read output dir") {
        let p = entry.expect("dir entry").path();
        if p.is_file() {
            let name = p.file_name().unwrap().to_string_lossy().into_owned();
            out.insert(name, fs::read(&p).expect("read file"));
        }
    }
    out
}

/// Compare two output folders file by file, ignoring the files that legitimately
/// differ between runs (timestamps in `events.jsonl`, and `run.json` which
/// records wall-clock timings and peak memory).
pub fn diff_output_dirs(a: &Path, b: &Path) -> Vec<String> {
    const VOLATILE: &[&str] = &["events.jsonl", "run.json", "bench.json"];
    let (fa, fb) = (read_output_dir(a), read_output_dir(b));
    let mut diffs = Vec::new();
    let mut names: Vec<&String> = fa.keys().collect();
    names.sort();
    for name in names {
        if VOLATILE.contains(&name.as_str()) {
            continue;
        }
        match fb.get(name) {
            None => diffs.push(format!("{name}: only in {}", a.display())),
            Some(other) if other != &fa[name] => diffs.push(format!("{name}: content differs")),
            Some(_) => {}
        }
    }
    for name in fb.keys() {
        if !fa.contains_key(name) && !VOLATILE.contains(&name.as_str()) {
            diffs.push(format!("{name}: only in {}", b.display()));
        }
    }
    diffs.sort();
    diffs
}

/// The pairwise SNP counts of a finished run, keyed by sorted genome-name pair.
pub fn reported_pairs(results: &RunResults) -> HashMap<(String, String), u64> {
    let mut out = HashMap::new();
    for (_id, m) in results.matrices.iter() {
        let n = m.ids.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let a = &results.names[m.ids[i] as usize];
                let b = &results.names[m.ids[j] as usize];
                let key = if a <= b {
                    (a.clone(), b.clone())
                } else {
                    (b.clone(), a.clone())
                };
                out.insert(key, m.snps[i][j]);
            }
        }
    }
    out
}

/// The truth pairwise SNP counts written by the simulator.
pub fn truth_pairs(truth_dir: &Path) -> HashMap<(String, String), u64> {
    let txt =
        fs::read_to_string(truth_dir.join("pairwise_snps.tsv")).expect("truth pairwise_snps.tsv");
    let mut out = HashMap::new();
    for (i, line) in txt.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let (a, b) = if f[0] <= f[1] {
            (f[0].to_string(), f[1].to_string())
        } else {
            (f[1].to_string(), f[0].to_string())
        };
        out.insert((a, b), f[2].parse().expect("snp count"));
    }
    out
}
