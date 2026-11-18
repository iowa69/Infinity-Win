//! M-Tato's genome simulator: deterministic assemblies with known truth
//! (section 10.1 of the outline).
//!
//! The simulator exists so that the engine's accuracy can be *measured* rather
//! than argued about. It grows a pedigree of nearly identical bacterial genomes
//! from a reference, records every mutation it makes, and writes both the
//! genomes and the record out. Tests T1-T5 of section 10.3 and the accuracy
//! columns of section 11 all read the files this crate produces.
//!
//! # Quick start
//!
//! ```no_run
//! use mtato_sim::{simulate, RefSpec, SimConfig};
//!
//! let cfg = SimConfig {
//!     seed: 20260118,
//!     reference: RefSpec::Random { length: 2_000_000, gc: 0.5 },
//!     clusters: 1,
//!     leaves_per_cluster: 50,
//!     mean_snps_per_branch: 3.0,
//!     ..SimConfig::default()
//! };
//! let out = simulate(&cfg, std::path::Path::new("/tmp/t1")).unwrap();
//! assert_eq!(out.genomes.len(), 50);
//! ```
//!
//! # What it builds
//!
//! 1. **A reference.** Either a FASTA the user supplies, whose contigs are
//!    concatenated into one linear sequence, or an i.i.d. random sequence of a
//!    chosen length and GC fraction.
//! 2. **`m` cluster roots.** Each is the reference with `divergence` (default 2%)
//!    of its positions substituted, at positions drawn once and never reused.
//! 3. **A pedigree per cluster.** A random binary tree with `n` leaves. Every
//!    branch draws its SNP count from a Poisson distribution with mean
//!    `mean_snps_per_branch`, and that same number is both how many
//!    substitutions the branch applies and what the truth Newick says the branch
//!    length is. The two cannot drift apart because there is only one number.
//! 4. **Optional branch events.** Small indels of 1-10 bases, and gains and
//!    losses of 500-5000 base segments.
//! 5. **Optional fragmentation.** Each genome is cut into `contigs` pieces and
//!    the pieces are shuffled; with `random_orientation` each one is also
//!    reverse-complemented with probability one half. Both are off by default,
//!    so a plain run emits one forward contig per genome.
//!
//! # The invariant that makes the truth exact
//!
//! Inside a cluster, **no two events ever touch the same reference position**:
//! not two SNPs, not a SNP and a deletion, not two insertion anchors, and not
//! any event and any event of any other branch of the same cluster. A branch
//! SNP also never lands on a position an ancestor deleted, or on one of the
//! substitutions that made the cluster root.
//!
//! Two consequences, both of which the outline depends on:
//!
//! * The true SNP distance between two genomes is exactly the sum of the branch
//!   lengths on the path between them in the truth tree. Nothing cancels out,
//!   nothing is overwritten, so T1 can demand an exact match.
//! * A leaf's mutation log is a set of edits at pairwise distinct reference
//!   coordinates. Two leaves differ at exactly the union of their private SNP
//!   positions, which is the symmetric difference of their logs.
//!
//! Segment losses are the one place where truth and observation legitimately
//! diverge: a lost stretch is not a SNP difference, it is missing coverage, and
//! the substitutions the cluster root carried inside it disappear with it. Those
//! root substitutions are shared by every member of the cluster, so they cancel
//! out of every within-cluster distance anyway. This is why T2 and T4 allow a
//! pair to be off by up to two.
//!
//! # Determinism
//!
//! Same seed, same configuration, same reference file: byte-identical FASTA and
//! byte-identical truth files. Nothing is read from the clock, the generator is
//! the crate's own ([`rng`], xoshiro256\*\* seeded from SplitMix64) rather than
//! an external one whose version could change the output, no hash map is ever
//! iterated on the way to a file, and every derived stream is keyed by what it
//! is for. Two useful consequences:
//!
//! * Changing `contigs` or `random_orientation` changes how genomes are cut, not
//!   one base of what they contain.
//! * Changing the seed of one cluster's pedigree does not move another cluster.
//!
//! The one thing not fixed by the seed alone is the reference file: a different
//! FASTA gives a different dataset, as it must.
//!
//! # Output layout
//!
//! Written under the output directory given to [`simulate`]:
//!
//! ```text
//! reference.fasta                  the sequence the cluster roots were derived from
//! genomes/<genome>.fasta           one file per leaf genome
//! roots/<root>.fasta               one file per cluster root (an ancestor, not a leaf)
//! truth/manifest.json              format version, config and seed, reference, genome list
//! truth/clusters.tsv               genome, cluster_id, file, length, contigs
//! truth/pairwise_snps.tsv          genome_a, genome_b, snps
//! truth/tree_c<NN>.nwk             true tree of cluster NN, branch lengths in SNPs
//! truth/mutations_<name>.tsv       kind, position, ref, alt - one file per leaf and per root
//! ```
//!
//! Every TSV has a header row, uses `\n` line endings on every platform - which
//! is part of what makes a run byte-identical on Windows and on Linux - and ends
//! with a single newline. Names are built as `<name_prefix>_c<cluster>_s<leaf>`
//! with both numbers zero-padded, so name order and cluster/leaf order agree and
//! a lexicographic sort is also a numeric one. Cluster roots are `root_c<NN>`.
//!
//! A FASTA record is named after its genome when the file holds one contig, and
//! `<genome>_ctg<i>` (1-based) when it holds more.
//!
//! ## `clusters.tsv`
//!
//! | Column | Meaning |
//! |---|---|
//! | `genome` | Leaf genome name, also the FASTA file stem |
//! | `cluster_id` | 0-based cluster the genome belongs to |
//! | `file` | Its FASTA path, relative to the output directory |
//! | `length` | Total bases over all contigs |
//! | `contigs` | Number of records in that file |
//!
//! The first two columns are the assignment section 10.1 asks for; the rest lets
//! a script sanity-check a run without opening the FASTA files. The same
//! information is in `manifest.json`, and [`Truth::read_from`] refuses a
//! directory where the two disagree.
//!
//! ## `pairwise_snps.tsv`
//!
//! | Column | Meaning |
//! |---|---|
//! | `genome_a`, `genome_b` | Two leaves of the same cluster, `genome_a < genome_b` |
//! | `snps` | True SNP count between them: the path length in the truth tree |
//!
//! One row per unordered pair *inside* a cluster, rows sorted by the pair.
//! Cross-cluster pairs are deliberately absent: the engine never measures them
//! (section 6.4), so a truth value would only invite a comparison that cannot be
//! made. [`Truth::pairwise_snps`] returns `None` for them and `Some(0)` for a
//! genome against itself.
//!
//! ## `tree_c<NN>.nwk`
//!
//! Newick, written by [`mtato_core::tree::to_newick`], rooted at the cluster
//! root, with integer branch lengths in SNPs. `NN` is the cluster id padded to
//! at least two digits. The distance between two leaves is the sum of the branch
//! lengths on the path between them, which [`pairwise_from_newick`] recomputes
//! from the file alone - the independent check that truth test 3 is made of.
//!
//! ## `mutations_<name>.tsv`
//!
//! | Column | Meaning |
//! |---|---|
//! | `kind` | `SNP`, `INS`, `DEL`, `GAIN` or `LOSS` |
//! | `position` | 1-based reference coordinate, see below |
//! | `ref` | Bases before the event, or `-` when there were none |
//! | `alt` | Bases after the event, or `-` when there are none |
//!
//! For a leaf this is the log of the path from its cluster root to it; for a
//! `root_c<NN>` file it is the path from the reference to that root, which is all
//! substitutions. Rows are sorted by position.
//!
//! Positions are coordinates on the reference, the one sequence every genome of
//! a cluster descends from. Section 5's rule applies: 0-based in the program
//! ([`Mutation::pos`]), 1-based in the files. For `SNP`, `DEL` and `LOSS` the
//! position is the first affected base. For `INS` and `GAIN` it is the base
//! immediately *after* the inserted sequence, which is what anchors an insertion
//! to a coordinate it does not itself occupy.
//!
//! `ref` for a `LOSS` is the sequence as it stood in that lineage, so it carries
//! any substitution the cluster root introduced inside it - that is the string a
//! search of the root's FASTA will find. The length change an event caused is
//! `len(alt) - len(ref)`, so the total length of a genome is its root's length
//! plus the sum of the deltas in its log.
//!
//! Because of the exclusivity invariant above, the log is also a replay script:
//! applying its rows left to right, in reference coordinates, to the cluster
//! root's sequence reproduces the genome before fragmentation.
//!
//! ## `manifest.json`
//!
//! See [`Manifest`]. It carries [`FORMAT_VERSION`], the crate version, the whole
//! [`SimConfig`] with its seed, the reference's provenance and length, and one
//! entry per leaf genome and per cluster root with its file path, length and
//! contig count. Everything except the mutation rows, the Newick strings and the
//! pairwise table is in this file, so a dataset can be regenerated exactly from
//! it.
//!
//! # Consuming the truth from Rust
//!
//! [`simulate`] returns the truth in memory as well as on disk, so an
//! integration test does not have to read back what it just wrote:
//!
//! ```no_run
//! # use mtato_sim::{simulate, SimConfig};
//! let cfg = SimConfig::default();
//! let out = simulate(&cfg, std::path::Path::new("/tmp/x")).unwrap();
//! let names = out.truth.genome_names();
//! let d = out.truth.pairwise_snps(&names[0], &names[1]);
//! # let _ = d;
//! ```
//!
//! [`Truth::read_from`] is the exact inverse of [`Truth::write_to`], and
//! [`Truth`] compares equal to itself after a write-read round trip.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use mtato_core::fasta::{read_genome, Contig, GenomeData};
use mtato_core::genome::stats_of;
use mtato_core::seq::revcomp;

mod config;
mod pedigree;
mod positions;
mod rng;
mod spine;
mod treegen;
mod truth;

pub use config::{RefSpec, SegmentEvents, SimConfig};
pub use truth::{
    mutations_file_name, pairwise_from_newick, tree_file_name, GenomeTruth, Manifest, Mutation,
    MutationKind, RefSource, ReferenceTruth, RootTruth, Truth, FORMAT_VERSION, GENOMES_DIR,
    MUTATIONS_HEADER, ROOTS_DIR, TRUTH_DIR,
};

use pedigree::build_cluster;
use rng::{distinct_sorted, stream, Rng, TAG_FRAGMENT};

/// Record name of the generated reference, and the stem of `reference.fasta`.
const REFERENCE_NAME: &str = "reference";
/// File name of the generated reference, relative to the output directory.
const REFERENCE_FILE: &str = "reference.fasta";

/// What [`simulate`] produced.
#[derive(Debug, Clone, PartialEq)]
pub struct SimOutput {
    /// The output directory the files were written to.
    pub outdir: PathBuf,
    /// `reference.fasta`.
    pub reference: PathBuf,
    /// One FASTA per cluster root, ascending by cluster id. These are ancestors,
    /// not simulated genomes; they are written so a benchmark can measure a leaf
    /// against the sequence its mutation log is relative to.
    pub roots: Vec<PathBuf>,
    /// One FASTA per leaf genome, ascending by genome name. This is the list to
    /// hand to the engine.
    pub genomes: Vec<PathBuf>,
    /// The whole truth, in memory.
    pub truth: Truth,
}

/// Run a simulation and write everything it produces under `outdir`.
///
/// The directory is created if needed; existing files with the same names are
/// overwritten and stale files from an earlier, differently configured run are
/// left alone, so prefer an empty directory.
///
/// Errors are plain messages aimed at whoever typed the command line: a bad
/// configuration, an unreadable reference, a reference too short for the
/// requested contig count, or an I/O failure.
pub fn simulate(cfg: &SimConfig, outdir: &Path) -> Result<SimOutput, String> {
    cfg.validate()?;
    let reference = load_reference(cfg)?;
    let ref_seq: &[u8] = &reference.data.contigs[0].seq;
    if ref_seq.len() < cfg.contigs {
        return Err(format!(
            "the reference is {} bp and cannot be cut into {} contigs",
            ref_seq.len(),
            cfg.contigs
        ));
    }
    let novel_gc = cfg.gc.unwrap_or(reference.gc);

    for sub in ["", GENOMES_DIR, ROOTS_DIR, TRUTH_DIR] {
        let dir = if sub.is_empty() {
            outdir.to_path_buf()
        } else {
            outdir.join(sub)
        };
        fs::create_dir_all(&dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
    }
    let reference_path = outdir.join(REFERENCE_FILE);
    write_fasta(
        &reference_path,
        &reference.truth.name,
        &[ref_seq.to_vec()],
        cfg.fasta_width,
    )?;

    let cluster_width = digits(cfg.clusters - 1).max(2);
    let leaf_width = digits(cfg.leaves_per_cluster - 1).max(3);
    let mut writer = LeafWriter {
        cfg,
        outdir,
        genomes: Vec::with_capacity(cfg.genome_count()),
        genome_paths: Vec::with_capacity(cfg.genome_count()),
        roots: Vec::with_capacity(cfg.clusters),
        root_paths: Vec::with_capacity(cfg.clusters),
    };
    let mut trees = BTreeMap::new();
    let mut pairwise: Vec<(String, String, u64)> = Vec::new();

    for cluster in 0..cfg.clusters as u32 {
        let leaf_names: Vec<String> = (0..cfg.leaves_per_cluster)
            .map(|leaf| leaf_name(cfg, cluster, leaf, cluster_width, leaf_width))
            .collect();
        let build = build_cluster(
            cluster,
            cfg,
            ref_seq,
            &leaf_names,
            novel_gc,
            &mut |leaf, name, seq, mutations| writer.emit(cluster, leaf, name, seq, mutations),
        )?;

        let root_name = root_name(cluster, cluster_width);
        let root_file = format!("{ROOTS_DIR}/{root_name}.fasta");
        let root_path = outdir.join(&root_file);
        let root_length = build.root_seq.len() as u64;
        write_fasta(&root_path, &root_name, &[build.root_seq], cfg.fasta_width)?;
        writer.root_paths.push(root_path);
        writer.roots.push(RootTruth {
            cluster,
            name: root_name,
            file: root_file,
            length: root_length,
            substitutions: build.root_mutations.len() as u64,
            mutations: build.root_mutations,
        });

        trees.insert(cluster, build.pedigree.to_newick(&leaf_names));
        let index = build.pedigree.distance_index();
        let leaf_nodes = build.pedigree.leaf_nodes();
        for i in 0..leaf_names.len() {
            for j in (i + 1)..leaf_names.len() {
                let snps = index.distance(leaf_nodes[i], leaf_nodes[j]);
                pairwise.push((leaf_names[i].clone(), leaf_names[j].clone(), snps));
            }
        }
    }

    // The pedigree walk emits leaves in tree order; the truth is keyed by name.
    let genomes = std::mem::take(&mut writer.genomes);
    let paths = std::mem::take(&mut writer.genome_paths);
    let mut emitted: Vec<(GenomeTruth, PathBuf)> = genomes.into_iter().zip(paths).collect();
    emitted.sort_by(|a, b| a.0.name.cmp(&b.0.name));
    let (genomes, genome_paths): (Vec<GenomeTruth>, Vec<PathBuf>) = emitted.into_iter().unzip();
    pairwise.sort();

    let truth = Truth {
        config: cfg.clone(),
        reference: reference.truth,
        roots: writer.roots,
        genomes,
        trees,
        pairwise,
    };
    truth.write_to(outdir)?;

    Ok(SimOutput {
        outdir: outdir.to_path_buf(),
        reference: reference_path,
        roots: writer.root_paths,
        genomes: genome_paths,
        truth,
    })
}

/// Collects leaf genomes as the pedigree walk produces them, so that only one
/// rendered genome per cluster is ever in memory.
struct LeafWriter<'a> {
    cfg: &'a SimConfig,
    outdir: &'a Path,
    genomes: Vec<GenomeTruth>,
    genome_paths: Vec<PathBuf>,
    roots: Vec<RootTruth>,
    root_paths: Vec<PathBuf>,
}

impl LeafWriter<'_> {
    /// Fragment one leaf, write its FASTA and record its truth.
    fn emit(
        &mut self,
        cluster: u32,
        leaf: usize,
        name: &str,
        seq: Vec<u8>,
        mutations: Vec<Mutation>,
    ) -> Result<(), String> {
        // Keyed by (cluster, leaf), not by emission order, so the way a genome is
        // cut depends only on which genome it is.
        let key = u64::from(cluster) * self.cfg.leaves_per_cluster as u64 + leaf as u64;
        let mut rng = stream(self.cfg.seed, TAG_FRAGMENT, key);
        let contigs = fragment(seq, self.cfg.contigs, self.cfg.random_orientation, &mut rng);
        let file = format!("{GENOMES_DIR}/{name}.fasta");
        let path = self.outdir.join(&file);
        write_fasta(&path, name, &contigs, self.cfg.fasta_width)?;
        let length = contigs.iter().map(|c| c.len() as u64).sum();
        self.genomes.push(GenomeTruth {
            name: name.to_string(),
            cluster,
            file,
            length,
            contigs: contigs.len() as u32,
            mutations,
        });
        self.genome_paths.push(path);
        Ok(())
    }
}

/// The reference sequence, its measured GC and its truth record.
struct LoadedReference {
    /// Always a single contig: a supplied FASTA is concatenated, because every
    /// simulated genome is cut into [`SimConfig::contigs`] contigs anyway and a
    /// second coordinate system would only make the truth harder to read.
    data: GenomeData,
    gc: f64,
    truth: ReferenceTruth,
}

fn load_reference(cfg: &SimConfig) -> Result<LoadedReference, String> {
    match &cfg.reference {
        RefSpec::Random { length, gc } => {
            let seq = pedigree::random_reference(cfg.seed, *length, *gc);
            finish_reference(seq, REFERENCE_NAME.to_string(), RefSource::Random, None)
        }
        RefSpec::File(path) => {
            let data = read_genome(path)
                .map_err(|e| format!("cannot read reference {}: {e}", path.display()))?;
            let mut seq: Vec<u8> = Vec::with_capacity(data.total_length() as usize);
            for contig in &data.contigs {
                seq.extend_from_slice(&contig.seq);
            }
            if seq.is_empty() {
                return Err(format!("reference {} contains no sequence", path.display()));
            }
            let name = data.name.clone();
            finish_reference(seq, name, RefSource::File, Some(path.clone()))
        }
    }
}

fn finish_reference(
    seq: Vec<u8>,
    name: String,
    source: RefSource,
    source_path: Option<PathBuf>,
) -> Result<LoadedReference, String> {
    let length = seq.len() as u64;
    let data = GenomeData {
        name: name.clone(),
        contigs: vec![Contig {
            name: name.clone(),
            seq,
        }],
    };
    let gc = stats_of(&data).gc;
    Ok(LoadedReference {
        data,
        gc,
        truth: ReferenceTruth {
            name,
            file: REFERENCE_FILE.to_string(),
            length,
            gc,
            source,
            source_path,
        },
    })
}

/// Cut a genome into `count` contigs, in random order and random orientation.
///
/// A genome shorter than `count` bases is cut into as many contigs as it has
/// bases, since a contig cannot be empty. With `count == 1` and no orientation
/// randomness the sequence is returned untouched and the generator is not
/// consumed at all.
fn fragment(seq: Vec<u8>, count: usize, orient: bool, rng: &mut Rng) -> Vec<Vec<u8>> {
    let len = seq.len();
    let n = count.min(len).max(1);
    if n == 1 && !orient {
        return vec![seq];
    }
    let cuts = if n > 1 {
        distinct_sorted(rng, 1, len, n - 1)
    } else {
        Vec::new()
    };
    let mut pieces: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut start = 0usize;
    for &cut in &cuts {
        pieces.push(seq[start..cut].to_vec());
        start = cut;
    }
    pieces.push(seq[start..].to_vec());
    rng.shuffle(&mut pieces);
    if orient {
        for piece in pieces.iter_mut() {
            if rng.coin() {
                *piece = revcomp(piece);
            }
        }
    }
    pieces
}

/// Write a multi-record FASTA. Line endings are `\n` on every platform, which
/// is what makes "byte-identical for a given seed" hold on Windows too.
fn write_fasta(path: &Path, name: &str, contigs: &[Vec<u8>], width: usize) -> Result<(), String> {
    let width = width.max(1);
    let mut buf: Vec<u8> = Vec::new();
    for (i, seq) in contigs.iter().enumerate() {
        buf.push(b'>');
        buf.extend_from_slice(name.as_bytes());
        if contigs.len() > 1 {
            buf.extend_from_slice(format!("_ctg{}", i + 1).as_bytes());
        }
        buf.push(b'\n');
        if seq.is_empty() {
            buf.push(b'\n');
            continue;
        }
        for line in seq.chunks(width) {
            buf.extend_from_slice(line);
            buf.push(b'\n');
        }
    }
    fs::write(path, &buf).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

/// Name of a leaf genome. Both numbers are zero-padded to a fixed width so that
/// sorting names as strings also sorts them as numbers.
fn leaf_name(
    cfg: &SimConfig,
    cluster: u32,
    leaf: usize,
    cluster_width: usize,
    leaf_width: usize,
) -> String {
    format!(
        "{}_c{:0cw$}_s{:0lw$}",
        cfg.name_prefix,
        cluster,
        leaf,
        cw = cluster_width,
        lw = leaf_width
    )
}

/// Name of a cluster root.
fn root_name(cluster: u32, cluster_width: usize) -> String {
    format!("root_c{:0cw$}", cluster, cw = cluster_width)
}

/// How many digits `v` needs, at least 1.
fn digits(v: usize) -> usize {
    v.to_string().len()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(tag: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("mtato-sim-lib-{}-{}", std::process::id(), tag));
        let _ = fs::remove_dir_all(&p);
        p
    }

    fn small_cfg(seed: u64) -> SimConfig {
        SimConfig {
            seed,
            reference: RefSpec::Random {
                length: 20_000,
                gc: 0.5,
            },
            clusters: 1,
            leaves_per_cluster: 4,
            mean_snps_per_branch: 3.0,
            ..SimConfig::default()
        }
    }

    #[test]
    fn fragment_keeps_every_base_and_the_piece_count() {
        let seq: Vec<u8> = (0..1000).map(|i| b"ACGT"[i % 4]).collect();
        let mut rng = Rng::from_seed(1);
        let pieces = fragment(seq.clone(), 7, true, &mut rng);
        assert_eq!(pieces.len(), 7);
        assert_eq!(pieces.iter().map(|p| p.len()).sum::<usize>(), 1000);
        assert!(pieces.iter().all(|p| !p.is_empty()));
        // Every piece is a slice of the genome or of its reverse complement.
        let rc = revcomp(&seq);
        for p in &pieces {
            let forward = seq.windows(p.len()).any(|w| w == p.as_slice());
            let backward = rc.windows(p.len()).any(|w| w == p.as_slice());
            assert!(
                forward || backward,
                "piece of {} bases is not from the genome",
                p.len()
            );
        }
    }

    #[test]
    fn fragment_without_orientation_keeps_every_piece_forward() {
        // A random sequence, not a periodic one: with "ACGTACGT..." a piece
        // would occur at many offsets and the reassembly below would be ambiguous.
        let seq = pedigree::random_reference(21, 100, 0.5);
        let mut rng = Rng::from_seed(2);
        let pieces = fragment(seq.clone(), 4, false, &mut rng);
        assert_eq!(pieces.len(), 4);
        // Shuffled, so the order is not the genome's, but no piece is flipped.
        for p in &pieces {
            assert!(seq.windows(p.len()).any(|w| w == p.as_slice()));
        }
        let mut rejoined: Vec<Vec<u8>> = pieces;
        rejoined.sort_by_key(|p| {
            seq.windows(p.len())
                .position(|w| w == p.as_slice())
                .unwrap()
        });
        let flat: Vec<u8> = rejoined.into_iter().flatten().collect();
        assert_eq!(flat, seq, "put back in order, the pieces are the genome");
    }

    #[test]
    fn fragment_is_deterministic_for_a_seed() {
        let seq: Vec<u8> = (0..500).map(|i| b"ACGT"[i % 4]).collect();
        let a = fragment(seq.clone(), 10, true, &mut Rng::from_seed(3));
        let b = fragment(seq.clone(), 10, true, &mut Rng::from_seed(3));
        let c = fragment(seq.clone(), 10, true, &mut Rng::from_seed(4));
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn fragment_caps_at_one_contig_per_base() {
        let seq = b"ACGT".to_vec();
        let mut rng = Rng::from_seed(5);
        let pieces = fragment(seq, 99, false, &mut rng);
        assert_eq!(pieces.len(), 4);
        assert_eq!(pieces.iter().map(|p| p.len()).sum::<usize>(), 4);
    }

    #[test]
    fn fragment_of_one_contig_consumes_nothing() {
        let seq = b"ACGTACGT".to_vec();
        let mut rng = Rng::from_seed(6);
        assert_eq!(fragment(seq.clone(), 1, false, &mut rng), vec![seq.clone()]);
        assert_eq!(rng.next_u64(), Rng::from_seed(6).next_u64());
    }

    #[test]
    fn fasta_records_are_named_after_the_genome_only_when_single() {
        let dir = tmpdir("fasta");
        fs::create_dir_all(&dir).unwrap();
        let one = dir.join("one.fasta");
        write_fasta(&one, "g1", &[vec![b'A'; 130]], 60).unwrap();
        let text = fs::read_to_string(&one).unwrap();
        assert!(text.starts_with(">g1\n"));
        assert_eq!(text.lines().count(), 1 + 3);
        assert!(text.ends_with('\n'));
        assert!(!text.contains("\r"));

        let many = dir.join("many.fasta");
        write_fasta(&many, "g2", &[vec![b'C'; 10], vec![b'G'; 5]], 4).unwrap();
        let text = fs::read_to_string(&many).unwrap();
        assert_eq!(text, ">g2_ctg1\nCCCC\nCCCC\nCC\n>g2_ctg2\nGGGG\nG\n");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn names_sort_numerically() {
        let cfg = SimConfig {
            clusters: 12,
            leaves_per_cluster: 1200,
            ..SimConfig::default()
        };
        let cw = digits(cfg.clusters - 1).max(2);
        let lw = digits(cfg.leaves_per_cluster - 1).max(3);
        let names: Vec<String> = (0..cfg.clusters as u32)
            .flat_map(|c| (0..cfg.leaves_per_cluster).map(move |l| (c, l)))
            .map(|(c, l)| leaf_name(&cfg, c, l, cw, lw))
            .collect();
        let mut sorted = names.clone();
        sorted.sort();
        assert_eq!(names, sorted, "names must already be in ascending order");
        assert_eq!(names[0], "sim_c00_s0000");
        assert_eq!(names[1199], "sim_c00_s1199");
        assert_eq!(root_name(3, cw), "root_c03");
    }

    #[test]
    fn a_written_genome_reads_back_through_the_engine() {
        let dir = tmpdir("readback");
        let cfg = small_cfg(17);
        let out = simulate(&cfg, &dir).unwrap();
        assert_eq!(out.genomes.len(), 4);
        for path in &out.genomes {
            let data = read_genome(path).unwrap();
            assert_eq!(data.contigs.len(), 1);
            let stats = stats_of(&data);
            assert_eq!(stats.length, 20_000);
            assert_eq!(stats.n_bases, 0);
            assert!((stats.gc - 0.5).abs() < 0.05);
        }
        // Names derived from the paths match the truth.
        let from_paths: Vec<String> = out
            .genomes
            .iter()
            .map(|p| mtato_core::fasta::genome_name(p))
            .collect();
        assert_eq!(from_paths, out.truth.genome_names());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_supplied_reference_is_concatenated_and_reported() {
        let dir = tmpdir("ref-file");
        fs::create_dir_all(&dir).unwrap();
        let ref_path = dir.join("myref.fasta");
        fs::write(
            &ref_path,
            b">c1 first\nACGTACGTAC\n>c2 second\nGGGGCCCCAAAA\n",
        )
        .unwrap();
        let out_dir = dir.join("out");
        let cfg = SimConfig {
            seed: 3,
            reference: RefSpec::File(ref_path.clone()),
            divergence: 0.0,
            leaves_per_cluster: 2,
            mean_snps_per_branch: 0.0,
            // Otherwise the single contig may come back reverse-complemented.
            random_orientation: false,
            ..SimConfig::default()
        };
        let out = simulate(&cfg, &out_dir).unwrap();
        assert_eq!(out.truth.reference.length, 22);
        assert_eq!(out.truth.reference.name, "myref");
        assert_eq!(out.truth.reference.source, RefSource::File);
        assert_eq!(
            out.truth.reference.source_path.as_deref(),
            Some(ref_path.as_path())
        );
        // Zero divergence and zero SNPs: every genome is the reference itself.
        let data = read_genome(&out.genomes[0]).unwrap();
        assert_eq!(data.contigs[0].seq, b"ACGTACGTACGGGGCCCCAAAA");
        assert_eq!(
            out.truth
                .pairwise_snps(&out.truth.genome_names()[0], &out.truth.genome_names()[1]),
            Some(0)
        );
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_missing_reference_is_reported_clearly() {
        let dir = tmpdir("missing-ref");
        let cfg = SimConfig {
            reference: RefSpec::File(PathBuf::from("/nonexistent/ref.fasta")),
            leaves_per_cluster: 2,
            ..SimConfig::default()
        };
        let err = simulate(&cfg, &dir).unwrap_err();
        assert!(err.contains("cannot read reference"), "{err}");
    }

    #[test]
    fn a_reference_shorter_than_the_contig_count_is_refused() {
        let dir = tmpdir("short-ref");
        let cfg = SimConfig {
            reference: RefSpec::Random {
                length: 40,
                gc: 0.5,
            },
            contigs: 3,
            leaves_per_cluster: 2,
            ..SimConfig::default()
        };
        // 40 bases is enough for 3 contigs, so this one succeeds.
        assert!(simulate(&cfg, &dir).is_ok());
        let too_many = SimConfig {
            contigs: 100,
            ..cfg
        };
        let err = simulate(&too_many, &dir.join("b")).unwrap_err();
        assert!(err.contains("cannot be cut into"), "{err}");
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn an_invalid_config_is_refused_before_any_work() {
        let dir = tmpdir("invalid");
        let cfg = SimConfig {
            divergence: 1.5,
            ..SimConfig::default()
        };
        assert!(simulate(&cfg, &dir).is_err());
        assert!(!dir.exists());
    }
}
