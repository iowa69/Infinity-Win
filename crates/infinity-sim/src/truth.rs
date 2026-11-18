//! The truth a simulation leaves behind, and its on-disk formats.
//!
//! # Layout
//!
//! Everything below is relative to the simulation output directory, the one
//! passed to [`crate::simulate`]:
//!
//! ```text
//! reference.fasta                  the sequence the cluster roots were derived from
//! genomes/<genome>.fasta           one file per leaf genome
//! roots/<root>.fasta               one file per cluster root (not a leaf genome)
//! truth/manifest.json              format version, config, seed, reference, genome list
//! truth/clusters.tsv               genome -> cluster, plus file, length and contigs
//! truth/pairwise_snps.tsv          true SNP count for every pair inside a cluster
//! truth/tree_c<NN>.nwk             true tree of cluster NN, branch lengths in SNPs
//! truth/mutations_<name>.tsv       event log of one leaf, or of one cluster root
//! ```
//!
//! # Formats
//!
//! All TSV files have a header row, are UTF-8, use `\n` line endings on every
//! platform - which is what makes a run byte-identical on Windows and on Linux -
//! and end with a single newline.
//!
//! * `clusters.tsv` - `genome`, `cluster_id`, `file`, `length`, `contigs`. The
//!   first two columns are the assignment section 10.1 asks for; the rest is
//!   there so a script never has to open the FASTA files to sanity-check a run.
//! * `pairwise_snps.tsv` - `genome_a`, `genome_b`, `snps`. One row per unordered
//!   pair of leaves *of the same cluster*, `genome_a < genome_b` in byte order,
//!   rows sorted by that pair. Cross-cluster pairs are absent on purpose: the
//!   engine never measures them (section 6.4), so a truth value for them would
//!   only invite a comparison that cannot be made.
//! * `tree_c<NN>.nwk` - Newick with integer branch lengths in SNPs, written by
//!   [`infinity_core::tree::to_newick`]. `NN` is the cluster id padded to at least
//!   two digits. The tree is rooted at the cluster root, so the distance between
//!   two leaves is the sum of the branch lengths on the path between them.
//! * `mutations_<name>.tsv` - `kind`, `position`, `ref`, `alt`. One row per event
//!   on the path from the cluster root to that leaf (or, for a `root_c<NN>` file,
//!   from the reference to that root). Rows are sorted by position. An empty
//!   `ref` or `alt` is written as `-`.
//! * `manifest.json` - see [`Manifest`]. It carries the whole [`SimConfig`]
//!   including the seed, so a dataset can be regenerated exactly.
//!
//! # Positions
//!
//! Following section 5, positions are 0-based inside the program and 1-based in
//! the files; [`Mutation::pos`] is 0-based and the conversion happens once, in
//! the writers and readers below. All of them are coordinates on the *reference*,
//! which is the one sequence every genome of a cluster descends from.
//!
//! For `SNP`, `DEL` and `LOSS`, `position` is the first affected reference base.
//! For `INS` and `GAIN`, it is the reference base immediately *after* the
//! inserted sequence, which is what anchors an insertion to a coordinate it does
//! not itself occupy.
//!
//! # The exclusivity invariant
//!
//! Inside a cluster, no two events on the tree ever touch the same reference
//! position: not two SNPs, not a SNP and a deletion, not two insertions. That is
//! what makes the truth honest - a pairwise SNP count is then exactly the number
//! of positions where the two genomes differ, which is exactly the path length in
//! the tree, with nothing cancelling out. The occupancy bookkeeping inside the
//! simulator enforces it, and truth test 2 asserts it from these very files.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::config::SimConfig;

/// Version of the truth layout. Bumped when a file's columns change.
pub const FORMAT_VERSION: u32 = 1;

/// Name of the truth subdirectory inside the output directory.
pub const TRUTH_DIR: &str = "truth";
/// Name of the leaf genome subdirectory.
pub const GENOMES_DIR: &str = "genomes";
/// Name of the cluster root subdirectory.
pub const ROOTS_DIR: &str = "roots";

/// What kind of event a [`Mutation`] records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum MutationKind {
    /// One substituted base.
    #[serde(rename = "SNP")]
    Snp,
    /// A short inserted segment, 1-10 bases by default.
    #[serde(rename = "INS")]
    Insertion,
    /// A short deleted segment, 1-10 bases by default.
    #[serde(rename = "DEL")]
    Deletion,
    /// A gained segment, 500-5000 bases by default.
    #[serde(rename = "GAIN")]
    Gain,
    /// A lost segment, 500-5000 bases by default.
    #[serde(rename = "LOSS")]
    Loss,
}

impl MutationKind {
    /// The label used in the TSV files.
    pub fn as_str(self) -> &'static str {
        match self {
            MutationKind::Snp => "SNP",
            MutationKind::Insertion => "INS",
            MutationKind::Deletion => "DEL",
            MutationKind::Gain => "GAIN",
            MutationKind::Loss => "LOSS",
        }
    }

    /// Parse a TSV label.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "SNP" => Ok(MutationKind::Snp),
            "INS" => Ok(MutationKind::Insertion),
            "DEL" => Ok(MutationKind::Deletion),
            "GAIN" => Ok(MutationKind::Gain),
            "LOSS" => Ok(MutationKind::Loss),
            other => Err(format!("unknown mutation kind {other:?}")),
        }
    }

    /// True for the two events that add sequence.
    pub fn is_insertion(self) -> bool {
        matches!(self, MutationKind::Insertion | MutationKind::Gain)
    }

    /// True for the two events that remove sequence.
    pub fn is_deletion(self) -> bool {
        matches!(self, MutationKind::Deletion | MutationKind::Loss)
    }
}

/// One event in a genome's history.
///
/// `ref_bases` is what was there before and `alt_bases` what is there after, so
/// the length change an event caused is `alt_bases.len() - ref_bases.len()`. For
/// a `LOSS` the recorded `ref_bases` is the sequence as it was in that lineage,
/// which carries any substitution an ancestor introduced inside it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Mutation {
    /// What happened.
    pub kind: MutationKind,
    /// 0-based reference coordinate, see the module docs.
    pub pos: usize,
    /// Bases before the event; empty for an insertion or a gain.
    pub ref_bases: Vec<u8>,
    /// Bases after the event; empty for a deletion or a loss.
    pub alt_bases: Vec<u8>,
}

impl Mutation {
    /// A substitution at a 0-based reference position.
    pub fn snp(pos: usize, from: u8, to: u8) -> Self {
        Self {
            kind: MutationKind::Snp,
            pos,
            ref_bases: vec![from],
            alt_bases: vec![to],
        }
    }

    /// Change in genome length this event caused, in bases.
    pub fn length_delta(&self) -> i64 {
        self.alt_bases.len() as i64 - self.ref_bases.len() as i64
    }

    /// One TSV row, positions converted to 1-based. No header, no newline.
    pub fn to_tsv_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}",
            self.kind.as_str(),
            self.pos + 1,
            bases_field(&self.ref_bases),
            bases_field(&self.alt_bases)
        )
    }

    /// Parse one TSV row written by [`Mutation::to_tsv_row`].
    pub fn from_tsv_row(line: &str) -> Result<Self, String> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 4 {
            return Err(format!(
                "mutation row needs 4 columns, got {}: {line:?}",
                f.len()
            ));
        }
        let pos = f[1]
            .parse::<usize>()
            .map_err(|e| format!("bad mutation position {:?}: {e}", f[1]))?;
        if pos == 0 {
            return Err(format!(
                "mutation position must be 1-based, got 0 in {line:?}"
            ));
        }
        Ok(Self {
            kind: MutationKind::parse(f[0])?,
            pos: pos - 1,
            ref_bases: parse_bases_field(f[2]),
            alt_bases: parse_bases_field(f[3]),
        })
    }
}

/// The TSV header of a mutations file.
pub const MUTATIONS_HEADER: &str = "kind\tposition\tref\talt";

fn bases_field(b: &[u8]) -> String {
    if b.is_empty() {
        "-".to_string()
    } else {
        String::from_utf8_lossy(b).into_owned()
    }
}

fn parse_bases_field(s: &str) -> Vec<u8> {
    if s == "-" {
        Vec::new()
    } else {
        s.as_bytes().to_vec()
    }
}

/// Where the reference sequence came from.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RefSource {
    /// Drawn by the simulator.
    Random,
    /// Read from a user-supplied FASTA.
    File,
}

/// The reference every cluster root was derived from.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceTruth {
    /// Record name used in `reference.fasta`.
    pub name: String,
    /// Path relative to the output directory.
    pub file: String,
    /// Total bases.
    pub length: u64,
    /// Measured GC fraction.
    pub gc: f64,
    /// Drawn or supplied.
    pub source: RefSource,
    /// The file it was read from, when it was supplied.
    pub source_path: Option<PathBuf>,
}

/// One cluster root: the common ancestor of a cluster's leaves.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RootTruth {
    /// Cluster id, 0-based.
    pub cluster: u32,
    /// Genome-style name, `<prefix-less>root_c<NN>`.
    pub name: String,
    /// Path of the root's FASTA, relative to the output directory.
    pub file: String,
    /// Total bases.
    pub length: u64,
    /// Number of substitutions that separate it from the reference.
    pub substitutions: u64,
    /// The substitutions themselves, sorted by position.
    pub mutations: Vec<Mutation>,
}

/// One simulated leaf genome.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenomeTruth {
    /// Genome name, also the FASTA file stem.
    pub name: String,
    /// Cluster it belongs to, 0-based.
    pub cluster: u32,
    /// Path of its FASTA, relative to the output directory.
    pub file: String,
    /// Total bases over all contigs.
    pub length: u64,
    /// Number of contigs.
    pub contigs: u32,
    /// Events on the path from its cluster root to this genome, sorted by
    /// position. Empty for a genome identical to its root.
    pub mutations: Vec<Mutation>,
}

/// The `manifest.json` schema.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    /// Layout version, [`FORMAT_VERSION`].
    pub format_version: u32,
    /// Version of the crate that wrote it.
    pub infinity_version: String,
    /// The full configuration, seed included.
    pub config: SimConfig,
    /// The reference sequence.
    pub reference: ReferenceTruth,
    /// Cluster roots, ascending by cluster id.
    pub roots: Vec<RootTruth>,
    /// Leaf genomes, ascending by name.
    pub genomes: Vec<GenomeTruth>,
    /// Number of rows in `pairwise_snps.tsv`.
    pub pairwise_pairs: usize,
    /// Cluster ids that have a `tree_c<NN>.nwk`.
    pub clusters: Vec<u32>,
}

/// Everything a simulation knows that the engine does not.
///
/// Write it with [`Truth::write_to`] and read it back with [`Truth::read_from`];
/// the two are exact inverses, so a Rust integration test can take the truth
/// straight from [`crate::SimOutput`] without going through the files, and a
/// Python benchmark can take it from the files without going through Rust.
#[derive(Debug, Clone, PartialEq)]
pub struct Truth {
    /// The configuration that produced this truth, seed included.
    pub config: SimConfig,
    /// The reference sequence.
    pub reference: ReferenceTruth,
    /// One entry per cluster, ascending by cluster id.
    pub roots: Vec<RootTruth>,
    /// One entry per leaf genome, ascending by name.
    pub genomes: Vec<GenomeTruth>,
    /// Newick of each cluster's true tree, keyed by cluster id. Branch lengths
    /// are SNP counts.
    pub trees: BTreeMap<u32, String>,
    /// True SNP count for every within-cluster pair, as `(genome_a, genome_b,
    /// snps)` with `genome_a < genome_b`, sorted by that pair.
    pub pairwise: Vec<(String, String, u64)>,
}

/// File name of a cluster's truth tree, relative to the truth directory.
pub fn tree_file_name(cluster: u32) -> String {
    format!("tree_c{cluster:02}.nwk")
}

/// File name of a genome's or root's mutation log, relative to the truth
/// directory.
pub fn mutations_file_name(name: &str) -> String {
    format!("mutations_{name}.tsv")
}

impl Truth {
    /// Write every truth file under `outdir/truth`.
    ///
    /// `outdir` is the simulation output directory, not the truth directory
    /// itself: the FASTA paths recorded in the manifest are relative to it. The
    /// genome and reference FASTA files are *not* written here - they are large
    /// and [`crate::simulate`] streams them out as it goes.
    pub fn write_to(&self, outdir: &Path) -> Result<(), String> {
        let dir = outdir.join(TRUTH_DIR);
        fs::create_dir_all(&dir).map_err(|e| format!("cannot create {}: {e}", dir.display()))?;

        let manifest = Manifest {
            format_version: FORMAT_VERSION,
            infinity_version: infinity_core::VERSION.to_string(),
            config: self.config.clone(),
            reference: self.reference.clone(),
            roots: self.roots.clone(),
            genomes: self.genomes.clone(),
            pairwise_pairs: self.pairwise.len(),
            clusters: self.trees.keys().copied().collect(),
        };
        let json = serde_json::to_string_pretty(&manifest)
            .map_err(|e| format!("cannot serialise the manifest: {e}"))?;
        write_text(&dir.join("manifest.json"), &json)?;

        let mut clusters_tsv = String::from("genome\tcluster_id\tfile\tlength\tcontigs");
        for g in &self.genomes {
            let _ = write!(
                clusters_tsv,
                "\n{}\t{}\t{}\t{}\t{}",
                g.name, g.cluster, g.file, g.length, g.contigs
            );
        }
        clusters_tsv.push('\n');
        write_text(&dir.join("clusters.tsv"), &clusters_tsv)?;

        let mut pairs = String::from("genome_a\tgenome_b\tsnps");
        for (a, b, s) in &self.pairwise {
            let _ = write!(pairs, "\n{a}\t{b}\t{s}");
        }
        pairs.push('\n');
        write_text(&dir.join("pairwise_snps.tsv"), &pairs)?;

        for (&cluster, newick) in &self.trees {
            write_text(&dir.join(tree_file_name(cluster)), newick)?;
        }
        for g in &self.genomes {
            write_mutations(&dir, &g.name, &g.mutations)?;
        }
        for r in &self.roots {
            write_mutations(&dir, &r.name, &r.mutations)?;
        }
        Ok(())
    }

    /// Read a truth directory written by [`Truth::write_to`].
    ///
    /// `outdir` is the simulation output directory. Fails when a file is
    /// missing, when the format version is not the one this crate writes, or
    /// when `clusters.tsv` disagrees with the manifest.
    pub fn read_from(outdir: &Path) -> Result<Self, String> {
        let dir = outdir.join(TRUTH_DIR);
        let json = read_text(&dir.join("manifest.json"))?;
        let manifest: Manifest = serde_json::from_str(&json)
            .map_err(|e| format!("cannot parse {}: {e}", dir.join("manifest.json").display()))?;
        if manifest.format_version != FORMAT_VERSION {
            return Err(format!(
                "truth format version {} is not supported by this build (writes {})",
                manifest.format_version, FORMAT_VERSION
            ));
        }

        let pairs_tsv = read_text(&dir.join("pairwise_snps.tsv"))?;
        let pairwise = parse_pairwise(&pairs_tsv)?;
        if pairwise.len() != manifest.pairwise_pairs {
            return Err(format!(
                "pairwise_snps.tsv has {} rows but the manifest announces {}",
                pairwise.len(),
                manifest.pairwise_pairs
            ));
        }

        let mut trees = BTreeMap::new();
        for &cluster in &manifest.clusters {
            let path = dir.join(tree_file_name(cluster));
            trees.insert(cluster, read_text(&path)?);
        }

        let mut genomes = manifest.genomes.clone();
        for g in &mut genomes {
            g.mutations = read_mutations(&dir, &g.name)?;
        }
        let mut roots = manifest.roots.clone();
        for r in &mut roots {
            r.mutations = read_mutations(&dir, &r.name)?;
        }

        let truth = Self {
            config: manifest.config,
            reference: manifest.reference,
            roots,
            genomes,
            trees,
            pairwise,
        };
        truth.check_clusters_tsv(&dir)?;
        Ok(truth)
    }

    /// `clusters.tsv` duplicates what the manifest already says, so a reader
    /// that trusts one of them can be wrong. Checking costs one small file and
    /// catches a truncated or hand-edited run.
    fn check_clusters_tsv(&self, dir: &Path) -> Result<(), String> {
        let text = read_text(&dir.join("clusters.tsv"))?;
        let mut lines = text.lines();
        let header = lines.next().unwrap_or("");
        if header != "genome\tcluster_id\tfile\tlength\tcontigs" {
            return Err(format!("clusters.tsv has an unexpected header {header:?}"));
        }
        let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
        if rows.len() != self.genomes.len() {
            return Err(format!(
                "clusters.tsv has {} rows, the manifest announces {}",
                rows.len(),
                self.genomes.len()
            ));
        }
        for (row, g) in rows.iter().zip(&self.genomes) {
            let ok = row.len() == 5
                && row[0] == g.name
                && row[1] == g.cluster.to_string()
                && row[2] == g.file
                && row[3] == g.length.to_string()
                && row[4] == g.contigs.to_string();
            if !ok {
                return Err(format!(
                    "clusters.tsv row {row:?} does not match the manifest"
                ));
            }
        }
        Ok(())
    }

    /// The cluster a genome belongs to.
    pub fn cluster_of(&self, genome: &str) -> Option<u32> {
        self.genomes
            .iter()
            .find(|g| g.name == genome)
            .map(|g| g.cluster)
    }

    /// Genome name to cluster id, for consumers that want the map directly.
    pub fn clusters(&self) -> BTreeMap<String, u32> {
        self.genomes
            .iter()
            .map(|g| (g.name.clone(), g.cluster))
            .collect()
    }

    /// Cluster ids present in the truth, ascending.
    pub fn cluster_ids(&self) -> Vec<u32> {
        let mut ids: Vec<u32> = self.genomes.iter().map(|g| g.cluster).collect();
        ids.sort_unstable();
        ids.dedup();
        ids
    }

    /// The names of every leaf genome, ascending.
    pub fn genome_names(&self) -> Vec<String> {
        self.genomes.iter().map(|g| g.name.clone()).collect()
    }

    /// One genome's record.
    pub fn genome(&self, name: &str) -> Option<&GenomeTruth> {
        self.genomes.iter().find(|g| g.name == name)
    }

    /// One cluster root's record.
    pub fn root(&self, cluster: u32) -> Option<&RootTruth> {
        self.roots.iter().find(|r| r.cluster == cluster)
    }

    /// The mutations of one leaf genome.
    pub fn mutations_of(&self, name: &str) -> Option<&[Mutation]> {
        self.genome(name).map(|g| g.mutations.as_slice())
    }

    /// The true Newick of one cluster.
    pub fn tree(&self, cluster: u32) -> Option<&str> {
        self.trees.get(&cluster).map(|s| s.as_str())
    }

    /// The true SNP count between two genomes, in either order. `None` when the
    /// pair is not in the truth, which is the case across clusters.
    pub fn pairwise_snps(&self, a: &str, b: &str) -> Option<u64> {
        if a == b {
            return Some(0);
        }
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        self.pairwise
            .iter()
            .find(|(x, y, _)| x == lo && y == hi)
            .map(|&(_, _, s)| s)
    }

    /// Every leaf of one cluster, ascending by name.
    pub fn members(&self, cluster: u32) -> Vec<String> {
        self.genomes
            .iter()
            .filter(|g| g.cluster == cluster)
            .map(|g| g.name.clone())
            .collect()
    }
}

/// Parse `pairwise_snps.tsv`.
pub fn parse_pairwise(text: &str) -> Result<Vec<(String, String, u64)>, String> {
    let mut lines = text.lines();
    let header = lines.next().unwrap_or("");
    if header != "genome_a\tgenome_b\tsnps" {
        return Err(format!(
            "pairwise_snps.tsv has an unexpected header {header:?}"
        ));
    }
    let mut out = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 3 {
            return Err(format!("pairwise row needs 3 columns, got {line:?}"));
        }
        let snps = f[2]
            .parse::<u64>()
            .map_err(|e| format!("bad snp count {:?}: {e}", f[2]))?;
        out.push((f[0].to_string(), f[1].to_string(), snps));
    }
    Ok(out)
}

/// Pairwise path lengths read back out of a truth Newick string.
///
/// This is the independent route to the same numbers
/// [`crate::simulate`] writes into `pairwise_snps.tsv`: it parses the tree with
/// [`infinity_core::tree::parse_newick`] and sums branch lengths over the path
/// between every pair of leaves. Truth test 3 is exactly the statement that the
/// two agree.
pub fn pairwise_from_newick(newick: &str) -> Result<Vec<(String, String, u64)>, String> {
    let tree = infinity_core::tree::parse_newick(newick)?;
    let names = tree.leaf_names();
    let n = names.len();
    // (parent, branch length) per node index, from a flat pre-order walk.
    let mut parent: Vec<usize> = Vec::new();
    let mut blen: Vec<f64> = Vec::new();
    let mut leaf_id: BTreeMap<String, usize> = BTreeMap::new();
    flatten(&tree, usize::MAX, &mut parent, &mut blen, &mut leaf_id);

    let mut depth = vec![0.0f64; parent.len()];
    let mut edges = vec![0u32; parent.len()];
    for i in 0..parent.len() {
        if parent[i] != usize::MAX {
            depth[i] = depth[parent[i]] + blen[i];
            edges[i] = edges[parent[i]] + 1;
        }
    }
    let mut out = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let a = leaf_id[&names[i]];
            let b = leaf_id[&names[j]];
            let (mut x, mut y) = (a, b);
            while edges[x] > edges[y] {
                x = parent[x];
            }
            while edges[y] > edges[x] {
                y = parent[y];
            }
            while x != y {
                x = parent[x];
                y = parent[y];
            }
            let d = depth[a] + depth[b] - 2.0 * depth[x];
            // Branch lengths are SNP counts, so the sum is an integer; anything
            // else means the tree did not come from this simulator.
            let rounded = d.round();
            if (d - rounded).abs() > 1e-9 {
                return Err(format!(
                    "branch lengths between {} and {} sum to {d}, which is not a SNP count",
                    names[i], names[j]
                ));
            }
            let (lo, hi) = if names[i] < names[j] {
                (names[i].clone(), names[j].clone())
            } else {
                (names[j].clone(), names[i].clone())
            };
            out.push((lo, hi, rounded as u64));
        }
    }
    out.sort();
    Ok(out)
}

fn flatten(
    t: &infinity_core::tree::Tree,
    parent: usize,
    parents: &mut Vec<usize>,
    blen: &mut Vec<f64>,
    leaf_id: &mut BTreeMap<String, usize>,
) -> usize {
    let id = parents.len();
    parents.push(parent);
    blen.push(t.len());
    match t {
        infinity_core::tree::Tree::Leaf { name, .. } => {
            leaf_id.insert(name.clone(), id);
        }
        infinity_core::tree::Tree::Internal { children, .. } => {
            for c in children {
                flatten(c, id, parents, blen, leaf_id);
            }
        }
    }
    id
}

fn write_mutations(dir: &Path, name: &str, mutations: &[Mutation]) -> Result<(), String> {
    let mut text = String::from(MUTATIONS_HEADER);
    for m in mutations {
        text.push('\n');
        text.push_str(&m.to_tsv_row());
    }
    text.push('\n');
    write_text(&dir.join(mutations_file_name(name)), &text)
}

fn read_mutations(dir: &Path, name: &str) -> Result<Vec<Mutation>, String> {
    let path = dir.join(mutations_file_name(name));
    let text = read_text(&path)?;
    let mut lines = text.lines();
    let header = lines.next().unwrap_or("");
    if header != MUTATIONS_HEADER {
        return Err(format!(
            "{} has an unexpected header {header:?}",
            path.display()
        ));
    }
    lines.map(Mutation::from_tsv_row).collect()
}

fn write_text(path: &Path, text: &str) -> Result<(), String> {
    fs::write(path, text).map_err(|e| format!("cannot write {}: {e}", path.display()))
}

fn read_text(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("cannot read {}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kind_labels_round_trip() {
        for k in [
            MutationKind::Snp,
            MutationKind::Insertion,
            MutationKind::Deletion,
            MutationKind::Gain,
            MutationKind::Loss,
        ] {
            assert_eq!(MutationKind::parse(k.as_str()).unwrap(), k);
        }
        assert!(MutationKind::parse("snp").is_err());
        assert!(MutationKind::Insertion.is_insertion());
        assert!(MutationKind::Gain.is_insertion());
        assert!(MutationKind::Deletion.is_deletion());
        assert!(MutationKind::Loss.is_deletion());
        assert!(!MutationKind::Snp.is_deletion());
    }

    #[test]
    fn mutation_tsv_row_round_trips() {
        let cases = vec![
            Mutation::snp(999, b'A', b'T'),
            Mutation {
                kind: MutationKind::Insertion,
                pos: 0,
                ref_bases: Vec::new(),
                alt_bases: b"ACGTN".to_vec(),
            },
            Mutation {
                kind: MutationKind::Loss,
                pos: 42,
                ref_bases: vec![b'C'; 700],
                alt_bases: Vec::new(),
            },
        ];
        for m in cases {
            let row = m.to_tsv_row();
            assert_eq!(row.split('\t').count(), 4);
            assert_eq!(Mutation::from_tsv_row(&row).unwrap(), m);
        }
        assert_eq!(Mutation::snp(0, b'A', b'T').to_tsv_row(), "SNP\t1\tA\tT");
        assert!(Mutation::from_tsv_row("SNP\t0\tA\tT").is_err());
        assert!(Mutation::from_tsv_row("SNP\tx\tA\tT").is_err());
        assert!(Mutation::from_tsv_row("SNP\t1\tA").is_err());
    }

    #[test]
    fn length_delta_accounts_for_both_sides() {
        assert_eq!(Mutation::snp(0, b'A', b'T').length_delta(), 0);
        assert_eq!(
            Mutation {
                kind: MutationKind::Gain,
                pos: 0,
                ref_bases: Vec::new(),
                alt_bases: vec![b'A'; 1200]
            }
            .length_delta(),
            1200
        );
        assert_eq!(
            Mutation {
                kind: MutationKind::Deletion,
                pos: 0,
                ref_bases: vec![b'A'; 7],
                alt_bases: Vec::new()
            }
            .length_delta(),
            -7
        );
    }

    #[test]
    fn file_names_are_stable() {
        assert_eq!(tree_file_name(0), "tree_c00.nwk");
        assert_eq!(tree_file_name(7), "tree_c07.nwk");
        assert_eq!(tree_file_name(123), "tree_c123.nwk");
        assert_eq!(
            mutations_file_name("sim_c00_s001"),
            "mutations_sim_c00_s001.tsv"
        );
    }

    #[test]
    fn pairwise_from_newick_sums_the_path() {
        let nwk = "((a:1,b:2):3,c:4);";
        let got = pairwise_from_newick(nwk).unwrap();
        assert_eq!(
            got,
            vec![
                ("a".to_string(), "b".to_string(), 3),
                ("a".to_string(), "c".to_string(), 8),
                ("b".to_string(), "c".to_string(), 9),
            ]
        );
    }

    #[test]
    fn pairwise_from_newick_orders_by_name() {
        let nwk = "((zz:1,aa:2):0,mm:5);";
        let got = pairwise_from_newick(nwk).unwrap();
        assert_eq!(
            got,
            vec![
                ("aa".to_string(), "mm".to_string(), 7),
                ("aa".to_string(), "zz".to_string(), 3),
                ("mm".to_string(), "zz".to_string(), 6),
            ]
        );
        assert!(got.iter().all(|(a, b, _)| a < b));
        let mut sorted = got.clone();
        sorted.sort();
        assert_eq!(got, sorted, "rows must come back sorted by the pair");
    }

    #[test]
    fn pairwise_from_newick_rejects_fractional_lengths() {
        assert!(pairwise_from_newick("(a:1.5,b:0.25);").is_err());
        assert!(pairwise_from_newick("not a tree").is_err());
    }

    #[test]
    fn pairwise_from_newick_handles_one_leaf() {
        assert!(pairwise_from_newick("lonely:7;").unwrap().is_empty());
    }

    #[test]
    fn parse_pairwise_checks_its_header() {
        assert!(parse_pairwise("x\ty\tz\n").is_err());
        assert!(parse_pairwise("").is_err());
        let got = parse_pairwise("genome_a\tgenome_b\tsnps\np\tq\t12").unwrap();
        assert_eq!(got, vec![("p".to_string(), "q".to_string(), 12)]);
        assert!(parse_pairwise("genome_a\tgenome_b\tsnps\np\tq").is_err());
    }

    fn sample_truth() -> Truth {
        let config = SimConfig {
            seed: 77,
            clusters: 2,
            leaves_per_cluster: 2,
            ..SimConfig::default()
        };
        let genomes = vec![
            GenomeTruth {
                name: "sim_c00_s000".into(),
                cluster: 0,
                file: "genomes/sim_c00_s000.fasta".into(),
                length: 1000,
                contigs: 2,
                mutations: vec![
                    Mutation::snp(10, b'A', b'G'),
                    Mutation {
                        kind: MutationKind::Gain,
                        pos: 500,
                        ref_bases: Vec::new(),
                        alt_bases: vec![b'C'; 600],
                    },
                ],
            },
            GenomeTruth {
                name: "sim_c00_s001".into(),
                cluster: 0,
                file: "genomes/sim_c00_s001.fasta".into(),
                length: 1600,
                contigs: 2,
                mutations: Vec::new(),
            },
            GenomeTruth {
                name: "sim_c01_s000".into(),
                cluster: 1,
                file: "genomes/sim_c01_s000.fasta".into(),
                length: 1000,
                contigs: 2,
                mutations: vec![Mutation {
                    kind: MutationKind::Loss,
                    pos: 20,
                    ref_bases: vec![b'T'; 500],
                    alt_bases: Vec::new(),
                }],
            },
            GenomeTruth {
                name: "sim_c01_s001".into(),
                cluster: 1,
                file: "genomes/sim_c01_s001.fasta".into(),
                length: 1000,
                contigs: 2,
                mutations: Vec::new(),
            },
        ];
        let roots = vec![
            RootTruth {
                cluster: 0,
                name: "root_c00".into(),
                file: "roots/root_c00.fasta".into(),
                length: 1000,
                substitutions: 20,
                mutations: vec![Mutation::snp(3, b'A', b'T')],
            },
            RootTruth {
                cluster: 1,
                name: "root_c01".into(),
                file: "roots/root_c01.fasta".into(),
                length: 1000,
                substitutions: 20,
                mutations: Vec::new(),
            },
        ];
        Truth {
            config,
            reference: ReferenceTruth {
                name: "reference".into(),
                file: "reference.fasta".into(),
                length: 1000,
                gc: 0.5,
                source: RefSource::Random,
                source_path: None,
            },
            roots,
            genomes,
            trees: BTreeMap::from([
                (0, "(sim_c00_s000:2,sim_c00_s001:3);".to_string()),
                (1, "(sim_c01_s000:0,sim_c01_s001:1);".to_string()),
            ]),
            pairwise: vec![
                ("sim_c00_s000".into(), "sim_c00_s001".into(), 5),
                ("sim_c01_s000".into(), "sim_c01_s001".into(), 1),
            ],
        }
    }

    fn tmpdir(tag: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("infinity-sim-truth-{}-{}", std::process::id(), tag));
        let _ = fs::remove_dir_all(&p);
        p
    }

    #[test]
    fn write_then_read_is_the_identity() {
        let t = sample_truth();
        let dir = tmpdir("roundtrip");
        t.write_to(&dir).unwrap();
        let back = Truth::read_from(&dir).unwrap();
        assert_eq!(back, t);
        assert_eq!(back.mutations_of("sim_c00_s000").unwrap().len(), 2);
        assert_eq!(back.roots[0].mutations.len(), 1);
        assert_eq!(back.cluster_of("sim_c01_s000"), Some(1));
        assert_eq!(back.pairwise_snps("sim_c00_s001", "sim_c00_s000"), Some(5));
        assert_eq!(back.pairwise_snps("sim_c00_s000", "sim_c00_s000"), Some(0));
        assert_eq!(back.pairwise_snps("sim_c00_s000", "sim_c01_s000"), None);
        assert_eq!(back.members(1), vec!["sim_c01_s000", "sim_c01_s001"]);
        assert_eq!(back.cluster_ids(), vec![0, 1]);
        assert_eq!(back.genome_names().len(), 4);
        assert_eq!(back.tree(0).unwrap(), "(sim_c00_s000:2,sim_c00_s001:3);");
        assert_eq!(back.clusters().len(), 4);
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_tampered_clusters_tsv_is_caught() {
        let t = sample_truth();
        let dir = tmpdir("tampered");
        t.write_to(&dir).unwrap();
        let path = dir.join(TRUTH_DIR).join("clusters.tsv");
        let text = fs::read_to_string(&path).unwrap();
        fs::write(&path, text.replacen("sim_c00_s000", "sim_c99_s999", 1)).unwrap();
        assert!(Truth::read_from(&dir).is_err());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_missing_mutation_file_is_caught() {
        let t = sample_truth();
        let dir = tmpdir("missing");
        t.write_to(&dir).unwrap();
        fs::remove_file(dir.join(TRUTH_DIR).join("mutations_sim_c01_s000.tsv")).unwrap();
        let err = Truth::read_from(&dir).unwrap_err();
        assert!(err.contains("cannot read"), "{err}");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn an_unknown_format_version_is_refused() {
        let t = sample_truth();
        let dir = tmpdir("version");
        t.write_to(&dir).unwrap();
        let path = dir.join(TRUTH_DIR).join("manifest.json");
        let text = fs::read_to_string(&path).unwrap();
        fs::write(
            &path,
            text.replacen("\"format_version\": 1", "\"format_version\": 99", 1),
        )
        .unwrap();
        let err = Truth::read_from(&dir).unwrap_err();
        assert!(err.contains("format version"), "{err}");
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn manifest_carries_the_seed_and_config() {
        let t = sample_truth();
        let dir = tmpdir("manifest");
        t.write_to(&dir).unwrap();
        let text = fs::read_to_string(dir.join(TRUTH_DIR).join("manifest.json")).unwrap();
        assert!(text.contains("\"seed\": 77"));
        assert!(text.contains("\"mean_snps_per_branch\": 3.0"));
        assert!(text.contains("\"format_version\": 1"));
        fs::remove_dir_all(&dir).ok();
    }
}
