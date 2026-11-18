//! Genome statistics and the scan pass (6.2).

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::fasta::{read_genome, Contig, GenomeData};
use crate::params::Parameters;
use crate::seq::{code, scan_windows, Anchor, Thresholds};

/// Per-genome numbers recorded in `genomes.tsv` (section 9).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenomeStats {
    /// Total bases, including Ns.
    pub length: u64,
    pub contigs: u32,
    pub n50: u64,
    /// (G + C) / total length, in `[0, 1]`. Zero when the genome is empty.
    pub gc: f64,
    /// Longest contig.
    pub longest: u64,
    /// Bases that are N.
    pub n_bases: u64,
}

impl Default for GenomeStats {
    fn default() -> Self {
        Self {
            length: 0,
            contigs: 0,
            n50: 0,
            gc: 0.0,
            longest: 0,
            n_bases: 0,
        }
    }
}

/// N50 from contig lengths: the length of the contig that carries the running
/// total past half the assembly.
pub fn n50(lengths: &[u64]) -> u64 {
    if lengths.is_empty() {
        return 0;
    }
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let total: u64 = sorted.iter().sum();
    let half = total / 2;
    let mut acc = 0u64;
    for &len in &sorted {
        acc += len;
        if acc > half {
            return len;
        }
    }
    sorted.last().copied().unwrap_or(0)
}

/// Compute the statistics of an already-loaded genome.
pub fn stats_of(data: &GenomeData) -> GenomeStats {
    let mut length = 0u64;
    let mut gc_count = 0u64;
    let mut n_bases = 0u64;
    let mut longest = 0u64;
    let mut lengths = Vec::with_capacity(data.contigs.len());
    for c in &data.contigs {
        let l = c.len();
        lengths.push(l);
        length += l;
        longest = longest.max(l);
        for &b in &c.seq {
            match b {
                b'G' | b'C' => gc_count += 1,
                b'N' => n_bases += 1,
                _ => {}
            }
        }
    }
    GenomeStats {
        length,
        contigs: data.contigs.len() as u32,
        n50: n50(&lengths),
        gc: if length == 0 {
            0.0
        } else {
            gc_count as f64 / length as f64
        },
        longest,
        n_bases,
    }
}

/// What the scan pass produces for one genome (6.2).
#[derive(Debug, Clone)]
pub struct Scan {
    pub stats: GenomeStats,
    /// Every level-3 anchor, ordered by contig then position.
    pub anchors: Vec<Anchor>,
    /// Sorted, de-duplicated anchor hashes at `sketch_level`.
    pub sketch: Vec<u64>,
    /// Length of each contig, in input order.
    pub contig_lens: Vec<u64>,
}

/// Scan one already-loaded genome.
pub fn scan_data(data: &GenomeData, params: &Parameters) -> Scan {
    let th = Thresholds::new(&params.densities);
    let stats = stats_of(data);
    let mut anchors: Vec<Anchor> = Vec::new();
    let mut sketch: Vec<u64> = Vec::new();
    let mut contig_lens = Vec::with_capacity(data.contigs.len());
    let sketch_idx = params.sketch_level as usize;

    for (ci, contig) in data.contigs.iter().enumerate() {
        contig_lens.push(contig.len());
        let cid = ci as u32;
        let mut local = Vec::new();
        scan_windows(&contig.seq, params.k, |pos, h, plus| {
            if let Some(level) = th.level(h) {
                let a = Anchor {
                    h,
                    contig: cid,
                    pos: pos as u64,
                    plus,
                    level,
                };
                if level <= sketch_idx as u8 {
                    sketch.push(h);
                }
                local.push(a);
            }
        });
        anchors.extend(local);
    }
    sketch.sort_unstable();
    sketch.dedup();
    Scan {
        stats,
        anchors,
        sketch,
        contig_lens,
    }
}

/// Read and scan a genome file.
pub fn scan_file(path: &Path, params: &Parameters) -> std::io::Result<(GenomeData, Scan)> {
    let data = read_genome(path)?;
    let scan = scan_data(&data, params);
    Ok((data, scan))
}

/// Everything about a genome that the whole run keeps in memory.
///
/// The anchors and the sequence are deliberately *not* here: for thousands of
/// genomes they are recomputed from the file during refinement, which is what
/// keeps the run inside the memory budget of 6.10.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GenomeRecord {
    pub id: u32,
    pub name: String,
    pub path: PathBuf,
    pub stats: GenomeStats,
    /// Sorted, de-duplicated sketch hashes.
    pub sketch: Vec<u64>,
    /// Cluster assignment; `u32::MAX` before clustering.
    pub cluster: u32,
    pub is_pivot: bool,
    /// Best focus reached for this genome: 0, 1 or 2.
    pub focus: u8,
    pub snps_vs_pivot: u64,
    pub indels_vs_pivot: u64,
    pub uncovered_bp: u64,
    pub novel_bp: u64,
}

impl GenomeRecord {
    pub fn new(id: u32, name: String, path: PathBuf, stats: GenomeStats, sketch: Vec<u64>) -> Self {
        Self {
            id,
            name,
            path,
            stats,
            sketch,
            cluster: u32::MAX,
            is_pivot: false,
            focus: 0,
            snps_vs_pivot: 0,
            indels_vs_pivot: 0,
            uncovered_bp: 0,
            novel_bp: 0,
        }
    }
}

/// The contig boundaries of a scanned genome, so a position in concatenated
/// coordinates can be turned back into `(contig, offset)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigIndex {
    /// `offsets[i]` is the start of contig `i` in concatenated coordinates, and
    /// `offsets[contigs]` is the total length.
    pub offsets: Vec<u64>,
}

impl ContigIndex {
    pub fn new(contig_lens: &[u64]) -> Self {
        let mut offsets = Vec::with_capacity(contig_lens.len() + 1);
        let mut acc = 0u64;
        offsets.push(0);
        for &l in contig_lens {
            acc += l;
            offsets.push(acc);
        }
        Self { offsets }
    }

    pub fn total(&self) -> u64 {
        *self.offsets.last().unwrap_or(&0)
    }

    pub fn contig_of(&self, pos: u64) -> Option<usize> {
        let idx = self.offsets.partition_point(|&o| o <= pos);
        if idx == 0 || idx >= self.offsets.len() {
            return None;
        }
        Some(idx - 1)
    }
}

/// Concatenate contigs into one byte string with an index. Contigs are joined
/// with a single `N` separator that is not part of any contig's coordinates;
/// the separator stops anchors and blocks from spanning a contig boundary,
/// which would be wrong because the true distance between two contigs is
/// unknown.
pub fn concat_with_index(data: &GenomeData) -> (Vec<u8>, ContigIndex) {
    let mut lens: Vec<u64> = Vec::with_capacity(data.contigs.len());
    let mut total = 0usize;
    for c in &data.contigs {
        lens.push(c.len());
        total += c.seq.len() + if lens.len() > 1 { 1 } else { 0 };
    }
    let mut buf = Vec::with_capacity(total);
    for (i, c) in data.contigs.iter().enumerate() {
        if i > 0 {
            buf.push(b'N');
        }
        buf.extend_from_slice(&c.seq);
    }
    // Offsets must point into `buf`, which has separators, so recompute.
    let mut offsets = Vec::with_capacity(data.contigs.len() + 1);
    let mut acc = 0u64;
    offsets.push(0);
    for (i, c) in data.contigs.iter().enumerate() {
        if i > 0 {
            acc += 1;
        }
        acc += c.len();
        offsets.push(acc);
    }
    (buf, ContigIndex { offsets })
}

/// The contig slice for contig `i`.
pub fn contig_slice<'a>(buf: &'a [u8], index: &ContigIndex, i: usize) -> Option<&'a [u8]> {
    let s = *index.offsets.get(i)?;
    let e = *index.offsets.get(i + 1)?;
    // The separator before contig i (i > 0) belongs to no contig.
    let s = if i > 0 { s + 1 } else { s };
    if s > e || e > buf.len() as u64 {
        return None;
    }
    Some(&buf[s as usize..e as usize])
}

/// Convenience: total number of N-free k-mers, used by tests.
pub fn count_kmers(data: &GenomeData, k: usize) -> u64 {
    let mut n = 0u64;
    for c in &data.contigs {
        n += c.seq.len().saturating_sub(k).saturating_add(1) as u64;
    }
    n
}

/// Read contigs without scanning (used by the refine pass).
pub fn read_contigs(path: &Path) -> std::io::Result<Vec<Contig>> {
    Ok(read_genome(path)?.contigs)
}

/// Fraction of bases that are usable (not N), for logging.
pub fn usable_fraction(stats: &GenomeStats) -> f64 {
    if stats.length == 0 {
        return 0.0;
    }
    1.0 - (stats.n_bases as f64 / stats.length as f64)
}

/// Sanity helper used by tests: every base code is 0..=4.
pub fn assert_normalised(data: &GenomeData) {
    for c in &data.contigs {
        for &b in &c.seq {
            debug_assert!(code(b) <= 4);
        }
    }
}

/// Re-export so downstream crates do not need to reach into `seq`.
pub use crate::seq::Thresholds as ScanThresholds;

/// Anchors of one contig only; handy for tests and for the GUI map panel.
pub fn contig_anchor_positions(scan: &Scan, contig: u32, max_level: u8) -> Vec<u64> {
    scan.anchors
        .iter()
        .filter(|a| a.contig == contig && a.level <= max_level)
        .map(|a| a.pos)
        .collect()
}

/// Anchor positions per contig for one level, as `Vec<Vec<u64>>`.
pub fn anchors_by_contig(scan: &Scan, level: u8, n_contigs: usize) -> Vec<Vec<u64>> {
    let mut out = vec![Vec::new(); n_contigs];
    for a in &scan.anchors {
        if a.level <= level {
            if let Some(v) = out.get_mut(a.contig as usize) {
                v.push(a.pos);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::Contig;

    fn genome(contigs: &[&[u8]]) -> GenomeData {
        GenomeData {
            name: "g".to_string(),
            contigs: contigs
                .iter()
                .enumerate()
                .map(|(i, s)| Contig {
                    name: format!("c{i}"),
                    seq: s.to_vec(),
                })
                .collect(),
        }
    }

    #[test]
    fn n50_is_the_standard_definition() {
        assert_eq!(n50(&[2, 2, 2, 3, 3, 3, 4]), 3);
        assert_eq!(n50(&[100]), 100);
        assert_eq!(n50(&[]), 0);
    }

    #[test]
    fn stats_count_length_contigs_gc_and_ns() {
        let g = genome(&[b"ACGT", b"NNNN", b"GGCC"]);
        let s = stats_of(&g);
        assert_eq!(s.length, 12);
        assert_eq!(s.contigs, 3);
        assert_eq!(s.n_bases, 4);
        assert_eq!(s.longest, 4);
        assert!((s.gc - 6.0 / 12.0).abs() < 1e-12);
    }

    #[test]
    fn scan_produces_sorted_deduped_sketch() {
        let mut p = Parameters::default();
        p.k = 21;
        p.densities = [64, 16, 4];
        p.sketch_level = 2;
        let seq: Vec<u8> = (0..5000u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let g = genome(&[&seq]);
        let scan = scan_data(&g, &p);
        assert!(
            scan.sketch.windows(2).all(|w| w[0] < w[1]),
            "sorted+deduped"
        );
        assert!(scan.anchors.windows(2).all(|w| w[0].pos < w[1].pos));
        // anchors are nested: level 1 subset of level 2 subset of level 3
        let l1 = scan.anchors.iter().filter(|a| a.level <= 1).count();
        let l2 = scan.anchors.iter().filter(|a| a.level <= 2).count();
        let l3 = scan.anchors.len();
        assert!(l1 <= l2 && l2 <= l3);
    }

    #[test]
    fn concat_index_round_trips_contig_slices() {
        let g = genome(&[b"AAAA", b"CCCC", b"GGGG"]);
        let (buf, idx) = concat_with_index(&g);
        assert_eq!(buf.len(), 14); // 12 bases + 2 separators
        assert_eq!(contig_slice(&buf, &idx, 0), Some(&b"AAAA"[..]));
        assert_eq!(contig_slice(&buf, &idx, 1), Some(&b"CCCC"[..]));
        assert_eq!(contig_slice(&buf, &idx, 2), Some(&b"GGGG"[..]));
        assert_eq!(contig_slice(&buf, &idx, 3), None);
        assert_eq!(idx.total(), 14);
    }

    #[test]
    fn separator_blocks_cross_contig_anchors() {
        // A k-mer spanning the separator must contain an N and be skipped, so
        // no anchor may sit within k bases of a contig junction on both sides.
        let mut p = Parameters::default();
        p.k = 5;
        p.densities = [4, 2, 1];
        let g = genome(&[b"AAAAAAAAAA", b"CCCCCCCCCC"]);
        let (buf, _idx) = concat_with_index(&g);
        let mut hits = 0;
        scan_windows(&buf, p.k, |_pos, _h, _plus| hits += 1);
        let per_contig = 10usize - p.k + 1; // 6 windows each
        assert_eq!(
            hits,
            2 * per_contig,
            "the two separator windows are skipped"
        );
    }

    #[test]
    fn contig_index_finds_the_contig_of_a_position() {
        let g = genome(&[b"AAAA", b"CCCC", b"GGGG"]);
        let (_buf, idx) = concat_with_index(&g);
        assert_eq!(idx.contig_of(0), Some(0));
        assert_eq!(idx.contig_of(3), Some(0));
        assert_eq!(idx.contig_of(5), Some(1));
        assert_eq!(idx.contig_of(13), Some(2));
        assert_eq!(idx.contig_of(14), None);
    }

    #[test]
    fn count_kmers_ignores_junctions() {
        let g = genome(&[b"AAAAAA", b"CCCCCC"]);
        assert_eq!(count_kmers(&g, 3), 4 + 4);
    }
}
