//! The files in the output folder (section 9).
//!
//! Coordinates in every file are 1-based, because that is what the rest of
//! genomics uses; inside the program they are 0-based and half-open (section 5).
//! The conversion happens here and nowhere else.
//!
//! Every writer is fed from [`crate::engine::RunResults`], which is already in a
//! deterministic order, so the same run always produces byte-identical files
//! whatever the thread count (section 0 rule 3).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::engine::RunResults;
use crate::summary::RunSummary;
use crate::tree::to_newick;
use crate::variant::VariantKind;

/// Write every output file of section 9 *except* `run.json` into `dir`.
///
/// `run.json` is deliberately left out: [`crate::engine::run`] calls this
/// first, fills in the last of the timings, and only then calls
/// [`write_run_json`], so that `run.json` reports complete timings instead of
/// `null` for the stages that were still running when it was serialised.
///
/// `events.jsonl` is not written here: the event sink appends to it live, so
/// that a consumer can follow the run while it is still going.
pub fn write_tables(dir: &Path, results: &RunResults) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    write_genomes(dir, results)?;
    write_clusters(dir, results)?;
    write_tree(dir, results)?;
    for c in &results.clusters {
        write_cluster_tree(dir, results, c.id)?;
        write_cluster_distances(dir, results, c.id)?;
        write_cluster_variants(dir, results, c.id)?;
    }
    Ok(())
}

/// Write every output file of section 9 into `dir`: the tables written by
/// [`write_tables`], then `run.json`.
///
/// Callers that need the timings inside `run.json` to be complete should call
/// the two halves separately instead; see [`write_tables`].
///
/// `events.jsonl` is not written here: the event sink appends to it live, so
/// that a consumer can follow the run while it is still going.
pub fn write_all(dir: &Path, results: &RunResults, summary: &RunSummary) -> std::io::Result<()> {
    write_tables(dir, results)?;
    write_run_json(dir, results, summary)
}

fn writer(dir: &Path, name: &str) -> std::io::Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(dir.join(name))?))
}

/// `genomes.tsv`: one row per genome.
pub fn write_genomes(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    let mut w = writer(dir, "genomes.tsv")?;
    writeln!(
        w,
        "id\tname\tpath\tlength\tcontigs\tn50\tgc\tcluster\tis_pivot\tfocus\t\
         snps_vs_pivot\tindels_vs_pivot\tuncovered_bp\tnovel_bp"
    )?;
    for g in &r.records {
        let cluster = if g.cluster == u32::MAX {
            "NA".to_string()
        } else {
            g.cluster.to_string()
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.id,
            g.name,
            g.path.display(),
            g.stats.length,
            g.stats.contigs,
            g.stats.n50,
            g.stats.gc,
            cluster,
            g.is_pivot,
            g.focus,
            g.snps_vs_pivot,
            g.indels_vs_pivot,
            g.uncovered_bp,
            g.novel_bp
        )?;
    }
    w.flush()
}

/// `clusters.tsv`: one row per cluster.
pub fn write_clusters(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    let mut w = writer(dir, "clusters.tsv")?;
    writeln!(w, "cluster_id\tpivot\tpivot_name\tsize\tfocus")?;
    let mut ids: Vec<u32> = r.clusters.iter().map(|c| c.id).collect();
    ids.sort_unstable();
    for id in ids {
        let c = match r.clusters.iter().find(|c| c.id == id) {
            Some(c) => c,
            None => continue,
        };
        let focus = r.cluster_focus.get(&id).copied().unwrap_or(0);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            c.id,
            c.pivot,
            r.names[c.pivot as usize],
            c.size(),
            focus
        )?;
    }
    w.flush()
}

/// `tree.nwk`: the full tree.
pub fn write_tree(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    let mut w = writer(dir, "tree.nwk")?;
    writeln!(w, "{}", to_newick(&r.full_tree))?;
    w.flush()
}

/// `cluster_<id>.nwk`: the tree of one cluster, in SNP counts.
pub fn write_cluster_tree(dir: &Path, r: &RunResults, cluster: u32) -> std::io::Result<()> {
    let Some(t) = r.cluster_snp_trees.get(&cluster) else {
        return Ok(());
    };
    let mut w = writer(dir, &format!("cluster_{cluster}.nwk"))?;
    writeln!(w, "{}", to_newick(t))?;
    w.flush()
}

/// `cluster_<id>.dist.tsv`: pairwise SNP counts and comparable bp.
pub fn write_cluster_distances(dir: &Path, r: &RunResults, cluster: u32) -> std::io::Result<()> {
    let Some(m) = r.matrices.get(&cluster) else {
        return Ok(());
    };
    let mut w = writer(dir, &format!("cluster_{cluster}.dist.tsv"))?;
    writeln!(
        w,
        "genome_a\tname_a\tgenome_b\tname_b\tsnps\tindels\tcomparable_bp\tsubs_per_site"
    )?;
    let n = m.ids.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let na = r.names.get(m.ids[i] as usize).cloned().unwrap_or_default();
            let nb = r.names.get(m.ids[j] as usize).cloned().unwrap_or_default();
            let d = if m.comparable[i][j] == 0 {
                "NA".to_string()
            } else {
                format!("{:.12}", m.dist[i][j])
            };
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                m.ids[i], na, m.ids[j], nb, m.snps[i][j], m.indels[i][j], m.comparable[i][j], d
            )?;
        }
    }
    w.flush()
}

/// `cluster_<id>.variants.tsv`: every difference between a member and its pivot.
pub fn write_cluster_variants(dir: &Path, r: &RunResults, cluster: u32) -> std::io::Result<()> {
    let cl = match r.clusters.iter().find(|c| c.id == cluster) {
        Some(c) => c,
        None => return Ok(()),
    };
    let mut w = writer(dir, &format!("cluster_{cluster}.variants.tsv"))?;
    writeln!(
        w,
        "genome\tgenome_name\tpivot_contig\tpivot_pos\tref\talt\ttype\t\
         genome_contig\tgenome_pos\tstrand\tflags"
    )?;
    // Contig names come from the run results, so writing this file never has to
    // open a FASTA file again.
    let empty: Vec<String> = Vec::new();
    let pivot_contig_names = r.pivot_contig_names.get(&cluster).unwrap_or(&empty);

    // Members in ascending id order, so the file is byte-identical run to run.
    let mut ids: Vec<u32> = r
        .comparisons
        .iter()
        .filter(|(mid, _)| cl.members.contains(mid))
        .map(|(mid, _)| *mid)
        .collect();
    ids.sort_unstable();

    for mid in ids {
        let Some(cmp) = r.comparisons.get(&mid) else {
            continue;
        };
        let gname = r.names.get(mid as usize).cloned().unwrap_or_default();
        let genome_contigs = r.contig_names.get(&mid).unwrap_or(&empty);
        for v in &cmp.variants {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                mid,
                gname,
                name_or_index(pivot_contig_names, v.pivot_contig),
                one_based(v.pivot_start, v.kind),
                or_dash(&v.ref_allele),
                or_dash(&v.alt_allele),
                v.kind.as_str(),
                name_or_index(genome_contigs, v.genome_contig),
                // The same convention on both sides of the row. Writing
                // `genome_start + 1` unconditionally made an insertion's
                // member coordinate one greater than its pivot coordinate even
                // though, for a forward-strand insertion, they are the same
                // number: both name the base the insertion sits after.
                one_based(v.genome_start, v.kind),
                v.strand,
                v.flags_string()
            )?;
        }
    }
    w.flush()
}

/// Internal 0-based start to the 1-based position written to files.
///
/// For an insertion the reference interval is empty, so the position reported
/// is that of the base the insertion sits *after*, which is the convention that
/// makes an insertion and its matching deletion report the same coordinate.
fn one_based(start: u64, kind: VariantKind) -> u64 {
    match kind {
        VariantKind::Insertion => start,
        _ => start + 1,
    }
}

fn or_dash(s: &str) -> &str {
    if s.is_empty() {
        "-"
    } else {
        s
    }
}

/// A contig's name when it is known, otherwise its index. An unnamed contig is
/// not an error worth failing a run over, and the index is still unambiguous.
fn name_or_index(names: &[String], index: u32) -> String {
    names
        .get(index as usize)
        .cloned()
        .unwrap_or_else(|| index.to_string())
}

/// `run.json`: parameters, version, timings, counters and peak memory.
pub fn write_run_json(dir: &Path, r: &RunResults, summary: &RunSummary) -> std::io::Result<()> {
    let mut s = summary.clone();
    s.timings = r.timings.clone();
    s.counters = r.counters.clone();
    s.notes = r.notes.clone();
    s.genomes = r.records.len() as u32;
    s.clusters = r.clusters.len() as u32;
    s.peak_memory_bytes = crate::summary::peak_memory_bytes();
    let json = serde_json::to_string_pretty(&s)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    let mut w = writer(dir, "run.json")?;
    writeln!(w, "{json}")?;
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insertion_positions_are_not_shifted_by_one() {
        // An insertion has an empty reference interval, so adding one would
        // point at the base after the insertion site.
        assert_eq!(one_based(10, VariantKind::Insertion), 10);
        assert_eq!(one_based(10, VariantKind::Snp), 11);
        assert_eq!(one_based(10, VariantKind::Deletion), 11);
        assert_eq!(one_based(10, VariantKind::Complex), 11);
    }

    #[test]
    fn empty_alleles_are_written_as_a_dash() {
        assert_eq!(or_dash(""), "-");
        assert_eq!(or_dash("A"), "A");
    }
}
