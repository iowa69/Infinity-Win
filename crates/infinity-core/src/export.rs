//! Result files for people and for other tools: square SNP matrices, SNP
//! clusters at fixed thresholds, a VCF and a core-SNP alignment per cluster, a
//! tree picture ([`svg`]), and one Excel workbook holding the tables
//! ([`workbook`]).
//!
//! Every writer reads only [`RunResults`] and walks it in sorted order, so the
//! files are byte-identical whatever the thread count (section 0 rule 3).
//! Positions are 1-based, as everywhere in the output folder.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::engine::RunResults;
use crate::filter::PASS;
use crate::profile::DistanceMatrix;
use crate::variant::VariantKind;

fn writer(dir: &Path, name: &str) -> std::io::Result<BufWriter<File>> {
    Ok(BufWriter::new(File::create(dir.join(name))?))
}

fn sorted_cluster_ids(r: &RunResults) -> Vec<u32> {
    let mut ids: Vec<u32> = r.clusters.iter().map(|c| c.id).collect();
    ids.sort_unstable();
    ids
}

fn name(r: &RunResults, id: u32) -> String {
    r.names
        .get(id as usize)
        .cloned()
        .unwrap_or_else(|| id.to_string())
}

/// Write every file of this module.
pub fn write_exports(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    for id in sorted_cluster_ids(r) {
        if let Some(m) = r.matrices.get(&id) {
            if r.cluster_focus.get(&id).copied().unwrap_or(0) >= 2 {
                write_square(dir, &format!("cluster_{id}.snp_matrix.tsv"), r, m, false)?;
                if let Some(raw) = r.raw_matrices.get(&id) {
                    write_square(
                        dir,
                        &format!("cluster_{id}.raw_snp_matrix.tsv"),
                        r,
                        raw,
                        false,
                    )?;
                }
                write_vcf(dir, r, id)?;
                write_core_alignment(dir, r, id)?;
            }
        }
    }
    write_snp_clusters(dir, r)?;
    write_tree_svg(dir, r)?;
    write_xlsx(dir, r).map_err(|e| std::io::Error::other(format!("results.xlsx: {e}")))?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Square matrices
// ---------------------------------------------------------------------------

fn write_square(
    dir: &Path,
    file: &str,
    r: &RunResults,
    m: &DistanceMatrix,
    subs: bool,
) -> std::io::Result<()> {
    let mut w = writer(dir, file)?;
    let names: Vec<String> = m.ids.iter().map(|&i| name(r, i)).collect();
    writeln!(w, "genome\t{}", names.join("\t"))?;
    for (i, n) in names.iter().enumerate() {
        let row: Vec<String> = (0..names.len())
            .map(|j| {
                if i != j && m.comparable[i][j] == 0 {
                    "NA".to_string()
                } else if subs {
                    format!("{:.9}", m.dist[i][j])
                } else {
                    m.snps[i][j].to_string()
                }
            })
            .collect();
        writeln!(w, "{n}\t{}", row.join("\t"))?;
    }
    w.flush()
}

// ---------------------------------------------------------------------------
// SNP clusters
// ---------------------------------------------------------------------------

/// Genome name -> (cluster, per-threshold group label), for exact clusters.
fn snp_cluster_labels(r: &RunResults) -> BTreeMap<u32, Vec<String>> {
    let mut out = BTreeMap::new();
    for id in sorted_cluster_ids(r) {
        let (Some(m), Some(groups)) = (r.matrices.get(&id), r.snp_groups.get(&id)) else {
            continue;
        };
        for (row, &gid) in m.ids.iter().enumerate() {
            let labels = groups
                .iter()
                .map(|g| format!("c{id}.{}", g.get(row).copied().unwrap_or(0) + 1))
                .collect();
            out.insert(gid, labels);
        }
    }
    out
}

pub fn write_snp_clusters(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    let mut w = writer(dir, "snp_clusters.tsv")?;
    let th: Vec<String> = r
        .params
        .snp_thresholds
        .iter()
        .map(|t| format!("snp_le_{t}"))
        .collect();
    writeln!(w, "id\tname\tcluster\t{}", th.join("\t"))?;
    let labels = snp_cluster_labels(r);
    for g in &r.records {
        let cl = if g.cluster == u32::MAX {
            "NA".to_string()
        } else {
            g.cluster.to_string()
        };
        let cols: Vec<String> = match labels.get(&g.id) {
            Some(v) => v.clone(),
            None => vec!["NA".to_string(); th.len()],
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            g.id,
            name(r, g.id),
            cl,
            cols.join("\t")
        )?;
    }
    w.flush()
}

// ---------------------------------------------------------------------------
// VCF and core SNP alignment
// ---------------------------------------------------------------------------

/// One SNP site of a cluster: pivot contig, 0-based position, the pivot base,
/// and per member the alt base and whether it passed the SNP rules.
struct Site {
    contig: u32,
    pos: u64,
    ref_base: String,
    calls: BTreeMap<u32, (String, bool, String)>,
}

fn cluster_sites(r: &RunResults, cluster: u32) -> (Vec<u32>, Vec<Site>) {
    let Some(cl) = r.clusters.iter().find(|c| c.id == cluster) else {
        return (Vec::new(), Vec::new());
    };
    let mut members: Vec<u32> = cl.members.clone();
    members.sort_unstable();
    // pivot first, then the members in id order
    members.retain(|&m| m != cl.pivot);
    members.insert(0, cl.pivot);
    let mut sites: BTreeMap<(u32, u64), Site> = BTreeMap::new();
    for &mid in &members {
        let Some(cmp) = r.comparisons.get(&mid) else {
            continue;
        };
        let reasons = r.filter_reasons.get(&mid);
        for (i, v) in cmp.variants.iter().enumerate() {
            if v.kind != VariantKind::Snp {
                continue;
            }
            let why = reasons
                .and_then(|w| w.get(i))
                .cloned()
                .unwrap_or_else(|| PASS.to_string());
            let site = sites
                .entry((v.pivot_contig, v.pivot_start))
                .or_insert(Site {
                    contig: v.pivot_contig,
                    pos: v.pivot_start,
                    ref_base: v.ref_allele.clone(),
                    calls: BTreeMap::new(),
                });
            site.calls
                .insert(mid, (v.alt_allele.clone(), why == PASS, why));
        }
    }
    (members, sites.into_values().collect())
}

fn covered(r: &RunResults, cluster: u32, mid: u32, contig: u32, pos: u64) -> bool {
    let Some(cmp) = r.comparisons.get(&mid) else {
        // the pivot covers itself
        return r.clusters.iter().any(|c| c.id == cluster && c.pivot == mid);
    };
    let lens = r.pivot_contig_lens.get(&cluster);
    let off: u64 = lens
        .map(|l| l.iter().take(contig as usize).sum())
        .unwrap_or(0);
    !crate::interval::contains(&cmp.uncovered, off + pos)
}

pub fn write_vcf(dir: &Path, r: &RunResults, cluster: u32) -> std::io::Result<()> {
    let (members, sites) = cluster_sites(r, cluster);
    if members.is_empty() {
        return Ok(());
    }
    let mut w = writer(dir, &format!("cluster_{cluster}.vcf"))?;
    let empty = Vec::new();
    let cnames = r.pivot_contig_names.get(&cluster).unwrap_or(&empty);
    let clens = r
        .pivot_contig_lens
        .get(&cluster)
        .cloned()
        .unwrap_or_default();
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=Infinity-Win {}", crate::VERSION)?;
    writeln!(
        w,
        "##reference={} (pivot of cluster {cluster})",
        name(r, members[0])
    )?;
    for (i, n) in cnames.iter().enumerate() {
        writeln!(
            w,
            "##contig=<ID={n},length={}>",
            clens.get(i).copied().unwrap_or(0)
        )?;
    }
    writeln!(
        w,
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Genomes carrying an alternative base\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=short_contig,Description=\"Pivot contig shorter than the minimum length\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=contig_end,Description=\"Near a pivot contig end\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=high_copy,Description=\"Pivot contig of high copy number\">"
    )?;
    writeln!(
        w,
        "##FILTER=<ID=not_core,Description=\"Not covered by every genome of the cluster\">"
    )?;
    writeln!(w, "##FILTER=<ID=flag,Description=\"Flagged by the comparison engine (repeat, end extension, ...)\">")?;
    writeln!(w, "##FILTER=<ID=snp_cluster,Description=\"Inside a window with more SNPs than chance allows\">")?;
    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haploid genotype; . = not covered\">"
    )?;
    writeln!(w, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"PASS or the reasons this call does not count\">")?;
    let sample_names: Vec<String> = members.iter().map(|&m| name(r, m)).collect();
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        sample_names.join("\t")
    )?;
    for s in &sites {
        let alts: Vec<String> = s
            .calls
            .values()
            .map(|c| c.0.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let site_filter: BTreeSet<String> = s
            .calls
            .values()
            .filter(|c| !c.1)
            .flat_map(|c| {
                c.2.split(';')
                    .map(|x| x.split(':').next().unwrap_or(x).to_string())
                    .collect::<Vec<_>>()
            })
            .collect();
        let filter = if site_filter.is_empty() {
            PASS.to_string()
        } else {
            site_filter.into_iter().collect::<Vec<_>>().join(";")
        };
        let gts: Vec<String> = members
            .iter()
            .map(|&m| match s.calls.get(&m) {
                Some((alt, _, why)) => {
                    let idx = alts.iter().position(|a| a == alt).unwrap_or(0) + 1;
                    format!("{idx}:{why}")
                }
                None => {
                    if covered(r, cluster, m, s.contig, s.pos) {
                        "0:PASS".to_string()
                    } else {
                        ".:not_covered".to_string()
                    }
                }
            })
            .collect();
        writeln!(
            w,
            "{}\t{}\t.\t{}\t{}\t.\t{}\tNS={}\tGT:FT\t{}",
            cnames
                .get(s.contig as usize)
                .cloned()
                .unwrap_or_else(|| s.contig.to_string()),
            s.pos + 1,
            s.ref_base,
            alts.join(","),
            filter,
            s.calls.len(),
            gts.join("\t")
        )?;
    }
    w.flush()
}

/// The core SNP alignment: one sequence per genome, one column per site where
/// at least one genome has a passing SNP. A genome whose call at a site did not
/// pass the rules gets `N` there. Ready for IQ-TREE (with `+ASC`), RAxML or
/// Gubbins-style tools.
pub fn write_core_alignment(dir: &Path, r: &RunResults, cluster: u32) -> std::io::Result<()> {
    let (members, sites) = cluster_sites(r, cluster);
    let sites: Vec<&Site> = sites
        .iter()
        .filter(|s| s.calls.values().any(|c| c.1))
        .collect();
    if members.is_empty() {
        return Ok(());
    }
    let mut w = writer(dir, &format!("cluster_{cluster}.core_snps.fasta"))?;
    for &m in &members {
        let seq: String = sites
            .iter()
            .map(|s| match s.calls.get(&m) {
                Some((alt, true, _)) => alt.chars().next().unwrap_or('N'),
                Some((_, false, _)) => 'N',
                None => s.ref_base.chars().next().unwrap_or('N'),
            })
            .collect();
        writeln!(w, ">{}", name(r, m))?;
        for chunk in seq.as_bytes().chunks(80) {
            w.write_all(chunk)?;
            w.write_all(b"\n")?;
        }
        if seq.is_empty() {
            w.write_all(b"\n")?;
        }
    }
    w.flush()
}

// ---------------------------------------------------------------------------
// Tree picture and Excel workbook
// ---------------------------------------------------------------------------

mod svg;
mod workbook;

pub use svg::{tree_svg, CLUSTER_COLOURS};
pub use workbook::write_xlsx;

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

pub fn write_tree_svg(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    std::fs::write(dir.join("tree.svg"), tree_svg(r))
}

/// A small finished run for the tests of this module and its children: two
/// clusters (a0..a3 and b0, b1), a long estimated backbone, one pair without
/// shared sequence (a2-a3), a divergent assembly (a3) and a cluster whose core
/// genome is 70% of its pivot (cluster 1).
#[cfg(test)]
fn fixture() -> RunResults {
    use crate::cluster::Cluster;
    use crate::genome::{GenomeRecord, GenomeStats};
    use crate::tree::parse_newick;

    let names = ["a0", "a1", "a2", "a3", "b0", "b1"];
    let mut r = RunResults {
        names: names.iter().map(|s| s.to_string()).collect(),
        ..RunResults::default()
    };
    for (i, n) in names.iter().enumerate() {
        let stats = GenomeStats {
            length: 1_000_000,
            ..GenomeStats::default()
        };
        let mut g = GenomeRecord::new(
            i as u32,
            n.to_string(),
            std::path::PathBuf::from(format!("{n}.fasta")),
            stats,
            Vec::new(),
        );
        g.cluster = if i < 4 { 0 } else { 1 };
        g.is_pivot = i == 0 || i == 4;
        g.focus = 2;
        r.records.push(g);
    }
    r.clusters = vec![
        Cluster {
            id: 0,
            pivot: 0,
            members: vec![0, 1, 2, 3],
        },
        Cluster {
            id: 1,
            pivot: 4,
            members: vec![4, 5],
        },
    ];
    let snps0 = vec![
        vec![0, 0, 3, 7],
        vec![0, 0, 3, 7],
        vec![3, 3, 0, 8],
        vec![7, 7, 8, 0],
    ];
    let mut comparable0 = vec![vec![900_000u64; 4]; 4];
    comparable0[2][3] = 0;
    comparable0[3][2] = 0;
    let snps1 = vec![vec![0, 3], vec![3, 0]];
    for (id, ids, snps, comparable, core) in [
        (0u32, vec![0u32, 1, 2, 3], snps0, comparable0, 950_000u64),
        (1, vec![4, 5], snps1, vec![vec![600_000u64; 2]; 2], 700_000),
    ] {
        r.snp_groups.insert(
            id,
            crate::filter::snp_groups(&snps, &r.params.snp_thresholds),
        );
        r.matrices.insert(
            id,
            DistanceMatrix {
                ids,
                snps,
                comparable,
                ..DistanceMatrix::default()
            },
        );
        r.cluster_focus.insert(id, 2);
        r.core_bp.insert(id, core);
        r.pivot_len.insert(id, 1_000_000);
    }
    r.cluster_snp_trees.insert(
        0,
        parse_newick("((a0:0,a1:0):1,a2:2,a3:6);").expect("newick"),
    );
    r.cluster_snp_trees
        .insert(1, parse_newick("(b0:1,b1:2);").expect("newick"));
    r.full_tree_snps =
        parse_newick("(((a0:0,a1:0):1,a2:2,a3:6):5000,(b0:1,b1:2):5000);").expect("newick");
    r.backbone_is_exact = false;
    r.outliers.push(crate::events::OutlierInfo {
        id: 3,
        name: "a3".into(),
        length: 1_000_000,
        gc: 0.5,
        nearest_dist: 0.03,
        nearest_name: "a2".into(),
        reasons: vec!["GC content far from the others".into()],
    });
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_bars_are_round_numbers() {
        use svg::nice;
        assert_eq!(nice(0.3), 1.0);
        assert_eq!(nice(3.0), 2.0);
        assert_eq!(nice(4.0), 5.0);
        assert_eq!(nice(80.0), 100.0);
        assert_eq!(nice(1234.0), 1000.0);
    }

    #[test]
    fn labels_are_escaped_for_svg() {
        assert_eq!(xml_escape("a<b&\"c\""), "a&lt;b&amp;&quot;c&quot;");
    }
}
