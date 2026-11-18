//! Result files for people and for other tools: square SNP matrices, SNP
//! clusters at fixed thresholds, a VCF and a core-SNP alignment per cluster, a
//! tree picture, and one Excel workbook holding the tables.
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
use crate::tree::Tree;
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
// Tree picture
// ---------------------------------------------------------------------------

/// Colour per cluster, the same palette order as the desktop program.
pub const CLUSTER_COLOURS: [&str; 10] = [
    "#e05d5d", "#5d8de0", "#7cc46b", "#d9a441", "#a46be0", "#3fb8b0", "#e07db8", "#9aa04a",
    "#6b7fe0", "#e0875d",
];

/// A rectangular phylogram of the full tree with branch lengths in SNPs, tip
/// labels coloured by cluster, and a scale bar in SNPs.
pub fn tree_svg(r: &RunResults) -> String {
    let t = &r.full_tree_snps;
    let leaves = t.leaf_count().max(1);
    let row = 16.0;
    let width = 1100.0;
    let label_w = 260.0;
    let top = 40.0;
    let left = 20.0;
    let depth = max_depth(t).max(1e-9);
    let sx = (width - left - label_w) / depth;
    let height = top + row * leaves as f64 + 70.0;
    let mut colour_of: BTreeMap<String, &str> = BTreeMap::new();
    for c in &r.clusters {
        let col = CLUSTER_COLOURS[c.id as usize % CLUSTER_COLOURS.len()];
        for &m in &c.members {
            colour_of.insert(name(r, m), col);
        }
    }
    let mut body = String::new();
    let mut next_row = 0usize;
    draw(
        t,
        0.0,
        &mut next_row,
        &mut body,
        left,
        top,
        sx,
        row,
        &colour_of,
    );
    // scale bar: a round number of SNPs about a fifth of the width
    let target = depth / 5.0;
    let bar = nice(target);
    let bar_px = bar * sx;
    let y = top + row * leaves as f64 + 30.0;
    let unit = if r.backbone_is_exact {
        "SNPs"
    } else {
        "SNPs (≈ between clusters)"
    };
    body.push_str(&format!(
        "<line x1=\"{left}\" y1=\"{y}\" x2=\"{x2}\" y2=\"{y}\" stroke=\"#222\" stroke-width=\"2\"/>\
         <text x=\"{left}\" y=\"{ty}\" font-size=\"12\" fill=\"#222\">{bar} {unit}</text>",
        x2 = left + bar_px,
        ty = y + 16.0,
        bar = fmt_num(bar),
    ));
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
         <svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" \
         viewBox=\"0 0 {width} {height}\" font-family=\"Helvetica, Arial, sans-serif\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n\
         <text x=\"{left}\" y=\"22\" font-size=\"14\" font-weight=\"bold\" fill=\"#222\">\
         Infinity-Win tree ({leaves} genomes, {nc} clusters, branch lengths in SNPs)</text>\n{body}\n</svg>\n",
        nc = r.clusters.len(),
    )
}

fn max_depth(t: &Tree) -> f64 {
    match t {
        Tree::Leaf { len, .. } => len.max(0.0),
        Tree::Internal { len, children } => {
            len.max(0.0) + children.iter().map(max_depth).fold(0.0, f64::max)
        }
    }
}

/// Draw `t` whose parent ends at x = `x0` (in SNPs); returns the node's y.
#[allow(clippy::too_many_arguments)]
fn draw(
    t: &Tree,
    x0: f64,
    next_row: &mut usize,
    out: &mut String,
    left: f64,
    top: f64,
    sx: f64,
    row: f64,
    colour_of: &BTreeMap<String, &str>,
) -> f64 {
    let x1 = x0 + t.len().max(0.0);
    let (px0, px1) = (left + x0 * sx, left + x1 * sx);
    match t {
        Tree::Leaf { name, .. } => {
            let y = top + row * (*next_row as f64) + row / 2.0;
            *next_row += 1;
            let col = colour_of.get(name).copied().unwrap_or("#222");
            out.push_str(&format!(
                "<line x1=\"{px0:.1}\" y1=\"{y:.1}\" x2=\"{px1:.1}\" y2=\"{y:.1}\" stroke=\"#333\" stroke-width=\"1.2\"/>\
                 <circle cx=\"{px1:.1}\" cy=\"{y:.1}\" r=\"2.5\" fill=\"{col}\"/>\
                 <text x=\"{tx:.1}\" y=\"{ty:.1}\" font-size=\"11\" fill=\"{col}\">{n}</text>\n",
                tx = px1 + 6.0,
                ty = y + 4.0,
                n = xml_escape(name),
            ));
            y
        }
        Tree::Internal { children, .. } => {
            let ys: Vec<f64> = children
                .iter()
                .map(|c| draw(c, x1, next_row, out, left, top, sx, row, colour_of))
                .collect();
            if ys.is_empty() {
                return top;
            }
            let (ymin, ymax) = (ys[0], ys[ys.len() - 1]);
            let y = (ymin + ymax) / 2.0;
            out.push_str(&format!(
                "<line x1=\"{px1:.1}\" y1=\"{ymin:.1}\" x2=\"{px1:.1}\" y2=\"{ymax:.1}\" stroke=\"#333\" stroke-width=\"1.2\"/>\
                 <line x1=\"{px0:.1}\" y1=\"{y:.1}\" x2=\"{px1:.1}\" y2=\"{y:.1}\" stroke=\"#333\" stroke-width=\"1.2\"/>\n"
            ));
            y
        }
    }
}

fn nice(x: f64) -> f64 {
    if x <= 0.0 || !x.is_finite() {
        return 1.0;
    }
    let p = 10f64.powf(x.log10().floor());
    let m = x / p;
    let n = if m < 1.5 {
        1.0
    } else if m < 3.5 {
        2.0
    } else if m < 7.5 {
        5.0
    } else {
        10.0
    };
    (n * p).max(1.0)
}

fn fmt_num(v: f64) -> String {
    if (v - v.round()).abs() < 1e-9 {
        format!("{}", v.round() as i64)
    } else {
        format!("{v}")
    }
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

pub fn write_tree_svg(dir: &Path, r: &RunResults) -> std::io::Result<()> {
    std::fs::write(dir.join("tree.svg"), tree_svg(r))
}

// ---------------------------------------------------------------------------
// Excel workbook
// ---------------------------------------------------------------------------

/// Excel's limits: rows per sheet and columns per sheet.
const XLSX_MAX_ROWS: u32 = 1_048_576;
const XLSX_MAX_COLS: usize = 16_384;

pub fn write_xlsx(dir: &Path, r: &RunResults) -> Result<(), rust_xlsxwriter::XlsxError> {
    use rust_xlsxwriter::{Format, Workbook};
    let mut wb = Workbook::new();
    let bold = Format::new().set_bold();

    // --- Summary ---
    let ws = wb.add_worksheet().set_name("Summary")?;
    let p = &r.params;
    let rows: Vec<(&str, String)> = vec![
        ("Program", format!("Infinity-Win {}", crate::VERSION)),
        ("Genomes", r.records.len().to_string()),
        ("Clusters", r.clusters.len().to_string()),
        ("Resolution", p.resolution.as_str().to_string()),
        (
            "Cluster distance (sketch)",
            format!(
                "{} (≈ {:.2}% identity)",
                p.cluster_dist,
                (1.0 - p.cluster_dist) * 100.0
            ),
        ),
        (
            "SNP thresholds",
            p.snp_thresholds
                .iter()
                .map(|t| t.to_string())
                .collect::<Vec<_>>()
                .join(", "),
        ),
        (
            "SNP rules",
            if p.snp_filter.enabled {
                "on".into()
            } else {
                "off".into()
            },
        ),
        ("Core genome only", p.snp_filter.core.to_string()),
        (
            "Minimum contig length (bp)",
            p.snp_filter.min_contig_len.to_string(),
        ),
        (
            "Contig-end margin (bp)",
            p.snp_filter.end_margin.to_string(),
        ),
        ("High-copy ratio", p.snp_filter.high_copy_ratio.to_string()),
        (
            "SNP-cluster window (bp)",
            p.snp_filter.scan_window.to_string(),
        ),
        ("SNP-cluster alpha", p.snp_filter.scan_alpha.to_string()),
        (
            "Between-cluster distances",
            if r.backbone_is_exact {
                "exact (SNPs)".into()
            } else {
                "sketch estimate".into()
            },
        ),
        ("k", p.k.to_string()),
    ];
    ws.write_string_with_format(0, 0, "Setting", &bold)?;
    ws.write_string_with_format(0, 1, "Value", &bold)?;
    for (i, (k, v)) in rows.iter().enumerate() {
        ws.write_string(i as u32 + 1, 0, *k)?;
        ws.write_string(i as u32 + 1, 1, v)?;
    }
    ws.set_column_width(0, 32)?;
    ws.set_column_width(1, 40)?;

    // --- Genomes ---
    let ws = wb.add_worksheet().set_name("Genomes")?;
    let head = [
        "id",
        "name",
        "length",
        "contigs",
        "N50",
        "GC",
        "cluster",
        "pivot",
        "focus",
        "SNPs vs pivot (all)",
        "SNPs vs pivot (pass)",
        "indels vs pivot",
        "uncovered bp",
        "novel bp",
    ];
    for (c, h) in head.iter().enumerate() {
        ws.write_string_with_format(0, c as u16, *h, &bold)?;
    }
    let labels = snp_cluster_labels(r);
    for (c, t) in p.snp_thresholds.iter().enumerate() {
        ws.write_string_with_format(
            0,
            (head.len() + c) as u16,
            format!("SNP cluster ≤{t}"),
            &bold,
        )?;
    }
    for (i, g) in r.records.iter().enumerate() {
        let row = i as u32 + 1;
        ws.write_number(row, 0, g.id as f64)?;
        ws.write_string(row, 1, name(r, g.id))?;
        ws.write_number(row, 2, g.stats.length as f64)?;
        ws.write_number(row, 3, g.stats.contigs as f64)?;
        ws.write_number(row, 4, g.stats.n50 as f64)?;
        ws.write_number(row, 5, g.stats.gc)?;
        if g.cluster != u32::MAX {
            ws.write_number(row, 6, g.cluster as f64)?;
        }
        ws.write_boolean(row, 7, g.is_pivot)?;
        ws.write_number(row, 8, g.focus as f64)?;
        ws.write_number(row, 9, g.snps_vs_pivot as f64)?;
        if let Some(&n) = r.pass_snps.get(&g.id) {
            ws.write_number(row, 10, n as f64)?;
        }
        ws.write_number(row, 11, g.indels_vs_pivot as f64)?;
        ws.write_number(row, 12, g.uncovered_bp as f64)?;
        ws.write_number(row, 13, g.novel_bp as f64)?;
        if let Some(l) = labels.get(&g.id) {
            for (c, v) in l.iter().enumerate() {
                ws.write_string(row, (head.len() + c) as u16, v)?;
            }
        }
    }
    ws.set_column_width(1, 28)?;

    // --- Clusters ---
    let ws = wb.add_worksheet().set_name("Clusters")?;
    for (c, h) in ["cluster", "pivot", "size", "focus", "core genome (bp)"]
        .iter()
        .enumerate()
    {
        ws.write_string_with_format(0, c as u16, *h, &bold)?;
    }
    for (i, id) in sorted_cluster_ids(r).iter().enumerate() {
        let Some(c) = r.clusters.iter().find(|c| c.id == *id) else {
            continue;
        };
        let row = i as u32 + 1;
        ws.write_number(row, 0, c.id as f64)?;
        ws.write_string(row, 1, name(r, c.pivot))?;
        ws.write_number(row, 2, c.size() as f64)?;
        ws.write_number(row, 3, r.cluster_focus.get(id).copied().unwrap_or(0) as f64)?;
        if let Some(&b) = r.core_bp.get(id) {
            ws.write_number(row, 4, b as f64)?;
        }
    }
    ws.set_column_width(1, 28)?;

    // --- one SNP matrix per exact cluster ---
    for id in sorted_cluster_ids(r) {
        if r.cluster_focus.get(&id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let Some(m) = r.matrices.get(&id) else {
            continue;
        };
        if m.ids.len() + 1 > XLSX_MAX_COLS || m.ids.len() < 2 {
            continue;
        }
        let ws = wb.add_worksheet().set_name(format!("Cluster {id} SNPs"))?;
        ws.write_string_with_format(0, 0, "SNPs", &bold)?;
        for (j, &g) in m.ids.iter().enumerate() {
            ws.write_string_with_format(0, j as u16 + 1, name(r, g), &bold)?;
            ws.write_string_with_format(j as u32 + 1, 0, name(r, g), &bold)?;
        }
        for i in 0..m.ids.len() {
            for j in 0..m.ids.len() {
                if i != j && m.comparable[i][j] == 0 {
                    ws.write_string(i as u32 + 1, j as u16 + 1, "NA")?;
                } else {
                    ws.write_number(i as u32 + 1, j as u16 + 1, m.snps[i][j] as f64)?;
                }
            }
        }
        ws.set_column_width(0, 24)?;
        ws.set_freeze_panes(1, 1)?;
    }

    // --- between-cluster SNPs ---
    if let Some(bb) = &r.backbone_exact {
        if bb.names.len() + 1 <= XLSX_MAX_COLS {
            let ws = wb.add_worksheet().set_name("Between clusters")?;
            ws.write_string_with_format(0, 0, "pivot SNPs", &bold)?;
            for (j, n) in bb.names.iter().enumerate() {
                ws.write_string_with_format(0, j as u16 + 1, n, &bold)?;
                ws.write_string_with_format(j as u32 + 1, 0, n, &bold)?;
                for i in 0..bb.names.len() {
                    ws.write_number(j as u32 + 1, i as u16 + 1, bb.snps[j][i] as f64)?;
                }
            }
            ws.set_column_width(0, 24)?;
        }
    }

    // --- variants ---
    let ws = wb.add_worksheet().set_name("Variants")?;
    let vh = [
        "genome",
        "cluster",
        "pivot contig",
        "pivot pos",
        "ref",
        "alt",
        "type",
        "filter",
    ];
    for (c, h) in vh.iter().enumerate() {
        ws.write_string_with_format(0, c as u16, *h, &bold)?;
    }
    let mut row = 1u32;
    'outer: for id in sorted_cluster_ids(r) {
        let Some(c) = r.clusters.iter().find(|c| c.id == id) else {
            continue;
        };
        let empty = Vec::new();
        let cnames = r.pivot_contig_names.get(&id).unwrap_or(&empty);
        let mut mids: Vec<u32> = c.members.clone();
        mids.sort_unstable();
        for mid in mids {
            let Some(cmp) = r.comparisons.get(&mid) else {
                continue;
            };
            let why = r.filter_reasons.get(&mid);
            for (i, v) in cmp.variants.iter().enumerate() {
                if row >= XLSX_MAX_ROWS {
                    break 'outer;
                }
                ws.write_string(row, 0, name(r, mid))?;
                ws.write_number(row, 1, id as f64)?;
                ws.write_string(
                    row,
                    2,
                    cnames
                        .get(v.pivot_contig as usize)
                        .cloned()
                        .unwrap_or_else(|| v.pivot_contig.to_string()),
                )?;
                ws.write_number(row, 3, (v.pivot_start + 1) as f64)?;
                ws.write_string(row, 4, truncate(&v.ref_allele))?;
                ws.write_string(row, 5, truncate(&v.alt_allele))?;
                ws.write_string(row, 6, v.kind.as_str())?;
                ws.write_string(
                    row,
                    7,
                    why.and_then(|w| w.get(i))
                        .map(String::as_str)
                        .unwrap_or(PASS),
                )?;
                row += 1;
            }
        }
    }
    ws.set_column_width(0, 24)?;
    ws.set_column_width(2, 36)?;

    wb.save(dir.join("results.xlsx"))
}

/// Excel cells hold at most 32,767 characters; long indel alleles are cut.
fn truncate(s: &str) -> String {
    if s.is_empty() {
        "-".to_string()
    } else if s.len() > 200 {
        format!("{}… ({} bp)", &s[..200], s.len())
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_bars_are_round_numbers() {
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
