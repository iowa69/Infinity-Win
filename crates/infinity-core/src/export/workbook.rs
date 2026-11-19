//! `results.xlsx`: one Excel workbook a microbiologist can open and read.
//!
//! Sheets, in order: **Summary** (the one that opens: the run in a few lines,
//! one row per cluster, the closest pairs and the possible problems),
//! **Genomes**, **Clusters**, **SNP groups**, one coloured SNP matrix per
//! exact cluster, **Between clusters** (when the pivots were compared exactly),
//! **Variants** and **Settings**.
//!
//! Everything is read from [`RunResults`] in sorted order. The only value that
//! depends on when the file was written is the run date (day resolution, and
//! `SOURCE_DATE_EPOCH` wins when set), which is also the document's creation
//! date, so two runs on the same day give the same bytes.

use std::collections::BTreeMap;
use std::path::Path;

use rust_xlsxwriter::{
    ConditionalFormat3ColorScale, ConditionalFormatType, DocProperties, ExcelDateTime, Format,
    FormatAlign, FormatBorder, Workbook, Worksheet, XlsxError,
};

use super::svg::fmt_snps;
use super::{name, snp_cluster_labels, sorted_cluster_ids};
use crate::engine::RunResults;
use crate::filter::PASS;
use crate::profile::DistanceMatrix;
use crate::tree::Tree;

/// Excel's limits: rows per sheet and columns per sheet.
const XLSX_MAX_ROWS: u32 = 1_048_576;
const XLSX_MAX_COLS: usize = 16_384;
/// Rows of the "Closest pairs" table on the Summary sheet.
pub(super) const CLOSEST_PAIRS_MAX: usize = 500;
/// A cluster whose core genome is below this fraction of its pivot is flagged.
pub(super) const CORE_WARN_FRACTION: f64 = 0.85;
/// Excel cells hold at most 32,767 characters.
const CELL_MAX_CHARS: usize = 32_000;

// Colours: one dark blue for headings, Excel's own green-yellow-red scale.
const NAVY: &str = "#1F3864";
const HEADER_FILL: &str = "#1F4E78";
const GRID: &str = "#D9D9D9";
const GREEN: &str = "#63BE7B";
const YELLOW: &str = "#FFEB84";
const RED: &str = "#F8696B";

/// The cell formats used across the workbook.
struct Styles {
    title: Format,
    subtitle: Format,
    section: Format,
    label: Format,
    header: Format,
    text: Format,
    id: Format,
    int: Format,
    dec: Format,
    pct: Format,
    note: Format,
    problem: Format,
    ok: Format,
    m_head: Format,
    m_row: Format,
    m_cell: Format,
    m_est: Format,
    m_diag: Format,
    m_na: Format,
}

impl Styles {
    fn new() -> Self {
        let cell = Format::new()
            .set_border(FormatBorder::Thin)
            .set_border_color(GRID);
        let m_cell = Format::new()
            .set_font_size(9)
            .set_num_format("#,##0")
            .set_align(FormatAlign::Center)
            .set_border(FormatBorder::Thin)
            .set_border_color("#FFFFFF");
        Self {
            title: Format::new()
                .set_bold()
                .set_font_size(16)
                .set_font_color(NAVY),
            subtitle: Format::new().set_italic().set_font_color("#595959"),
            section: Format::new()
                .set_bold()
                .set_font_size(12)
                .set_font_color(NAVY)
                .set_border_bottom(FormatBorder::Medium)
                .set_border_bottom_color(HEADER_FILL),
            label: Format::new().set_bold(),
            header: Format::new()
                .set_bold()
                .set_font_color("#FFFFFF")
                .set_background_color(HEADER_FILL)
                .set_border(FormatBorder::Thin)
                .set_border_color("#BFBFBF")
                .set_align(FormatAlign::Center)
                .set_align(FormatAlign::VerticalCenter)
                .set_text_wrap(),
            text: cell.clone(),
            id: cell.clone().set_align(FormatAlign::Center),
            int: cell.clone().set_num_format("#,##0"),
            dec: cell.clone().set_num_format("#,##0.0"),
            pct: cell.clone().set_num_format("0.0%"),
            note: Format::new()
                .set_italic()
                .set_font_size(9)
                .set_font_color("#595959"),
            problem: Format::new()
                .set_background_color("#FFF2CC")
                .set_font_color("#7F3F00")
                .set_border(FormatBorder::Thin)
                .set_border_color("#F4B183"),
            ok: Format::new()
                .set_background_color("#E2EFDA")
                .set_font_color("#375623")
                .set_border(FormatBorder::Thin)
                .set_border_color("#A9D08E"),
            m_head: Format::new()
                .set_bold()
                .set_font_size(9)
                .set_rotation(90)
                .set_align(FormatAlign::Center)
                .set_align(FormatAlign::Bottom)
                .set_background_color("#F2F2F2")
                .set_border(FormatBorder::Thin)
                .set_border_color(GRID),
            m_row: Format::new()
                .set_bold()
                .set_font_size(9)
                .set_background_color("#F2F2F2")
                .set_border(FormatBorder::Thin)
                .set_border_color(GRID),
            m_est: m_cell.clone().set_italic(),
            m_diag: m_cell.clone().set_background_color("#A6A6A6"),
            m_na: m_cell
                .clone()
                .set_background_color("#EDEDED")
                .set_font_color("#7F7F7F"),
            m_cell,
        }
    }
}

pub fn write_xlsx(dir: &Path, r: &RunResults) -> Result<(), XlsxError> {
    build(r)?.save(dir.join("results.xlsx"))
}

/// The whole workbook, ready to save.
pub(super) fn build(r: &RunResults) -> Result<Workbook, XlsxError> {
    let st = Styles::new();
    let mut wb = Workbook::new();
    let (y, m, d) = run_date();
    let date = ExcelDateTime::from_ymd(y, m, d)?;
    wb.set_properties(
        &DocProperties::new()
            .set_title("Infinity-Win results")
            .set_comment(format!("Written by Infinity-Win {}", crate::VERSION))
            .set_creation_datetime(&date),
    );

    summary_sheet(wb.add_worksheet(), r, &st, (y, m, d))?;
    genomes_sheet(wb.add_worksheet(), r, &st)?;
    clusters_sheet(wb.add_worksheet(), r, &st)?;
    snp_groups_sheet(wb.add_worksheet(), r, &st)?;

    for id in sorted_cluster_ids(r) {
        if r.cluster_focus.get(&id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let Some(m) = r.matrices.get(&id) else {
            continue;
        };
        let n = m.ids.len();
        if n < 2 || n + 1 > XLSX_MAX_COLS {
            continue;
        }
        let names: Vec<String> = m.ids.iter().map(|&g| name(r, g)).collect();
        let order = leaf_order(r.cluster_snp_trees.get(&id), &names);
        let grid = MatrixSheet {
            title: format!(
                "Cluster {id}: pairwise SNPs after the SNP rules (pivot {})",
                r.clusters
                    .iter()
                    .find(|c| c.id == id)
                    .map(|c| name(r, c.pivot))
                    .unwrap_or_default()
            ),
            names: &names,
            order: &order,
            min_is_zero: true,
            cell: &|i, j| {
                if i == j {
                    Cell::Diag
                } else if m.comparable[i][j] == 0 {
                    Cell::Na
                } else {
                    Cell::Val(m.snps[i][j], false)
                }
            },
        };
        write_matrix_sheet(
            wb.add_worksheet().set_name(format!("Cluster {id} SNPs"))?,
            &st,
            &grid,
        )?;
    }

    if let Some(bb) = &r.backbone_exact {
        let n = bb.names.len();
        if (2..XLSX_MAX_COLS).contains(&n) {
            let order = leaf_order(Some(&bb.tree), &bb.names);
            let estimated = |i: usize, j: usize| {
                bb.estimated
                    .iter()
                    .any(|&(a, b)| (a, b) == (i.min(j), i.max(j)))
            };
            let grid = MatrixSheet {
                title: "Between clusters: SNPs between the cluster pivots".to_string(),
                names: &bb.names,
                order: &order,
                min_is_zero: false,
                cell: &|i, j| {
                    if i == j {
                        Cell::Diag
                    } else {
                        Cell::Val(bb.snps[i][j], estimated(i, j))
                    }
                },
            };
            write_matrix_sheet(wb.add_worksheet().set_name("Between clusters")?, &st, &grid)?;
        }
    }

    variants_sheet(wb.add_worksheet(), r, &st)?;
    settings_sheet(wb.add_worksheet(), r, &st)?;
    Ok(wb)
}

// ---------------------------------------------------------------------------
// Numbers behind the Summary
// ---------------------------------------------------------------------------

/// Today's date (UTC), or that of `SOURCE_DATE_EPOCH` when it is set.
fn run_date() -> (u16, u8, u8) {
    let secs = std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|s| s.trim().parse::<i64>().ok())
        .unwrap_or_else(|| {
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_secs() as i64)
                .unwrap_or(0)
        });
    civil_from_days(secs.div_euclid(86_400))
}

/// Days since 1970-01-01 to (year, month, day) in the proleptic Gregorian
/// calendar (H. Hinnant's algorithm), clamped to Excel's years 1900-9999.
pub(super) fn civil_from_days(z: i64) -> (u16, u8, u8) {
    let z = z + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    if !(1900..=9999).contains(&y) {
        return (1970, 1, 1);
    }
    (y as u16, m as u8, d as u8)
}

/// Min, median and max SNPs over the comparable pairs of a matrix.
pub(super) fn pair_stats(m: &DistanceMatrix) -> Option<(u64, f64, u64)> {
    let n = m.ids.len();
    let mut v: Vec<u64> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if m.comparable[i][j] > 0 {
                v.push(m.snps[i][j]);
            }
        }
    }
    median_of(&mut v)
}

fn median_of(v: &mut [u64]) -> Option<(u64, f64, u64)> {
    if v.is_empty() {
        return None;
    }
    v.sort_unstable();
    let k = v.len();
    let med = if k % 2 == 1 {
        v[k / 2] as f64
    } else {
        (v[k / 2 - 1] + v[k / 2]) as f64 / 2.0
    };
    Some((v[0], med, v[k - 1]))
}

/// One SNP group: threshold, cluster, group index (0-based), member rows.
pub(super) struct Group {
    pub threshold: u64,
    pub cluster: u32,
    pub index: u32,
    pub rows: Vec<usize>,
}

/// Every SNP group of two or more genomes, by threshold, cluster and group.
pub(super) fn snp_groups(r: &RunResults) -> Vec<Group> {
    let mut out = Vec::new();
    for id in sorted_cluster_ids(r) {
        if r.cluster_focus.get(&id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let Some(groups) = r.snp_groups.get(&id) else {
            continue;
        };
        for (ti, g) in groups.iter().enumerate() {
            let Some(&t) = r.params.snp_thresholds.get(ti) else {
                continue;
            };
            let mut by: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
            for (row, &k) in g.iter().enumerate() {
                by.entry(k).or_default().push(row);
            }
            for (index, rows) in by {
                if rows.len() >= 2 {
                    out.push(Group {
                        threshold: t,
                        cluster: id,
                        index,
                        rows,
                    });
                }
            }
        }
    }
    out.sort_by_key(|g| (g.threshold, g.cluster, g.index));
    out
}

/// Largest SNP count between two comparable members of `rows`.
fn max_within(m: &DistanceMatrix, rows: &[usize]) -> Option<u64> {
    let mut best = None;
    for (a, &i) in rows.iter().enumerate() {
        for &j in &rows[a + 1..] {
            if m.comparable[i][j] > 0 {
                best = Some(best.unwrap_or(0).max(m.snps[i][j]));
            }
        }
    }
    best
}

/// A pair of genomes of one cluster.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Pair {
    pub snps: u64,
    pub cluster: u32,
    pub a: String,
    pub b: String,
    pub comparable: u64,
}

/// Every within-cluster pair at most `limit` SNPs apart, closest first.
pub(super) fn closest_pairs(r: &RunResults, limit: u64) -> Vec<Pair> {
    let mut out = Vec::new();
    for id in sorted_cluster_ids(r) {
        if r.cluster_focus.get(&id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let Some(m) = r.matrices.get(&id) else {
            continue;
        };
        let n = m.ids.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if m.comparable[i][j] > 0 && m.snps[i][j] <= limit {
                    let (a, b) = (name(r, m.ids[i]), name(r, m.ids[j]));
                    let (a, b) = if a <= b { (a, b) } else { (b, a) };
                    out.push(Pair {
                        snps: m.snps[i][j],
                        cluster: id,
                        a,
                        b,
                        comparable: m.comparable[i][j],
                    });
                }
            }
        }
    }
    out.sort_by(|x, y| (x.snps, x.cluster, &x.a, &x.b).cmp(&(y.snps, y.cluster, &y.a, &y.b)));
    out
}

/// The "Possible problems" lines: divergent assemblies, then clusters whose
/// core genome is a small part of the pivot.
pub(super) fn problems(r: &RunResults) -> Vec<String> {
    let mut out = Vec::new();
    for o in &r.outliers {
        let nearest = if o.nearest_name.is_empty() {
            String::new()
        } else {
            format!(
                "; closest genome {} at ≈{:.2}% identity",
                o.nearest_name,
                (1.0 - o.nearest_dist).clamp(0.0, 1.0) * 100.0
            )
        };
        out.push(format!(
            "Divergent assembly {}: {}{nearest}. Check it, or leave it out with --exclude.",
            o.name,
            if o.reasons.is_empty() {
                "flagged as unlike the others".to_string()
            } else {
                o.reasons.join("; ")
            }
        ));
    }
    for id in sorted_cluster_ids(r) {
        let Some(c) = r.clusters.iter().find(|c| c.id == id) else {
            continue;
        };
        if c.size() < 2 || r.cluster_focus.get(&id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let (Some(&core), Some(&plen)) = (r.core_bp.get(&id), r.pivot_len.get(&id)) else {
            continue;
        };
        if plen == 0 {
            continue;
        }
        let f = core as f64 / plen as f64;
        if f < CORE_WARN_FRACTION {
            out.push(format!(
                "Cluster {id}: the core genome is {:.1}% of the pivot length ({} of {} bp). \
                 The cluster may be too diverse: lower the cluster distance.",
                f * 100.0,
                core,
                plen
            ));
        }
    }
    out
}

/// Positions of `names` in the leaf order of `tree`; names not in the tree
/// follow in their own order.
pub(super) fn leaf_order(tree: Option<&Tree>, names: &[String]) -> Vec<usize> {
    let mut index: BTreeMap<&str, usize> = BTreeMap::new();
    for (i, n) in names.iter().enumerate() {
        index.entry(n.as_str()).or_insert(i);
    }
    let mut seen = vec![false; names.len()];
    let mut out = Vec::with_capacity(names.len());
    if let Some(t) = tree {
        for leaf in t.leaf_names() {
            if let Some(&i) = index.get(leaf.as_str()) {
                if !seen[i] {
                    seen[i] = true;
                    out.push(i);
                }
            }
        }
    }
    out.extend((0..names.len()).filter(|&i| !seen[i]));
    out
}

// ---------------------------------------------------------------------------
// Sheets
// ---------------------------------------------------------------------------

fn header_row(
    ws: &mut Worksheet,
    row: u32,
    heads: &[String],
    st: &Styles,
) -> Result<(), XlsxError> {
    for (c, h) in heads.iter().enumerate() {
        ws.write_string_with_format(row, c as u16, h, &st.header)?;
    }
    Ok(())
}

fn summary_sheet(
    ws: &mut Worksheet,
    r: &RunResults,
    st: &Styles,
    (y, m, d): (u16, u8, u8),
) -> Result<(), XlsxError> {
    ws.set_name("Summary")?;
    ws.set_screen_gridlines(false);
    ws.set_landscape();
    ws.set_print_fit_to_pages(1, 0);
    let p = &r.params;
    let thresholds = &p.snp_thresholds;
    let widths = [26.0, 30.0, 11.0, 14.0, 13.0, 10.0, 10.0, 10.0];
    for (c, w) in widths.iter().enumerate() {
        ws.set_column_width(c as u16, *w)?;
    }
    ws.set_column_range_width(8, 8 + thresholds.len().max(1) as u16, 11.0)?;

    ws.write_string_with_format(0, 0, "Infinity-Win results", &st.title)?;
    ws.set_row_height(0, 24)?;
    ws.write_string_with_format(
        1,
        0,
        format!(
            "Run of {y:04}-{m:02}-{d:02} · Infinity-Win {}",
            crate::VERSION
        ),
        &st.subtitle,
    )?;
    ws.set_freeze_panes(2, 0)?;

    // --- the run in a few lines ---
    let mut row = 3u32;
    let multi = r.clusters.iter().filter(|c| c.size() >= 2).count();
    let mut lines: Vec<(&str, String)> = vec![
        (
            "Genomes",
            if r.excluded.is_empty() {
                r.records.len().to_string()
            } else {
                format!(
                    "{} ({} left out on request: {})",
                    r.records.len(),
                    r.excluded.len(),
                    r.excluded.join(", ")
                )
            },
        ),
        (
            "Clusters",
            format!(
                "{} ({multi} with 2 or more genomes, {} single)",
                r.clusters.len(),
                r.clusters.len() - multi
            ),
        ),
        (
            "Cluster distance",
            format!(
                "{} — genomes join a cluster at ≈{:.2}% identity or more",
                p.cluster_dist,
                (1.0 - p.cluster_dist) * 100.0
            ),
        ),
        (
            "Resolution",
            match p.resolution {
                crate::params::Resolution::Exact => {
                    "exact (every difference resolved to bases, SNPs counted)".into()
                }
                crate::params::Resolution::Blocks => {
                    "blocks (where genomes differ, no SNP counts)".into()
                }
                crate::params::Resolution::Sketch => {
                    "sketch (clusters and the sketch tree only)".into()
                }
            },
        ),
        (
            "SNP thresholds",
            if thresholds.is_empty() {
                "none".into()
            } else {
                thresholds
                    .iter()
                    .map(|t| format!("≤{t}"))
                    .collect::<Vec<_>>()
                    .join(", ")
                    + " SNPs"
            },
        ),
        (
            "SNP rules",
            if p.snp_filter.enabled {
                format!(
                    "on: {}contigs ≥{} bp, {} bp from contig ends, SNP-dense windows of {} bp removed",
                    if p.snp_filter.core {
                        "core genome only, "
                    } else {
                        ""
                    },
                    p.snp_filter.min_contig_len,
                    p.snp_filter.end_margin,
                    p.snp_filter.scan_window
                )
            } else {
                "off: every difference counts".into()
            },
        ),
    ];
    if r.clusters.len() > 1 {
        lines.push((
            "Between clusters",
            if r.backbone_is_exact {
                "exact SNPs between the pivots (sheet Between clusters)".into()
            } else {
                "estimated from the sketch (too far apart to align, or not computed)".into()
            },
        ));
    }
    for (k, v) in &lines {
        ws.write_string_with_format(row, 0, *k, &st.label)?;
        ws.write_string(row, 1, v)?;
        row += 1;
    }

    // --- possible problems ---
    row += 1;
    ws.write_string_with_format(row, 0, "Possible problems", &st.section)?;
    for c in 1..8u16 {
        ws.write_blank(row, c, &st.section)?;
    }
    row += 1;
    let probs = problems(r);
    if probs.is_empty() {
        ws.merge_range(
            row,
            0,
            row,
            7,
            "None found: no divergent assemblies, every cluster's core genome is at least 85% of its pivot.",
            &st.ok,
        )?;
        row += 1;
    } else {
        for line in &probs {
            ws.merge_range(row, 0, row, 7, &truncate_cell(line), &st.problem)?;
            row += 1;
        }
    }

    // --- one row per cluster ---
    row += 1;
    ws.write_string_with_format(row, 0, "Clusters", &st.section)?;
    for c in 1..8u16 {
        ws.write_blank(row, c, &st.section)?;
    }
    row += 1;
    let mut heads: Vec<String> = [
        "Cluster",
        "Pivot",
        "Genomes",
        "Core genome (bp)",
        "Core, % of pivot",
        "Min SNPs",
        "Median SNPs",
        "Max SNPs",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    heads.extend(thresholds.iter().map(|t| format!("SNP groups ≤{t}")));
    header_row(ws, row, &heads, st)?;
    ws.set_row_height(row, 30)?;
    row += 1;
    let groups = snp_groups(r);
    for id in sorted_cluster_ids(r) {
        if row >= XLSX_MAX_ROWS - 1 {
            break;
        }
        let Some(c) = r.clusters.iter().find(|c| c.id == id) else {
            continue;
        };
        let exact = r.cluster_focus.get(&id).copied().unwrap_or(0) >= 2;
        ws.write_number_with_format(row, 0, id as f64, &st.id)?;
        ws.write_string_with_format(row, 1, name(r, c.pivot), &st.text)?;
        ws.write_number_with_format(row, 2, c.size() as f64, &st.int)?;
        match (r.core_bp.get(&id), r.pivot_len.get(&id)) {
            (Some(&core), plen) if c.size() >= 2 => {
                ws.write_number_with_format(row, 3, core as f64, &st.int)?;
                match plen {
                    Some(&pl) if pl > 0 => {
                        ws.write_number_with_format(row, 4, core as f64 / pl as f64, &st.pct)?;
                    }
                    _ => {
                        ws.write_blank(row, 4, &st.text)?;
                    }
                }
            }
            _ => {
                ws.write_blank(row, 3, &st.text)?;
                ws.write_blank(row, 4, &st.text)?;
            }
        }
        let stats = if exact {
            r.matrices.get(&id).and_then(pair_stats)
        } else {
            None
        };
        match stats {
            Some((lo, med, hi)) => {
                ws.write_number_with_format(row, 5, lo as f64, &st.int)?;
                let f = if med.fract() == 0.0 { &st.int } else { &st.dec };
                ws.write_number_with_format(row, 6, med, f)?;
                ws.write_number_with_format(row, 7, hi as f64, &st.int)?;
            }
            None => {
                for col in 5..8 {
                    ws.write_blank(row, col, &st.text)?;
                }
            }
        }
        for (ti, &t) in thresholds.iter().enumerate() {
            let col = 8 + ti as u16;
            if exact && c.size() >= 2 {
                let n = groups
                    .iter()
                    .filter(|g| g.cluster == id && g.threshold == t)
                    .count();
                ws.write_number_with_format(row, col, n as f64, &st.int)?;
            } else {
                ws.write_blank(row, col, &st.text)?;
            }
        }
        row += 1;
    }
    ws.write_string_with_format(
        row,
        0,
        "Pairwise SNPs after the SNP rules. SNP groups: genomes linked by chains of pairs within the threshold \
         (single linkage); only groups of 2 or more genomes are counted, listed on the sheet SNP groups.",
        &st.note,
    )?;
    row += 2;

    // --- closest pairs ---
    let Some(&limit) = thresholds.iter().max() else {
        return Ok(());
    };
    ws.write_string_with_format(
        row,
        0,
        format!("Closest pairs (≤{limit} SNPs)"),
        &st.section,
    )?;
    for c in 1..8u16 {
        ws.write_blank(row, c, &st.section)?;
    }
    row += 1;
    let pairs = closest_pairs(r, limit);
    if pairs.is_empty() {
        ws.write_string_with_format(
            row,
            0,
            format!("No two genomes of the same cluster are within {limit} SNPs."),
            &st.note,
        )?;
        return Ok(());
    }
    let heads: Vec<String> = ["Genome A", "Genome B", "Cluster", "SNPs", "Comparable bp"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    header_row(ws, row, &heads, st)?;
    row += 1;
    for pr in pairs.iter().take(CLOSEST_PAIRS_MAX) {
        if row >= XLSX_MAX_ROWS - 1 {
            break;
        }
        ws.write_string_with_format(row, 0, &pr.a, &st.text)?;
        ws.write_string_with_format(row, 1, &pr.b, &st.text)?;
        ws.write_number_with_format(row, 2, pr.cluster as f64, &st.id)?;
        ws.write_number_with_format(row, 3, pr.snps as f64, &st.int)?;
        ws.write_number_with_format(row, 4, pr.comparable as f64, &st.int)?;
        row += 1;
    }
    if pairs.len() > CLOSEST_PAIRS_MAX {
        ws.write_string_with_format(
            row,
            0,
            format!(
                "{} pairs are within {limit} SNPs; the {CLOSEST_PAIRS_MAX} closest are shown. \
                 Every pair is on the cluster sheets and in cluster_<id>.snp_matrix.tsv.",
                pairs.len()
            ),
            &st.note,
        )?;
    }
    Ok(())
}

fn genomes_sheet(ws: &mut Worksheet, r: &RunResults, st: &Styles) -> Result<(), XlsxError> {
    ws.set_name("Genomes")?;
    let p = &r.params;
    let mut heads: Vec<String> = [
        "id",
        "name",
        "length (bp)",
        "contigs",
        "N50 (bp)",
        "GC",
        "cluster",
        "pivot",
        "focus",
        "SNPs vs pivot (all)",
        "SNPs vs pivot (pass)",
        "indels vs pivot",
        "uncovered bp",
        "novel bp",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    let fixed = heads.len();
    heads.extend(p.snp_thresholds.iter().map(|t| format!("SNP group ≤{t}")));
    heads.push("warning".to_string());
    header_row(ws, 0, &heads, st)?;
    ws.set_row_height(0, 30)?;
    let labels = snp_cluster_labels(r);
    let warn: BTreeMap<u32, String> = r
        .outliers
        .iter()
        .map(|o| (o.id, format!("divergent: {}", o.reasons.join("; "))))
        .collect();
    let mut last = 0u32;
    for (i, g) in r.records.iter().enumerate() {
        let row = i as u32 + 1;
        if row >= XLSX_MAX_ROWS {
            break;
        }
        last = row;
        ws.write_number(row, 0, g.id as f64)?;
        ws.write_string(row, 1, name(r, g.id))?;
        ws.write_number_with_format(row, 2, g.stats.length as f64, &st.int)?;
        ws.write_number(row, 3, g.stats.contigs as f64)?;
        ws.write_number_with_format(row, 4, g.stats.n50 as f64, &st.int)?;
        ws.write_number_with_format(row, 5, g.stats.gc, &st.pct)?;
        if g.cluster != u32::MAX {
            ws.write_number(row, 6, g.cluster as f64)?;
        }
        if g.is_pivot {
            ws.write_string(row, 7, "yes")?;
        }
        ws.write_number(row, 8, g.focus as f64)?;
        ws.write_number_with_format(row, 9, g.snps_vs_pivot as f64, &st.int)?;
        if let Some(&n) = r.pass_snps.get(&g.id) {
            ws.write_number_with_format(row, 10, n as f64, &st.int)?;
        }
        ws.write_number_with_format(row, 11, g.indels_vs_pivot as f64, &st.int)?;
        ws.write_number_with_format(row, 12, g.uncovered_bp as f64, &st.int)?;
        ws.write_number_with_format(row, 13, g.novel_bp as f64, &st.int)?;
        if let Some(l) = labels.get(&g.id) {
            for (c, v) in l.iter().enumerate() {
                ws.write_string(row, (fixed + c) as u16, v)?;
            }
        }
        if let Some(w) = warn.get(&g.id) {
            ws.write_string(row, (heads.len() - 1) as u16, truncate_cell(w))?;
        }
    }
    ws.set_column_width(0, 6)?;
    ws.set_column_width(1, name_width(r.records.iter().map(|g| name(r, g.id))))?;
    ws.set_column_range_width(2, (fixed - 1) as u16, 12)?;
    ws.set_column_range_width(fixed as u16, (heads.len() - 2) as u16, 11)?;
    ws.set_column_width((heads.len() - 1) as u16, 40)?;
    ws.autofilter(0, 0, last.max(1), (heads.len() - 1) as u16)?;
    ws.set_freeze_panes(1, 2)?;
    Ok(())
}

fn clusters_sheet(ws: &mut Worksheet, r: &RunResults, st: &Styles) -> Result<(), XlsxError> {
    ws.set_name("Clusters")?;
    let heads: Vec<String> = [
        "cluster",
        "pivot",
        "genomes",
        "focus",
        "pivot length (bp)",
        "core genome (bp)",
        "core, % of pivot",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    header_row(ws, 0, &heads, st)?;
    ws.set_row_height(0, 30)?;
    let mut row = 0u32;
    for id in sorted_cluster_ids(r) {
        let Some(c) = r.clusters.iter().find(|c| c.id == id) else {
            continue;
        };
        row += 1;
        if row >= XLSX_MAX_ROWS {
            break;
        }
        ws.write_number(row, 0, c.id as f64)?;
        ws.write_string(row, 1, name(r, c.pivot))?;
        ws.write_number(row, 2, c.size() as f64)?;
        ws.write_number(
            row,
            3,
            r.cluster_focus.get(&id).copied().unwrap_or(0) as f64,
        )?;
        let plen = r
            .pivot_len
            .get(&id)
            .copied()
            .or_else(|| r.records.get(c.pivot as usize).map(|g| g.stats.length));
        if let Some(pl) = plen {
            ws.write_number_with_format(row, 4, pl as f64, &st.int)?;
        }
        if let Some(&b) = r.core_bp.get(&id) {
            ws.write_number_with_format(row, 5, b as f64, &st.int)?;
            if let Some(pl) = plen.filter(|&pl| pl > 0) {
                ws.write_number_with_format(row, 6, b as f64 / pl as f64, &st.pct)?;
            }
        }
    }
    ws.set_column_width(0, 8)?;
    ws.set_column_width(1, name_width(r.clusters.iter().map(|c| name(r, c.pivot))))?;
    ws.set_column_range_width(2, 6, 13)?;
    ws.autofilter(0, 0, row.max(1), 6)?;
    ws.set_freeze_panes(1, 0)?;
    Ok(())
}

fn snp_groups_sheet(ws: &mut Worksheet, r: &RunResults, st: &Styles) -> Result<(), XlsxError> {
    ws.set_name("SNP groups")?;
    let heads: Vec<String> = [
        "SNP threshold",
        "cluster",
        "group",
        "genomes",
        "max SNPs within group",
        "members",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    header_row(ws, 0, &heads, st)?;
    ws.set_row_height(0, 30)?;
    let groups = snp_groups(r);
    let mut row = 0u32;
    for g in &groups {
        let Some(m) = r.matrices.get(&g.cluster) else {
            continue;
        };
        row += 1;
        if row >= XLSX_MAX_ROWS - 1 {
            break;
        }
        let members: Vec<String> = g.rows.iter().map(|&i| name(r, m.ids[i])).collect();
        ws.write_number(row, 0, g.threshold as f64)?;
        ws.write_number(row, 1, g.cluster as f64)?;
        ws.write_string(row, 2, format!("c{}.{}", g.cluster, g.index + 1))?;
        ws.write_number(row, 3, g.rows.len() as f64)?;
        if let Some(x) = max_within(m, &g.rows) {
            ws.write_number(row, 4, x as f64)?;
        }
        ws.write_string(row, 5, truncate_cell(&members.join(", ")))?;
    }
    if groups.is_empty() {
        ws.write_string_with_format(
            1,
            0,
            "No two genomes are linked within any SNP threshold.",
            &st.note,
        )?;
    }
    ws.set_column_range_width(0, 4, 12)?;
    ws.set_column_width(5, 90)?;
    if row > 0 {
        ws.autofilter(0, 0, row, 5)?;
    }
    ws.set_freeze_panes(1, 0)?;
    Ok(())
}

/// One cell of a square matrix.
enum Cell {
    Diag,
    Na,
    /// SNPs, and whether the value is an estimate.
    Val(u64, bool),
}

struct MatrixSheet<'a> {
    title: String,
    names: &'a [String],
    /// Matrix rows in display order.
    order: &'a [usize],
    /// Anchor green at 0 SNPs (within clusters) rather than at the lowest value.
    min_is_zero: bool,
    cell: &'a dyn Fn(usize, usize) -> Cell,
}

/// A coloured square matrix: a title, a legend line, rotated genome names as
/// column headers, rows and columns in tree order, and a green-yellow-red
/// scale over the SNP counts.
fn write_matrix_sheet(ws: &mut Worksheet, st: &Styles, g: &MatrixSheet) -> Result<(), XlsxError> {
    const TOP: u32 = 3;
    let n = g.order.len();
    let mut vals: Vec<u64> = Vec::new();
    let mut any_na = false;
    let mut any_est = false;
    for (a, &i) in g.order.iter().enumerate() {
        for &j in &g.order[a + 1..] {
            match (g.cell)(i, j) {
                Cell::Val(v, est) => {
                    vals.push(v);
                    any_est |= est;
                }
                Cell::Na => any_na = true,
                Cell::Diag => {}
            }
        }
    }
    let stats = median_of(&mut vals);
    ws.write_string_with_format(0, 0, &g.title, &st.label)?;
    let mut legend = match stats {
        Some((lo, med, hi)) => format!(
            "Colour = SNP distance: green {} · yellow median ({}) · red maximum ({}).",
            if g.min_is_zero {
                "0".to_string()
            } else {
                format!("closest ({})", fmt_snps(lo as f64))
            },
            fmt_snps(med),
            fmt_snps(hi as f64),
        ),
        None => "Colour = SNP distance: green low, yellow median, red high.".to_string(),
    };
    legend.push_str(" Dark grey diagonal: the same genome.");
    if any_na {
        legend.push_str(" NA: no shared sequence to compare.");
    }
    if any_est {
        legend.push_str(" Italic: estimated from the sketch (too divergent to align).");
    }
    legend.push_str(" Rows and columns follow the tree, so close relatives form blocks.");
    ws.write_string_with_format(1, 0, legend, &st.note)?;

    // Headers.
    ws.write_string_with_format(TOP - 1, 0, "SNPs", &st.m_row)?;
    let longest = g
        .order
        .iter()
        .map(|&i| g.names[i].chars().count())
        .max()
        .unwrap_or(4);
    ws.set_row_height(TOP - 1, (longest as f64 * 5.6 + 10.0).clamp(30.0, 160.0))?;
    for (c, &j) in g.order.iter().enumerate() {
        ws.write_string_with_format(TOP - 1, c as u16 + 1, &g.names[j], &st.m_head)?;
    }
    for (rr, &i) in g.order.iter().enumerate() {
        let row = TOP + rr as u32;
        ws.write_string_with_format(row, 0, &g.names[i], &st.m_row)?;
        for (c, &j) in g.order.iter().enumerate() {
            let col = c as u16 + 1;
            match (g.cell)(i, j) {
                Cell::Diag => {
                    ws.write_blank(row, col, &st.m_diag)?;
                }
                Cell::Na => {
                    ws.write_string_with_format(row, col, "NA", &st.m_na)?;
                }
                Cell::Val(v, est) => {
                    let f = if est { &st.m_est } else { &st.m_cell };
                    ws.write_number_with_format(row, col, v as f64, f)?;
                }
            }
        }
    }
    let hi = stats.map(|s| s.2).unwrap_or(0);
    let digits = group_digits(hi);
    ws.set_column_width(0, name_width(g.names.iter().cloned()))?;
    ws.set_column_range_width(1, n as u16, (digits as f64 + 1.5).max(3.6))?;
    let mut scale = ConditionalFormat3ColorScale::new()
        .set_minimum_color(GREEN)
        .set_midpoint_color(YELLOW)
        .set_maximum_color(RED);
    if g.min_is_zero {
        scale = scale.set_minimum(ConditionalFormatType::Number, 0);
    }
    ws.add_conditional_format(TOP, 1, TOP + n as u32 - 1, n as u16, &scale)?;
    ws.set_freeze_panes(TOP, 1)?;
    Ok(())
}

fn group_digits(v: u64) -> usize {
    let d = v.max(1).ilog10() as usize + 1;
    d + d.saturating_sub(1) / 3
}

fn variants_sheet(ws: &mut Worksheet, r: &RunResults, st: &Styles) -> Result<(), XlsxError> {
    ws.set_name("Variants")?;
    let heads: Vec<String> = [
        "genome",
        "cluster",
        "pivot contig",
        "pivot pos",
        "ref",
        "alt",
        "type",
        "filter",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    header_row(ws, 0, &heads, st)?;
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
    ws.set_column_width(7, 30)?;
    ws.autofilter(0, 0, (row - 1).max(1), 7)?;
    ws.set_freeze_panes(1, 0)?;
    Ok(())
}

fn settings_sheet(ws: &mut Worksheet, r: &RunResults, st: &Styles) -> Result<(), XlsxError> {
    ws.set_name("Settings")?;
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
            if p.snp_filter.enabled { "on" } else { "off" }.into(),
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
    header_row(ws, 0, &["Setting".to_string(), "Value".to_string()], st)?;
    let mut row = 1u32;
    for (k, v) in &rows {
        ws.write_string_with_format(row, 0, *k, &st.text)?;
        ws.write_string_with_format(row, 1, v, &st.text)?;
        row += 1;
    }
    row += 1;
    header_row(
        ws,
        row,
        &[
            "Parameter (as in run.json)".to_string(),
            "Value".to_string(),
        ],
        st,
    )?;
    row += 1;
    let mut all = Vec::new();
    if let Ok(v) = serde_json::to_value(p) {
        flatten("", &v, &mut all);
    }
    // The thread count never changes a result, and the workbook must be the
    // same bytes whatever it is (section 0 rule 3); run.json records it.
    all.retain(|(k, _)| k != "threads");
    for (k, v) in &all {
        ws.write_string_with_format(row, 0, k, &st.text)?;
        ws.write_string_with_format(row, 1, truncate_cell(v), &st.text)?;
        row += 1;
    }
    if !r.excluded.is_empty() || !r.notes.is_empty() {
        row += 1;
        header_row(ws, row, &["Notes".to_string(), String::new()], st)?;
        row += 1;
        for e in &r.excluded {
            ws.write_string(row, 0, "left out on request")?;
            ws.write_string(row, 1, e)?;
            row += 1;
        }
        for n in &r.notes {
            if row >= XLSX_MAX_ROWS {
                break;
            }
            ws.write_string(row, 0, "note")?;
            ws.write_string(row, 1, truncate_cell(n))?;
            row += 1;
        }
    }
    ws.set_column_width(0, 34)?;
    ws.set_column_width(1, 60)?;
    ws.set_freeze_panes(1, 0)?;
    Ok(())
}

/// `run.json`-style parameter names (`snp_filter.core`) and their values.
fn flatten(prefix: &str, v: &serde_json::Value, out: &mut Vec<(String, String)>) {
    use serde_json::Value;
    match v {
        Value::Object(m) => {
            for (k, x) in m {
                let key = if prefix.is_empty() {
                    k.clone()
                } else {
                    format!("{prefix}.{k}")
                };
                flatten(&key, x, out);
            }
        }
        Value::Array(a) => out.push((
            prefix.to_string(),
            a.iter().map(scalar).collect::<Vec<_>>().join(", "),
        )),
        other => out.push((prefix.to_string(), scalar(other))),
    }
}

fn scalar(v: &serde_json::Value) -> String {
    match v {
        serde_json::Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

/// Column width for a column of genome names.
fn name_width(names: impl Iterator<Item = String>) -> f64 {
    let longest = names.map(|n| n.chars().count()).max().unwrap_or(8);
    (longest as f64 + 2.0).clamp(10.0, 45.0)
}

/// Long indel alleles are cut to keep the sheet readable.
fn truncate(s: &str) -> String {
    if s.is_empty() {
        "-".to_string()
    } else if s.len() > 200 {
        format!("{}… ({} bp)", &s[..200], s.len())
    } else {
        s.to_string()
    }
}

/// Keep a text within Excel's cell limit.
fn truncate_cell(s: &str) -> String {
    if s.len() <= CELL_MAX_CHARS {
        return s.to_string();
    }
    let mut end = CELL_MAX_CHARS;
    while !s.is_char_boundary(end) {
        end -= 1;
    }
    format!("{}…", &s[..end])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::export::fixture;
    use std::io::Read;

    /// The parts of an .xlsx file (a zip archive) as text, by name.
    fn unzip(bytes: &[u8]) -> BTreeMap<String, String> {
        let u16le = |p: usize| u16::from_le_bytes([bytes[p], bytes[p + 1]]) as usize;
        let u32le = |p: usize| {
            u32::from_le_bytes([bytes[p], bytes[p + 1], bytes[p + 2], bytes[p + 3]]) as usize
        };
        let eocd = (0..bytes.len() - 21)
            .rev()
            .find(|&i| bytes[i..i + 4] == [0x50, 0x4b, 0x05, 0x06])
            .expect("end of central directory");
        let mut p = u32le(eocd + 16);
        let mut out = BTreeMap::new();
        for _ in 0..u16le(eocd + 10) {
            let (method, csize) = (u16le(p + 10), u32le(p + 20));
            let (nlen, xlen, clen) = (u16le(p + 28), u16le(p + 30), u16le(p + 32));
            let name = String::from_utf8(bytes[p + 46..p + 46 + nlen].to_vec()).unwrap();
            let off = u32le(p + 42);
            let start = off + 30 + u16le(off + 26) + u16le(off + 28);
            let data = &bytes[start..start + csize];
            let mut text = String::new();
            if method == 8 {
                flate2::read::DeflateDecoder::new(data)
                    .read_to_string(&mut text)
                    .unwrap();
            } else {
                text = String::from_utf8_lossy(data).into_owned();
            }
            out.insert(name, text);
            p += 46 + nlen + xlen + clen;
        }
        out
    }

    fn sheet_names(r: &RunResults) -> Vec<String> {
        build(r)
            .unwrap()
            .worksheets()
            .iter()
            .map(|w| w.name())
            .collect()
    }

    #[test]
    fn sheets_come_in_reading_order_with_the_summary_first() {
        let mut r = fixture();
        assert_eq!(
            sheet_names(&r),
            [
                "Summary",
                "Genomes",
                "Clusters",
                "SNP groups",
                "Cluster 0 SNPs",
                "Cluster 1 SNPs",
                "Variants",
                "Settings"
            ]
        );
        r.backbone_exact = Some(crate::backbone::ExactBackbone {
            pivots: vec![0, 4],
            names: vec!["a0".into(), "b0".into()],
            snps: vec![vec![0, 5000], vec![5000, 0]],
            comparable: vec![vec![0, 0], vec![0, 0]],
            estimated: vec![(0, 1)],
            tree: crate::tree::parse_newick("(a0:2500,b0:2500);").unwrap(),
        });
        let names = sheet_names(&r);
        let at = |s: &str| names.iter().position(|n| n == s).unwrap();
        assert!(at("Cluster 1 SNPs") < at("Between clusters"));
        assert!(at("Between clusters") < at("Variants"));
    }

    #[test]
    fn matrices_are_coloured_ordered_by_the_tree_and_mark_missing_pairs() {
        let mut r = fixture();
        // A tree order that differs from the matrix order.
        r.cluster_snp_trees.insert(
            0,
            crate::tree::parse_newick("((a3:6,a2:2):1,(a1:0,a0:0):1);").unwrap(),
        );
        let parts = unzip(&build(&r).unwrap().save_to_buffer().unwrap());
        let book = &parts["xl/workbook.xml"];
        assert!(book.find("Summary").unwrap() < book.find("Genomes").unwrap());
        let summary = &parts["xl/worksheets/sheet1.xml"];
        assert!(
            summary.contains("showGridLines=\"0\""),
            "no gridlines on the Summary"
        );
        assert!(summary.contains("<pane "), "frozen header rows");
        let m0 = &parts["xl/worksheets/sheet5.xml"];
        assert!(
            m0.contains("<conditionalFormatting sqref=\"B4:E7\">"),
            "{m0}"
        );
        assert!(m0.contains("type=\"colorScale\""));
        assert!(
            m0.contains("<cfvo type=\"num\" val=\"0\"/>"),
            "green is 0 SNPs"
        );
        assert!(m0.contains("<cfvo type=\"percentile\" val=\"50\"/>"));
        assert!(m0.contains("FF63BE7B") && m0.contains("FFFFEB84") && m0.contains("FFF8696B"));
        assert!(m0.contains("<pane "), "frozen names");
        let strings = &parts["xl/sharedStrings.xml"];
        assert!(strings.contains(">NA<"), "a2-a3 share no sequence");
        // Row 4 (the first data row) is a3, the first leaf of the tree.
        let row4 = m0.split("<row r=\"4\"").nth(1).unwrap();
        let first_cell_string = row4.split("<v>").nth(1).unwrap().split('<').next().unwrap();
        let idx: usize = first_cell_string.parse().unwrap();
        let si: Vec<&str> = strings.split("<si>").skip(1).collect();
        assert!(si[idx].contains(">a3<"), "{}", si[idx]);
    }

    #[test]
    fn closest_pairs_are_sorted_and_skip_pairs_without_shared_sequence() {
        let r = fixture();
        let pairs = closest_pairs(&r, 20);
        let got: Vec<(String, String, u64)> = pairs
            .iter()
            .map(|p| (p.a.clone(), p.b.clone(), p.snps))
            .collect();
        assert_eq!(got[0], ("a0".into(), "a1".into(), 0));
        assert!(pairs.windows(2).all(|w| w[0].snps <= w[1].snps));
        assert!(!got.iter().any(|(a, b, _)| a == "a2" && b == "a3"));
        assert!(got.contains(&("b0".into(), "b1".into(), 3)));
        assert_eq!(pairs.len(), 6, "{got:?}");
        assert_eq!(closest_pairs(&r, 2).len(), 1);
    }

    #[test]
    fn problems_name_divergent_assemblies_and_thin_core_genomes() {
        let p = problems(&fixture());
        assert_eq!(p.len(), 2, "{p:?}");
        assert!(p[0].contains("a3") && p[0].contains("GC content") && p[0].contains("97.00%"));
        assert!(p[1].starts_with("Cluster 1") && p[1].contains("70.0%"));
        assert!(p[1].contains("lower the cluster distance"));
        assert!(problems(&RunResults::default()).is_empty());
    }

    #[test]
    fn snp_groups_keep_two_or_more_genomes() {
        let r = fixture();
        let g = snp_groups(&r);
        let at5: Vec<(u32, usize)> = g
            .iter()
            .filter(|g| g.threshold == 5)
            .map(|g| (g.cluster, g.rows.len()))
            .collect();
        assert_eq!(at5, vec![(0, 3), (1, 2)]);
        let m = &r.matrices[&0];
        assert_eq!(max_within(m, &[0, 1, 2]), Some(3));
        assert_eq!(max_within(m, &[2, 3]), None, "no comparable pair");
    }

    #[test]
    fn pair_statistics_use_comparable_pairs_only() {
        let r = fixture();
        // a0-a1 0, a0-a2 3, a0-a3 7, a1-a2 3, a1-a3 7 (a2-a3 has no data).
        assert_eq!(pair_stats(&r.matrices[&0]), Some((0, 3.0, 7)));
        assert_eq!(median_of(&mut [4, 1, 3, 2]), Some((1, 2.5, 4)));
        assert_eq!(median_of(&mut []), None);
    }

    #[test]
    fn leaf_order_follows_the_tree_then_the_rest() {
        let names: Vec<String> = ["a", "b", "c", "d"].iter().map(|s| s.to_string()).collect();
        let t = crate::tree::parse_newick("((c:1,a:1):1,b:1);").unwrap();
        assert_eq!(leaf_order(Some(&t), &names), vec![2, 0, 1, 3]);
        assert_eq!(leaf_order(None, &names), vec![0, 1, 2, 3]);
    }

    #[test]
    fn dates_are_civil() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(11_017), (2000, 3, 1));
        assert_eq!(civil_from_days(20_776), (2026, 11, 19));
        assert_eq!(civil_from_days(-1), (1969, 12, 31));
    }

    #[test]
    fn long_texts_fit_in_a_cell() {
        let long = "é".repeat(CELL_MAX_CHARS);
        let t = truncate_cell(&long);
        assert!(t.len() <= CELL_MAX_CHARS + 3);
        assert!(t.ends_with('…'));
        assert_eq!(truncate_cell("short"), "short");
    }

    #[test]
    fn the_workbook_is_written_to_the_output_folder() {
        let dir = std::env::temp_dir().join(format!("infinity-xlsx-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        write_xlsx(&dir, &fixture()).unwrap();
        let bytes = std::fs::read(dir.join("results.xlsx")).unwrap();
        assert_eq!(&bytes[..2], b"PK");
        std::fs::remove_dir_all(&dir).ok();
    }
}
