//! The tree picture, `tree.svg`: a rectangular phylogram meant for a report or
//! a paper.
//!
//! Drawing rules: one line width everywhere, square elbows, every vertical
//! connector spans exactly its first to its last child, and no line runs past
//! a node or a tip (butt caps; the vertical connector is extended by half a
//! line width so that the corners close). Tip labels sit right of the tips
//! with the genome's own branch length in SNPs (`+3`; `+0` means identical to
//! its sister). Clusters are marked by a coloured strip with a label, branches
//! whose length is an estimate are drawn in light grey, and the SNP groups of
//! the smallest threshold are shaded as soft bands behind their clade.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use super::name;
use crate::engine::RunResults;
use crate::tree::Tree;

/// Colour per cluster, the same palette order as the desktop program.
pub const CLUSTER_COLOURS: [&str; 10] = [
    "#e05d5d", "#5d8de0", "#7cc46b", "#d9a441", "#a46be0", "#3fb8b0", "#e07db8", "#9aa04a",
    "#6b7fe0", "#e0875d",
];

/// Height of one tip row, px.
const ROW: f64 = 18.0;
/// Tip label font size, px.
const FONT: f64 = 11.0;
/// Average advance of a glyph at [`FONT`] in Helvetica/Arial, px. Used only to
/// place the cluster strips right of the longest label.
const CHAR_W: f64 = 6.4;
/// Width of the tree itself (root to deepest tip), px.
const TREE_W: f64 = 560.0;
const MARGIN: f64 = 24.0;
/// Short stub drawn left of the root.
const ROOT_STUB: f64 = 10.0;
const STROKE: f64 = 1.3;
/// Shorten the branches between clusters when the deepest path inside a
/// cluster is less than this fraction of the deepest path overall...
const BREAK_RATIO: f64 = 0.25;
/// ...and give them this share of the tree width.
const BACKBONE_SHARE: f64 = 0.3;
/// Space above the first tip row: title, subtitle.
const HEADER: f64 = 64.0;
/// Space below the last tip row: scale bar and legend.
const FOOTER: f64 = 64.0;

const INK: &str = "#1f2328";
const MUTED: &str = "#6a737d";
const ESTIMATE: &str = "#a8adb3";
const BANDS: [&str; 2] = ["#fbefc5", "#dcebf7"];

/// One node of the laid-out tree.
struct Node {
    /// x of the parent's vertical connector, px.
    x_parent: f64,
    /// x of this node (the tip for a leaf), px.
    x: f64,
    /// y of this node, px.
    y: f64,
    /// Branch length into this node, SNPs.
    len: f64,
    /// First and last tip row below this node.
    first: usize,
    last: usize,
    parent: Option<usize>,
    children: Vec<usize>,
    leaf: Option<String>,
    /// The branch into this node lies inside one cluster.
    within: bool,
    /// The branch into this node is an estimate rather than an exact count.
    estimate: bool,
}

/// The tree in drawing order. Nodes are stored children first (post-order),
/// so walking them backwards visits every parent before its children.
struct Layout {
    nodes: Vec<Node>,
    root: usize,
    /// Tip names in drawing order (top to bottom).
    tips: Vec<String>,
}

fn row_y(row: usize) -> f64 {
    HEADER + ROW * row as f64 + ROW / 2.0
}

/// Rows and y positions; x comes later from [`place`].
fn layout(t: &Tree) -> Layout {
    fn walk(t: &Tree, nodes: &mut Vec<Node>, tips: &mut Vec<String>) -> usize {
        let len = t.len().max(0.0);
        let (first, last, y, children, leaf) = match t {
            Tree::Leaf { name, .. } => {
                let row = tips.len();
                tips.push(name.clone());
                (row, row, row_y(row), Vec::new(), Some(name.clone()))
            }
            Tree::Internal { children, .. } => {
                let kids: Vec<usize> = children.iter().map(|c| walk(c, nodes, tips)).collect();
                match (kids.first(), kids.last()) {
                    (Some(&a), Some(&b)) => (
                        nodes[a].first,
                        nodes[b].last,
                        (nodes[a].y + nodes[b].y) / 2.0,
                        kids,
                        None,
                    ),
                    _ => (tips.len(), tips.len(), row_y(tips.len()), kids, None),
                }
            }
        };
        nodes.push(Node {
            x_parent: 0.0,
            x: 0.0,
            y,
            len,
            first,
            last,
            parent: None,
            children,
            leaf,
            within: false,
            estimate: false,
        });
        let me = nodes.len() - 1;
        for k in nodes[me].children.clone() {
            nodes[k].parent = Some(me);
        }
        me
    }
    let mut nodes = Vec::new();
    let mut tips = Vec::new();
    let root = walk(t, &mut nodes, &mut tips);
    Layout { nodes, root, tips }
}

/// Mark each branch as inside a cluster or between clusters, and as an
/// estimate: between clusters unless the backbone is exact, inside clusters
/// that never reached exact SNPs.
fn classify(l: &mut Layout, r: &RunResults, cluster_of: &BTreeMap<String, u32>) {
    let mut in_tree: BTreeMap<u32, usize> = BTreeMap::new();
    for tip in &l.tips {
        if let Some(&c) = cluster_of.get(tip) {
            *in_tree.entry(c).or_insert(0) += 1;
        }
    }
    for node in &mut l.nodes {
        let (first, last) = (node.first, node.last);
        if first > last || last >= l.tips.len() {
            node.estimate = !r.backbone_is_exact;
            continue;
        }
        let mut uniform = cluster_of.get(&l.tips[first]).copied();
        for tip in &l.tips[first..=last] {
            if cluster_of.get(tip).copied() != uniform {
                uniform = None;
                break;
            }
        }
        let n = last - first + 1;
        // The branch above a whole cluster is the backbone's.
        node.within = matches!(uniform, Some(c) if n < in_tree.get(&c).copied().unwrap_or(0));
        node.estimate = match uniform {
            Some(c) if node.within => r.cluster_focus.get(&c).copied().unwrap_or(0) < 2,
            _ => !r.backbone_is_exact,
        };
    }
}

/// Deepest path from the root, counting only the branches `keep` selects.
fn depth_of(l: &Layout, keep: impl Fn(&Node) -> bool) -> f64 {
    let mut d = vec![0.0f64; l.nodes.len()];
    let mut best: f64 = 0.0;
    for i in (0..l.nodes.len()).rev() {
        let n = &l.nodes[i];
        let up = n.parent.map(|p| d[p]).unwrap_or(0.0);
        d[i] = up + if keep(n) { n.len } else { 0.0 };
        best = best.max(d[i]);
    }
    best
}

/// Horizontal positions: branches inside clusters at `sx_in` px per SNP,
/// between clusters at `sx_bb`.
fn place(l: &mut Layout, x0: f64, sx_in: f64, sx_bb: f64) {
    for i in (0..l.nodes.len()).rev() {
        let xp = l.nodes[i].parent.map(|p| l.nodes[p].x).unwrap_or(x0);
        let n = &mut l.nodes[i];
        n.x_parent = xp;
        n.x = xp + n.len * if n.within { sx_in } else { sx_bb };
    }
}

/// A count of SNPs for a label: whole numbers (with thousands separators)
/// when the value is one or is large, otherwise one decimal.
pub(super) fn fmt_snps(v: f64) -> String {
    let v = v.max(0.0);
    let r = v.round();
    if (v - r).abs() < 0.05 || v >= 100.0 {
        group_thousands(r as u64)
    } else {
        format!("{v:.1}")
    }
}

fn group_thousands(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

fn text_w(s: &str) -> f64 {
    s.chars().count() as f64 * CHAR_W
}

/// Maximal runs of consecutive rows with the same key, as `(key, first, last)`.
fn runs<K: Copy + PartialEq>(keys: &[Option<K>]) -> Vec<(K, usize, usize)> {
    let mut out: Vec<(K, usize, usize)> = Vec::new();
    for (i, k) in keys.iter().enumerate() {
        let Some(k) = *k else { continue };
        match out.last_mut() {
            Some((pk, _, last)) if *pk == k && *last + 1 == i => *last = i,
            _ => out.push((k, i, i)),
        }
    }
    out
}

/// x of the most recent common ancestor of tip rows `[a, b]`.
fn mrca_x(l: &Layout, a: usize, b: usize) -> f64 {
    l.nodes
        .iter()
        .filter(|n| n.first <= a && n.last >= b)
        .map(|n| n.x)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// SNP groups of at least two genomes at the smallest threshold, as sets of
/// tip names, in cluster then group order.
fn smallest_threshold_groups(r: &RunResults) -> Option<(u64, Vec<Vec<String>>)> {
    let (ti, &t) = r
        .params
        .snp_thresholds
        .iter()
        .enumerate()
        .min_by_key(|(_, &t)| t)?;
    let mut out = Vec::new();
    for id in super::sorted_cluster_ids(r) {
        if r.cluster_focus.get(&id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let (Some(m), Some(groups)) = (r.matrices.get(&id), r.snp_groups.get(&id)) else {
            continue;
        };
        let Some(g) = groups.get(ti) else {
            continue;
        };
        let mut by_group: BTreeMap<u32, Vec<String>> = BTreeMap::new();
        for (row, &gid) in m.ids.iter().enumerate() {
            if let Some(&k) = g.get(row) {
                by_group.entry(k).or_default().push(name(r, gid));
            }
        }
        out.extend(by_group.into_values().filter(|v| v.len() >= 2));
    }
    Some((t, out))
}

/// A rectangular phylogram of the full tree with branch lengths in SNPs.
pub fn tree_svg(r: &RunResults) -> String {
    let t = &r.full_tree_snps;
    let mut cluster_of: BTreeMap<String, u32> = BTreeMap::new();
    for c in &r.clusters {
        for &m in &c.members {
            cluster_of.insert(name(r, m), c.id);
        }
    }
    let tree_left = MARGIN + ROOT_STUB;
    let mut l = layout(t);
    classify(&mut l, r, &cluster_of);
    // When the branches between clusters dwarf those inside them, shorten the
    // former (drawn with a break and their length written on them) so that
    // the structure inside the clusters stays readable.
    let total = depth_of(&l, |_| true).max(1.0);
    let inside = depth_of(&l, |n| n.within);
    let between = depth_of(&l, |n| !n.within);
    let broken = inside >= 1.0 && between > 0.0 && inside / total < BREAK_RATIO;
    let (sx_in, sx_bb, scale_depth) = if broken {
        (
            TREE_W * (1.0 - BACKBONE_SHARE) / inside,
            TREE_W * BACKBONE_SHARE / between,
            inside,
        )
    } else {
        (TREE_W / total, TREE_W / total, total)
    };
    place(&mut l, tree_left, sx_in, sx_bb);
    let n_tips = l.tips.len();

    // Tip labels: the name, then the terminal branch in SNPs.
    let mut labels_end = tree_left + TREE_W;
    let mut tip_text = String::new();
    for node in &l.nodes {
        let Some(tip) = &node.leaf else { continue };
        let snps = format!(
            "{}{}",
            if node.estimate { "≈+" } else { "+" },
            fmt_snps(node.len)
        );
        let x = node.x + 5.0;
        let _ = writeln!(
            tip_text,
            "<text x=\"{x:.1}\" y=\"{y:.1}\">{n}<tspan dx=\"6\" fill=\"{MUTED}\" font-size=\"{fs}\">{snps}</tspan></text>",
            y = node.y + FONT * 0.36,
            n = super::xml_escape(tip),
            fs = FONT - 1.0,
        );
        labels_end = labels_end.max(x + text_w(tip) + 6.0 + text_w(&snps));
    }

    // Cluster strips right of the longest label.
    let strip_x = labels_end + 14.0;
    let mut strips = String::new();
    let mut cluster_label_w: f64 = 0.0;
    let size_of: BTreeMap<u32, usize> = r.clusters.iter().map(|c| (c.id, c.size())).collect();
    let cluster_keys: Vec<Option<u32>> =
        l.tips.iter().map(|t| cluster_of.get(t).copied()).collect();
    for (id, a, b) in runs(&cluster_keys) {
        let col = CLUSTER_COLOURS[id as usize % CLUSTER_COLOURS.len()];
        let size = size_of.get(&id).copied().unwrap_or(0);
        let y0 = HEADER + ROW * a as f64 + 2.0;
        let y1 = HEADER + ROW * (b + 1) as f64 - 2.0;
        let label = if size <= 1 {
            format!("Cluster {id}")
        } else {
            format!("Cluster {id} (n={size})")
        };
        cluster_label_w = cluster_label_w.max(text_w(&label));
        let _ = writeln!(
            strips,
            "<rect x=\"{strip_x:.1}\" y=\"{y0:.1}\" width=\"5\" height=\"{h:.1}\" fill=\"{col}\"/>\
             <text x=\"{tx:.1}\" y=\"{ty:.1}\" fill=\"{INK}\">{label}</text>",
            h = (y1 - y0).max(1.0),
            tx = strip_x + 11.0,
            ty = (y0 + y1) / 2.0 + FONT * 0.36,
        );
    }

    // SNP-group bands behind their clade.
    let mut bands = String::new();
    let mut band_threshold = None;
    if let Some((th, groups)) = smallest_threshold_groups(r) {
        let mut group_of: BTreeMap<&str, usize> = BTreeMap::new();
        for (k, g) in groups.iter().enumerate() {
            for n in g {
                group_of.insert(n.as_str(), k);
            }
        }
        let keys: Vec<Option<usize>> = l
            .tips
            .iter()
            .map(|t| group_of.get(t.as_str()).copied())
            .collect();
        // One colour per group (alternating in tree order), so that a group
        // split by an outsider keeps its colour.
        let mut colour_of_group: BTreeMap<usize, usize> = BTreeMap::new();
        for (g, a, b) in runs(&keys) {
            let next = colour_of_group.len();
            let k = *colour_of_group.entry(g).or_insert(next);
            let x0 = (mrca_x(&l, a, b) - 6.0).max(MARGIN);
            let y0 = HEADER + ROW * a as f64 + 1.0;
            let y1 = HEADER + ROW * (b + 1) as f64 - 1.0;
            let _ = writeln!(
                bands,
                "<rect x=\"{x0:.1}\" y=\"{y0:.1}\" width=\"{w:.1}\" height=\"{h:.1}\" rx=\"3\" fill=\"{fill}\"/>",
                w = labels_end + 4.0 - x0,
                h = y1 - y0,
                fill = BANDS[k % BANDS.len()],
            );
            band_threshold = Some(th);
        }
    }

    // Branches: horizontals from each parent connector to the node, then the
    // vertical connector of every internal node.
    let half = STROKE / 2.0;
    let (mut exact_d, mut est_d) = (String::new(), String::new());
    let root = &l.nodes[l.root];
    if n_tips > 0 {
        let stub_est =
            !root.children.is_empty() && root.children.iter().all(|&c| l.nodes[c].estimate);
        let d = if stub_est { &mut est_d } else { &mut exact_d };
        let _ = write!(d, "M{:.1},{:.1}H{:.1}", MARGIN, root.y, root.x_parent);
    }
    for (i, node) in l.nodes.iter().enumerate() {
        let d = if node.estimate {
            &mut est_d
        } else {
            &mut exact_d
        };
        // The incoming branch starts half a line width left of the parent's
        // connector so that it joins it without a gap (it is covered there).
        if node.x > node.x_parent {
            let start = if i == l.root {
                node.x_parent
            } else {
                node.x_parent - half
            };
            let _ = write!(d, "M{start:.1},{:.1}H{:.1}", node.y, node.x);
        }
        if node.children.len() >= 2 {
            let (Some(&a), Some(&b)) = (node.children.first(), node.children.last()) else {
                continue;
            };
            let all_est = node.children.iter().all(|&c| l.nodes[c].estimate);
            let d = if all_est { &mut est_d } else { &mut exact_d };
            let _ = write!(
                d,
                "M{x:.1},{y0:.1}V{y1:.1}",
                x = node.x,
                y0 = l.nodes[a].y - half,
                y1 = l.nodes[b].y + half,
            );
        }
    }

    // Shortened branches: a break mark and the true length written above.
    let mut breaks = String::new();
    if broken {
        for node in &l.nodes {
            let px = node.x - node.x_parent;
            if node.within || node.len < 1.0 || px < 24.0 {
                continue;
            }
            let xm = (node.x_parent + node.x) / 2.0;
            let y = node.y;
            let col = if node.estimate { ESTIMATE } else { INK };
            let _ = write!(
                breaks,
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"4\" height=\"{:.1}\" fill=\"#ffffff\"/>\
                 <path d=\"M{:.1},{:.1}L{:.1},{:.1}M{:.1},{:.1}L{:.1},{:.1}\" stroke=\"{col}\" stroke-width=\"{STROKE}\" fill=\"none\"/>",
                xm - 2.0,
                y - STROKE,
                2.0 * STROKE,
                xm - 5.0,
                y + 4.5,
                xm - 1.0,
                y - 4.5,
                xm + 1.0,
                y + 4.5,
                xm + 5.0,
                y - 4.5,
            );
            let label = format!(
                "{}{} SNPs",
                if node.estimate { "≈" } else { "" },
                fmt_snps(node.len)
            );
            let _ = writeln!(
                breaks,
                "<text x=\"{xm:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"{}\" fill=\"{MUTED}\">{label}</text>",
                y - 6.0,
                FONT - 1.5,
            );
        }
    }

    // Scale bar: a round number of SNPs about a fifth of the tree width
    // (inside clusters, when the branches between them are shortened).
    let bar = nice(scale_depth / 5.0);
    let bar_px = bar * sx_in;
    let y_bar = HEADER + ROW * n_tips as f64 + 24.0;
    let unit = if bar == 1.0 { "SNP" } else { "SNPs" };
    let scale = format!(
        "<path d=\"M{x0:.1},{y:.1}H{x1:.1}M{x0:.1},{ya:.1}V{yb:.1}M{x1:.1},{ya:.1}V{yb:.1}\" \
         stroke=\"{INK}\" stroke-width=\"{STROKE}\" fill=\"none\"/>\
         <text x=\"{xc:.1}\" y=\"{yt:.1}\" text-anchor=\"middle\" fill=\"{INK}\">{b} {unit}</text>",
        x0 = tree_left,
        x1 = tree_left + bar_px,
        y = y_bar,
        ya = y_bar - 4.0,
        yb = y_bar + 4.0,
        xc = tree_left + bar_px / 2.0,
        yt = y_bar + 17.0,
        b = fmt_snps(bar),
    );

    // One-line legend right of the scale bar.
    let mut legend = String::new();
    let mut lx = tree_left + bar_px.max(60.0) + 36.0;
    let ly = y_bar + 4.0;
    legend_item(
        &mut legend,
        &mut lx,
        ly,
        Swatch::None,
        "+n: SNPs on the genome's own branch (+0 = identical to its sister)",
    );
    if broken {
        legend_item(
            &mut legend,
            &mut lx,
            ly,
            Swatch::Break(INK),
            "shortened branch, length written on it",
        );
    }
    if l.nodes.iter().any(|n| n.estimate && n.x > n.x_parent) {
        legend_item(&mut legend, &mut lx, ly, Swatch::Line(ESTIMATE), "estimate");
    }
    if let Some(th) = band_threshold {
        legend_item(
            &mut legend,
            &mut lx,
            ly,
            Swatch::Band(BANDS[0]),
            &format!("SNP group ≤{th}"),
        );
    }

    let content_right = if r.clusters.is_empty() {
        labels_end
    } else {
        strip_x + 11.0 + cluster_label_w
    };
    let width = (content_right + MARGIN).max(lx + MARGIN).ceil();
    let height = (HEADER + ROW * n_tips as f64 + FOOTER).ceil();

    let n_multi = r.clusters.iter().filter(|c| c.size() >= 2).count();
    let subtitle = format!(
        "{n_tips} genomes · {nc} cluster{s} ({n_multi} with 2 or more genomes) · branch lengths in SNPs{bb} · Infinity-Win {v}",
        nc = r.clusters.len(),
        s = if r.clusters.len() == 1 { "" } else { "s" },
        bb = if r.backbone_is_exact || r.clusters.len() <= 1 {
            ""
        } else {
            ", between clusters estimated"
        },
        v = crate::VERSION,
    );
    let mut svg = String::new();
    let _ = writeln!(svg, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" \
         viewBox=\"0 0 {width} {height}\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"{FONT}\">"
    );
    let _ = writeln!(
        svg,
        "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>"
    );
    let _ = writeln!(
        svg,
        "<text x=\"{MARGIN}\" y=\"26\" font-size=\"15\" font-weight=\"bold\" fill=\"{INK}\">Core-genome SNP tree</text>"
    );
    let _ = writeln!(
        svg,
        "<text x=\"{MARGIN}\" y=\"44\" fill=\"{MUTED}\">{}</text>",
        super::xml_escape(&subtitle)
    );
    if n_tips == 0 {
        let _ = writeln!(
            svg,
            "<text x=\"{MARGIN}\" y=\"{}\" fill=\"{MUTED}\">No tree: no genomes were clustered.</text>",
            HEADER + ROW
        );
    } else {
        let _ = writeln!(svg, "<g id=\"snp-groups\">\n{bands}</g>");
        let _ = writeln!(svg, "<g id=\"clusters\">\n{strips}</g>");
        let _ = writeln!(
            svg,
            "<g id=\"branches\" fill=\"none\" stroke-width=\"{STROKE}\" stroke-linecap=\"butt\" stroke-linejoin=\"miter\">"
        );
        if !est_d.is_empty() {
            let _ = writeln!(svg, "<path stroke=\"{ESTIMATE}\" d=\"{est_d}\"/>");
        }
        if !exact_d.is_empty() {
            let _ = writeln!(svg, "<path stroke=\"{INK}\" d=\"{exact_d}\"/>");
        }
        let _ = writeln!(svg, "</g>");
        if !breaks.is_empty() {
            let _ = writeln!(svg, "<g id=\"breaks\">\n{breaks}</g>");
        }
        let _ = writeln!(svg, "<g id=\"tips\" fill=\"{INK}\">\n{tip_text}</g>");
        let _ = writeln!(svg, "<g id=\"scale\">{scale}</g>");
        let _ = writeln!(svg, "<g id=\"legend\">{legend}</g>");
    }
    svg.push_str("</svg>\n");
    svg
}

enum Swatch {
    None,
    Line(&'static str),
    Break(&'static str),
    Band(&'static str),
}

/// Append one legend entry at `lx` and move `lx` past it.
fn legend_item(out: &mut String, lx: &mut f64, ly: f64, swatch: Swatch, text: &str) {
    let tx = match swatch {
        Swatch::None => *lx,
        Swatch::Line(col) => {
            let _ = write!(
                out,
                "<path d=\"M{:.1},{:.1}H{:.1}\" stroke=\"{col}\" stroke-width=\"{STROKE}\"/>",
                *lx,
                ly - 4.0,
                *lx + 16.0
            );
            *lx + 22.0
        }
        Swatch::Break(col) => {
            let (x, y) = (*lx, ly - 4.0);
            let _ = write!(
                out,
                "<path d=\"M{:.1},{y:.1}H{:.1}M{:.1},{y:.1}H{:.1}M{:.1},{:.1}L{:.1},{:.1}M{:.1},{:.1}L{:.1},{:.1}\" \
                 stroke=\"{col}\" stroke-width=\"{STROKE}\" fill=\"none\"/>",
                x,
                x + 6.0,
                x + 10.0,
                x + 16.0,
                x + 3.0,
                y + 4.5,
                x + 7.0,
                y - 4.5,
                x + 9.0,
                y + 4.5,
                x + 13.0,
                y - 4.5,
            );
            *lx + 22.0
        }
        Swatch::Band(col) => {
            let _ = write!(
                out,
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"16\" height=\"10\" rx=\"2\" fill=\"{col}\"/>",
                *lx,
                ly - 9.0
            );
            *lx + 22.0
        }
    };
    let _ = write!(
        out,
        "<text x=\"{tx:.1}\" y=\"{ly:.1}\" fill=\"{MUTED}\">{}</text>",
        super::xml_escape(text)
    );
    *lx = tx + text_w(text) + 20.0;
}

/// A round number near `x`: 1, 2 or 5 times a power of ten, at least 1.
pub(super) fn nice(x: f64) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::parse_newick;

    fn laid_out(r: &RunResults) -> Layout {
        let mut cluster_of = BTreeMap::new();
        for c in &r.clusters {
            for &m in &c.members {
                cluster_of.insert(name(r, m), c.id);
            }
        }
        let mut l = layout(&r.full_tree_snps);
        classify(&mut l, r, &cluster_of);
        place(&mut l, MARGIN + ROOT_STUB, 10.0, 0.01);
        l
    }

    #[test]
    fn connectors_span_first_to_last_child_and_meet_the_children() {
        let r = crate::export::fixture();
        let l = laid_out(&r);
        assert_eq!(l.tips, ["a0", "a1", "a2", "a3", "b0", "b1"]);
        for n in &l.nodes {
            if let (Some(&a), Some(&b)) = (n.children.first(), n.children.last()) {
                assert_eq!(n.y, (l.nodes[a].y + l.nodes[b].y) / 2.0);
                assert_eq!((n.first, n.last), (l.nodes[a].first, l.nodes[b].last));
                for &c in &n.children {
                    assert_eq!(l.nodes[c].x_parent, n.x, "a child starts at its parent");
                    assert!(l.nodes[c].x >= n.x);
                }
            }
        }
    }

    #[test]
    fn branches_are_classified_inside_and_between_clusters() {
        let r = crate::export::fixture();
        let l = laid_out(&r);
        let leaf = |name: &str| {
            l.nodes
                .iter()
                .find(|n| n.leaf.as_deref() == Some(name))
                .unwrap()
        };
        // Inside cluster 0: exact SNPs.
        assert!(leaf("a2").within && !leaf("a2").estimate);
        // The branch above the whole of cluster 1 is the backbone's, which is
        // estimated in this run.
        let top_b = l
            .nodes
            .iter()
            .find(|n| n.first == 4 && n.last == 5)
            .unwrap();
        assert!(!top_b.within && top_b.estimate);
    }

    #[test]
    fn tips_carry_their_own_branch_in_snps() {
        let svg = tree_svg(&crate::export::fixture());
        assert!(svg.contains(">a0<tspan"), "{svg}");
        assert!(svg.contains(">+0</tspan>"), "a0 is identical to its sister");
        assert!(svg.contains(">+6</tspan>"), "a3 sits 6 SNPs out");
        assert!(svg.contains("Cluster 0 (n=4)"));
        assert!(svg.contains("Cluster 1 (n=2)"));
        assert!(svg.contains("stroke-linecap=\"butt\""));
        assert!(svg.starts_with("<?xml"));
        assert!(svg.trim_end().ends_with("</svg>"));
    }

    #[test]
    fn a_long_backbone_is_shortened_with_its_length_written_on_it() {
        let svg = tree_svg(&crate::export::fixture());
        assert!(svg.contains("id=\"breaks\""));
        assert!(svg.contains("≈5,000 SNPs"));
        assert!(svg.contains("shortened branch"));
        // The scale bar measures the inside of the clusters.
        assert!(svg.contains(">1 SNP</text>"), "{svg}");
    }

    #[test]
    fn a_short_exact_tree_is_drawn_to_scale() {
        let mut r = crate::export::fixture();
        r.full_tree_snps = parse_newick("(((a0:0,a1:0):1,a2:2,a3:6):4,(b0:1,b1:2):4);").unwrap();
        r.backbone_is_exact = true;
        let svg = tree_svg(&r);
        assert!(!svg.contains("id=\"breaks\""));
        assert!(!svg.contains(ESTIMATE), "nothing is an estimate");
        assert!(svg.contains(">2 SNPs</text>"), "{svg}");
    }

    #[test]
    fn snp_groups_of_the_smallest_threshold_are_shaded() {
        let svg = tree_svg(&crate::export::fixture());
        let bands = svg
            .split("<g id=\"snp-groups\">")
            .nth(1)
            .and_then(|s| s.split("</g>").next())
            .unwrap();
        // a0, a1, a2 within 5 SNPs; b0, b1 within 5 SNPs; a3 alone.
        assert_eq!(bands.matches("<rect").count(), 2, "{bands}");
        assert!(svg.contains("SNP group ≤5"));
    }

    #[test]
    fn an_empty_run_still_draws_a_page() {
        let svg = tree_svg(&RunResults::default());
        assert!(svg.contains("No tree"));
        assert!(svg.trim_end().ends_with("</svg>"));
    }

    #[test]
    fn snp_counts_read_well() {
        assert_eq!(fmt_snps(3.0), "3");
        assert_eq!(fmt_snps(0.01), "0");
        assert_eq!(fmt_snps(2.5), "2.5");
        assert_eq!(fmt_snps(40591.5), "40,592");
        assert_eq!(group_thousands(1_234_567), "1,234,567");
        assert_eq!(group_thousands(999), "999");
    }

    #[test]
    fn runs_group_consecutive_rows() {
        let keys = [Some(1), Some(1), None, Some(2), Some(1)];
        assert_eq!(runs(&keys), vec![(1, 0, 1), (2, 3, 3), (1, 4, 4)]);
    }
}
