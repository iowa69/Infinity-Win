//! The tree as a publication-style rectangular phylogram, built once as a
//! list of drawing primitives in pixels ([`Prim`]).
//!
//! One builder serves every picture of the tree: the main tree panel (in the
//! window's theme), the white tree window, and the SVG export
//! ([`crate::svg::from_scene`] writes the same list), so what is saved is what
//! is on screen. Pure data, no `egui::Ui`; [`paint`] is the only function that
//! touches a painter.
//!
//! # The conventions (ggtree, iTOL and FigTree defaults)
//!
//! * One line width for every branch, square elbows, and each vertical
//!   connector spans exactly from its first to its last child: nothing runs
//!   past a node or a tip.
//! * Exact SNP branches in ink, estimated ones (between clusters, or a cluster
//!   still being refined) in a lighter grey, both solid.
//! * Tip labels right after the tip, followed by the SNPs on that genome's own
//!   branch: `isolate_12  +3`, with `≈` when the number is an estimate.
//! * A collapsed cluster is a filled wedge with `c4 · 30 genomes · exact`.
//! * SNP groups are soft rounded bands behind their tips, labelled
//!   `G3 · 4 genomes ≤5 SNPs`.
//! * A footer under the tree holds the SNP scale bar(s) and a legend of only
//!   what the picture actually shows.
//!
//! Honesty rules carry over from the rest of the program (section 0 rule 8):
//! every estimated number carries `≈`, and grey branches say "estimate".

use std::collections::{BTreeMap, HashMap, HashSet};

use egui::{Color32, Pos2, Rect, Vec2};

use crate::fmt;
use crate::layout::{self, NodeKind, TreeLayout};

/// Every branch is drawn this wide, in points.
pub const BRANCH_W: f32 = 1.3;
/// Gap between a tip and its label.
const TIP_GAP: f32 = 5.0;
/// Below this row height the tip labels are dropped (they would overlap).
pub const MIN_ROW_FOR_LABELS: f32 = 6.0;
/// Line height of the footer (scale bars and legend).
const FOOTER_LINE: f32 = 22.0;
const FOOTER_PAD: f32 = 8.0;

// ---------------------------------------------------------------------------
// Colours
// ---------------------------------------------------------------------------

/// The colours one picture is drawn in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ink {
    pub dark: bool,
    pub paper: Color32,
    /// Exact branches.
    pub exact: Color32,
    /// Estimated branches.
    pub estimate: Color32,
    pub text: Color32,
    /// Tip SNP numbers, branch labels, secondary text.
    pub text_weak: Color32,
    /// Leader lines, frames, the footer rule.
    pub hairline: Color32,
    /// The band behind a selected genome's row.
    pub select_band: Color32,
    /// Selection marks drawn on top: ring, badge.
    pub select_ink: Color32,
    /// The ring around the search hit.
    pub search: Color32,
}

/// Paper: the white figure of the tree window and the SVG export.
pub const PAPER: Color32 = Color32::WHITE;

impl Ink {
    /// White paper, near-black ink: the figure.
    pub const PAPER: Ink = Ink {
        dark: false,
        paper: PAPER,
        exact: Color32::from_rgb(28, 28, 30),
        estimate: Color32::from_rgb(170, 172, 176),
        text: Color32::from_rgb(28, 28, 30),
        text_weak: Color32::from_rgb(112, 114, 120),
        hairline: Color32::from_rgb(214, 216, 220),
        select_band: Color32::from_rgba_premultiplied(120, 100, 42, 120),
        select_ink: Color32::from_rgb(196, 98, 0),
        search: Color32::from_rgb(214, 120, 0),
    };

    /// The tree panel, in the window's theme.
    pub fn screen(p: &crate::theme::Palette) -> Ink {
        if p.dark {
            Ink {
                dark: true,
                paper: p.canvas,
                exact: Color32::from_rgb(222, 228, 236),
                estimate: Color32::from_rgb(104, 116, 134),
                text: Color32::from_rgb(214, 221, 230),
                text_weak: Color32::from_rgb(132, 146, 164),
                hairline: p.border,
                select_band: crate::theme::with_alpha(p.accent, 44),
                select_ink: p.accent,
                search: Color32::from_rgb(250, 200, 60),
            }
        } else {
            Ink {
                dark: false,
                paper: p.canvas,
                exact: Color32::from_rgb(30, 34, 40),
                estimate: Color32::from_rgb(160, 168, 180),
                text: Color32::from_rgb(30, 36, 44),
                text_weak: Color32::from_rgb(108, 118, 132),
                hairline: p.border,
                select_band: crate::theme::with_alpha(p.accent, 40),
                select_ink: p.accent,
                search: Color32::from_rgb(214, 120, 0),
            }
        }
    }

    pub fn branch(&self, exact: bool) -> Color32 {
        if exact {
            self.exact
        } else {
            self.estimate
        }
    }

    /// The soft fill of SNP group `number` (1-based).
    pub fn group_fill(&self, number: usize) -> Color32 {
        let c = group_hue(number);
        crate::theme::with_alpha(c, if self.dark { 72 } else { 120 })
    }

    /// The label colour of SNP group `number`: the hue, deep enough to read.
    pub fn group_text(&self, number: usize) -> Color32 {
        let c = group_hue(number);
        if self.dark {
            crate::theme::mix(c, Color32::WHITE, 0.25)
        } else {
            crate::theme::mix(c, Color32::BLACK, 0.55)
        }
    }

    /// A cluster's colour: the theme's own (the paper uses the light one,
    /// which reads as text on white).
    pub fn cluster(&self, id: u32) -> Color32 {
        crate::colors::cluster_color(id, self.dark)
    }
}

/// Pale, clearly distinct hues for SNP groups (after ColorBrewer Set3).
const GROUP_HUES: [Color32; 10] = [
    Color32::from_rgb(141, 211, 199),
    Color32::from_rgb(251, 128, 114),
    Color32::from_rgb(128, 177, 211),
    Color32::from_rgb(253, 180, 98),
    Color32::from_rgb(179, 222, 105),
    Color32::from_rgb(252, 205, 229),
    Color32::from_rgb(190, 186, 218),
    Color32::from_rgb(255, 237, 111),
    Color32::from_rgb(188, 128, 189),
    Color32::from_rgb(204, 235, 197),
];

pub fn group_hue(number: usize) -> Color32 {
    GROUP_HUES[number.saturating_sub(1) % GROUP_HUES.len()]
}

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------

/// One thing to draw, in pixels.
#[derive(Debug, Clone, PartialEq)]
pub enum Prim {
    Line {
        a: Pos2,
        b: Pos2,
        width: f32,
        color: Color32,
    },
    Circle {
        c: Pos2,
        r: f32,
        fill: Color32,
        stroke: Option<(f32, Color32)>,
    },
    Rect {
        r: Rect,
        fill: Color32,
        rounding: f32,
    },
    Poly {
        pts: Vec<Pos2>,
        fill: Color32,
        stroke: (f32, Color32),
    },
    /// Text whose left edge (or centre) is at `pos.x`, centred vertically on
    /// `pos.y`, optionally followed by a smaller second part in another colour
    /// (the tip's SNP number).
    Text {
        pos: Pos2,
        align: HAlign,
        text: String,
        size: f32,
        color: Color32,
        bold: bool,
        tail: Option<Tail>,
    },
}

/// The second part of a tip label: `+3` after the name.
#[derive(Debug, Clone, PartialEq)]
pub struct Tail {
    pub text: String,
    pub size: f32,
    pub color: Color32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HAlign {
    Left,
    Center,
}

fn text(pos: Pos2, align: HAlign, s: impl Into<String>, size: f32, color: Color32) -> Prim {
    Prim::Text {
        pos,
        align,
        text: s.into(),
        size,
        color,
        bold: false,
        tail: None,
    }
}

/// Estimated width of a text, the same in the window and in the SVG so the
/// two lay out alike.
pub fn text_width(text: &str, size: f32) -> f32 {
    text.chars().count() as f32 * size * 0.56
}

/// The gap between a tip label and its SNP number.
pub fn tail_gap(size: f32) -> f32 {
    size * 0.55
}

/// Label font size for a row height.
pub fn label_size(px_per_row: f32) -> f32 {
    (px_per_row * 0.7).clamp(7.5, 13.0)
}

// ---------------------------------------------------------------------------
// What to draw
// ---------------------------------------------------------------------------

/// What a click on the picture hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hit {
    Genome(u32),
    /// A collapsed wedge, a cluster bracket or an opened cluster's root.
    Cluster(u32),
    /// The "n selected" badge on a collapsed cluster.
    Badge(u32),
    /// A SNP group band's label: index into [`Marks::groups`].
    Group(usize),
}

/// One cluster, as far as the picture needs it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterBrief {
    pub size: usize,
    pub focus: u8,
    pub approximate: bool,
}

impl ClusterBrief {
    pub fn exact(&self) -> bool {
        self.focus >= fmt::FOCUS_EXACT && !self.approximate
    }
}

/// A SNP group to shade.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupBand {
    /// 1-based: the `3` of `G3`, and the colour.
    pub number: usize,
    pub members: Vec<u32>,
}

/// Everything about the picture besides the tree's shape.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Marks {
    pub genome_cluster: HashMap<u32, u32>,
    pub clusters: BTreeMap<u32, ClusterBrief>,
    pub selection: HashSet<u32>,
    /// Selected genomes per collapsed cluster, for its badge.
    pub selected_in: HashMap<u32, usize>,
    pub groups: Vec<GroupBand>,
    /// The N of "≤ N SNPs".
    pub group_threshold: u32,
    pub search_hit: Option<u32>,
}

/// The shape and the marks: what the tree window keeps between frames.
#[derive(Debug, Clone, Default)]
pub struct Figure {
    pub layout: TreeLayout,
    pub marks: Marks,
}

/// Decorations the user can switch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PubOptions {
    /// The SNPs on each genome's own branch, after its name.
    pub tip_snps: bool,
    /// SNP counts above internal branches (only where at least 1).
    pub branch_labels: bool,
    pub cluster_brackets: bool,
    pub legend: bool,
    /// Tip labels in one column, with faint leader lines.
    pub aligned_tips: bool,
    /// A cluster-coloured strip before each tip label.
    pub color_strip: bool,
    /// Shade the SNP groups.
    pub groups: bool,
}

impl Default for PubOptions {
    fn default() -> Self {
        Self {
            tip_snps: true,
            branch_labels: false,
            cluster_brackets: true,
            legend: true,
            aligned_tips: false,
            color_strip: false,
            groups: true,
        }
    }
}

/// Layout units to pixels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geometry {
    pub origin: Pos2,
    pub px_per_row: f32,
    pub px_per_unit: f32,
    /// What is visible: anything wholly outside is skipped.
    pub clip: Rect,
}

impl Geometry {
    pub fn screen(&self, x: f64, y: f64) -> Pos2 {
        Pos2::new(
            self.origin.x + x as f32 * self.px_per_unit,
            self.origin.y + y as f32 * self.px_per_row,
        )
    }

    fn row_visible(&self, y: f32, slack: f32) -> bool {
        y >= self.clip.top() - slack && y <= self.clip.bottom() + slack
    }
}

/// What a scene actually drew, for a legend of only that.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Used {
    pub exact: bool,
    pub estimate: bool,
    pub groups: bool,
    pub selection: bool,
    pub wedges: bool,
    pub tip_snps: bool,
}

/// A built picture.
#[derive(Debug, Clone, Default)]
pub struct Scene {
    pub prims: Vec<Prim>,
    pub hits: Vec<(Hit, Rect)>,
    pub used: Used,
}

fn is_exact(focus: u8, approximate: bool) -> bool {
    focus >= fmt::FOCUS_EXACT && !approximate
}

// ---------------------------------------------------------------------------
// Tip text and the columns right of the tree
// ---------------------------------------------------------------------------

/// Children per node, and the y range of each node's children (in layout
/// rows, at `positions` when the tree is moving).
struct Kids {
    count: Vec<u32>,
    lo: Vec<f64>,
    hi: Vec<f64>,
}

fn kids(l: &TreeLayout, pos: &dyn Fn(usize) -> (f64, f64)) -> Kids {
    let n = l.nodes.len();
    let mut k = Kids {
        count: vec![0; n],
        lo: vec![f64::INFINITY; n],
        hi: vec![f64::NEG_INFINITY; n],
    };
    for node in &l.nodes {
        if let Some(p) = node.parent {
            let y = pos(node.idx).1;
            k.count[p] += 1;
            k.lo[p] = k.lo[p].min(y);
            k.hi[p] = k.hi[p].max(y);
        }
    }
    k
}

/// The SNPs on a tip's own branch: the drawn segment that ends at the tip,
/// through any single-child nodes (a one-genome cluster hangs off the
/// backbone that way). Returns the count and whether it is an estimate.
fn own_branch(l: &TreeLayout, k: &Kids, leaf: usize) -> (f64, bool) {
    let mut len = 0.0;
    let mut approx = false;
    let mut i = leaf;
    loop {
        let n = &l.nodes[i];
        len += n.len.max(0.0);
        approx |= n.edge_approximate;
        match n.parent {
            Some(p) if k.count[p] == 1 && l.nodes[p].parent.is_some() => i = p,
            _ => break,
        }
    }
    (len, approx)
}

/// `+3`, `+2.5`, `≈+1 200`.
pub fn tip_number(len: f64, approximate: bool) -> String {
    format!("{}+{}", if approximate { "≈" } else { "" }, fmt::snps(len))
}

/// "exact" or "≈ estimate", for a cluster's label.
pub fn cluster_state(c: Option<&ClusterBrief>) -> &'static str {
    match c {
        Some(c) if c.exact() => "exact",
        _ => "≈ estimate",
    }
}

fn wedge_label(id: u32, c: Option<&ClusterBrief>, members: usize) -> String {
    format!(
        "c{id} · {members} genome{} · {}",
        if members == 1 { "" } else { "s" },
        cluster_state(c)
    )
}

fn group_label(b: &GroupBand, threshold: u32) -> String {
    format!(
        "G{} · {} genomes ≤{} SNP{}",
        b.number,
        b.members.len(),
        threshold,
        if threshold == 1 { "" } else { "s" }
    )
}

/// Width of a collapsed wedge, in pixels.
fn wedge_px(g: &Geometry) -> f32 {
    (layout::COLLAPSED_SPREAD as f32 * g.px_per_unit).clamp(14.0, 40.0)
}

/// The text after a tip, and its width, for the label-room estimates.
fn tip_texts(
    l: &TreeLayout,
    k: &Kids,
    m: &Marks,
    o: &PubOptions,
    i: usize,
    size: f32,
) -> (String, Option<String>, f32) {
    let n = &l.nodes[i];
    match &n.kind {
        NodeKind::Genome { name, .. } => {
            let tail = o.tip_snps.then(|| {
                let (len, approx) = own_branch(l, k, i);
                tip_number(len, approx)
            });
            let w = text_width(name, size)
                + tail
                    .as_ref()
                    .map(|t| tail_gap(size) + text_width(t, size * 0.86))
                    .unwrap_or(0.0);
            (name.clone(), tail, w)
        }
        NodeKind::ClusterCollapsed { id, members, .. } => {
            let s = wedge_label(*id, m.clusters.get(id), *members);
            let badge = if m.selected_in.contains_key(id) {
                70.0
            } else {
                0.0
            };
            let w = text_width(&s, size) + badge;
            (s, None, w)
        }
        _ => (String::new(), None, 0.0),
    }
}

/// Where the tip labels start: at the tip, or in one aligned column.
fn label_x(l: &TreeLayout, g: &Geometry, o: &PubOptions, tip_x: f32, collapsed: bool) -> f32 {
    let wedge = if collapsed { wedge_px(g) } else { 0.0 };
    let x = if o.aligned_tips {
        g.screen(l.x_max, 0.0).x + TIP_GAP + 3.0
    } else {
        tip_x + wedge + TIP_GAP
    };
    x + if o.color_strip { 10.0 } else { 0.0 }
}

/// Groups that have at least one tip in this layout (opened clusters only).
fn visible_groups(l: &TreeLayout, m: &Marks) -> Vec<usize> {
    if m.groups.is_empty() {
        return Vec::new();
    }
    let shown: HashSet<u32> = l
        .leaves
        .iter()
        .filter_map(|&i| match &l.nodes[i].kind {
            NodeKind::Genome { id, .. } => Some(*id),
            _ => None,
        })
        .collect();
    (0..m.groups.len())
        .filter(|&i| m.groups[i].members.iter().any(|g| shown.contains(g)))
        .collect()
}

/// Room the picture needs right of the deepest tip at this row height: the
/// widest tip label (with its SNP number), the group labels and the cluster
/// brackets. Used to fit the tree into a panel.
pub fn label_room(l: &TreeLayout, m: &Marks, o: &PubOptions, px_per_row: f32) -> f32 {
    let size = label_size(px_per_row);
    let pos = |i: usize| (l.nodes[i].x, l.nodes[i].y);
    let k = kids(l, &pos);
    let labels = if px_per_row >= MIN_ROW_FOR_LABELS {
        l.leaves
            .iter()
            .map(|&i| {
                let w = tip_texts(l, &k, m, o, i, size).2;
                let wedge = if matches!(l.nodes[i].kind, NodeKind::ClusterCollapsed { .. }) {
                    30.0
                } else {
                    0.0
                };
                w + wedge
            })
            .fold(0.0f32, f32::max)
            + TIP_GAP
            + if o.color_strip { 10.0 } else { 0.0 }
    } else {
        6.0
    };
    labels + group_column(l, m, o, px_per_row) + bracket_column(l, o, px_per_row) + 8.0
}

/// Width of the group-label column (0 when no group band shows).
fn group_column(l: &TreeLayout, m: &Marks, o: &PubOptions, px_per_row: f32) -> f32 {
    if !o.groups {
        return 0.0;
    }
    let size = group_size(px_per_row);
    let w = visible_groups(l, m)
        .into_iter()
        .map(|i| text_width(&group_label(&m.groups[i], m.group_threshold), size))
        .fold(0.0f32, f32::max);
    if w > 0.0 {
        w + 18.0
    } else {
        0.0
    }
}

fn group_size(px_per_row: f32) -> f32 {
    (label_size(px_per_row) * 0.9).clamp(7.5, 11.0)
}

fn bracket_size(px_per_row: f32) -> f32 {
    (label_size(px_per_row) * 0.95).clamp(9.0, 12.5)
}

/// Width of the cluster-bracket column (0 when no cluster is opened).
fn bracket_column(l: &TreeLayout, o: &PubOptions, px_per_row: f32) -> f32 {
    if !o.cluster_brackets {
        return 0.0;
    }
    let size = bracket_size(px_per_row);
    l.nodes
        .iter()
        .filter_map(|n| match &n.kind {
            NodeKind::ClusterRoot { id, members, .. } => {
                Some(text_width(&format!("c{id} · {members}"), size) + 20.0)
            }
            _ => None,
        })
        .fold(0.0f32, f32::max)
}

/// The lowest common ancestor of some nodes.
fn mrca(l: &TreeLayout, nodes: &[usize]) -> usize {
    let Some(&first) = nodes.first() else {
        return l.root;
    };
    let mut chain: Vec<usize> = Vec::new();
    let mut i = Some(first);
    while let Some(n) = i {
        chain.push(n);
        i = l.nodes[n].parent;
    }
    let at: HashMap<usize, usize> = chain.iter().enumerate().map(|(k, &n)| (n, k)).collect();
    let mut deepest = 0usize;
    for &n in &nodes[1..] {
        let mut j = Some(n);
        while let Some(x) = j {
            if let Some(&k) = at.get(&x) {
                deepest = deepest.max(k);
                break;
            }
            j = l.nodes[x].parent;
        }
    }
    chain[deepest.min(chain.len() - 1)]
}

// ---------------------------------------------------------------------------
// The scene
// ---------------------------------------------------------------------------

/// The tree: selection and group bands, branches, tips with their labels and
/// SNP numbers, wedges, branch labels and cluster brackets. `positions`
/// overrides the node positions while the tree is moving (parallel to
/// `l.nodes`).
pub fn scene(
    l: &TreeLayout,
    m: &Marks,
    g: &Geometry,
    o: &PubOptions,
    ink: &Ink,
    positions: Option<&[(f64, f64)]>,
) -> Scene {
    let pos = |i: usize| -> (f64, f64) {
        positions
            .and_then(|p| p.get(i).copied())
            .unwrap_or((l.nodes[i].x, l.nodes[i].y))
    };
    let px = |i: usize| -> Pos2 {
        let (x, y) = pos(i);
        g.screen(x, y)
    };
    let k = kids(l, &pos);
    let row = g.px_per_row;
    let show_labels = row >= MIN_ROW_FOR_LABELS;
    let size = label_size(row);
    let mut s = Scene {
        prims: Vec::with_capacity(l.nodes.len() * 3),
        ..Default::default()
    };

    // The columns right of the tree. Measured over every tip, not just the
    // visible ones, so they do not shift while the user pans.
    let mut labels_right = f32::MIN;
    for &i in &l.leaves {
        let collapsed = matches!(l.nodes[i].kind, NodeKind::ClusterCollapsed { .. });
        let x0 = label_x(l, g, o, px(i).x, collapsed);
        let w = if show_labels {
            tip_texts(l, &k, m, o, i, size).2
        } else {
            0.0
        };
        labels_right = labels_right.max(x0 + w);
    }
    if labels_right == f32::MIN {
        labels_right = g.origin.x;
    }
    let group_x = labels_right + 10.0;
    let bracket_x = labels_right + group_column(l, m, o, row) + 10.0;

    // --- selection bands, behind everything ---------------------------------
    for &i in &l.leaves {
        let y = px(i).y;
        if !g.row_visible(y, row) {
            continue;
        }
        let strength = match &l.nodes[i].kind {
            NodeKind::Genome { id, .. } if m.selection.contains(id) => Some(1.0),
            NodeKind::ClusterCollapsed { id, .. } if m.selected_in.contains_key(id) => Some(0.6),
            _ => None,
        };
        if let Some(f) = strength {
            let h = row.max(2.0);
            let [r, gg, b, a] = ink.select_band.to_srgba_unmultiplied();
            s.prims.push(Prim::Rect {
                r: Rect::from_min_size(
                    Pos2::new(g.clip.left(), y - h / 2.0),
                    Vec2::new(g.clip.width(), h),
                ),
                fill: Color32::from_rgba_unmultiplied(r, gg, b, (f32::from(a) * f) as u8),
                rounding: 0.0,
            });
            s.used.selection = true;
        }
    }

    // --- SNP group bands -----------------------------------------------------
    if o.groups && !m.groups.is_empty() {
        let leaf_of: HashMap<u32, usize> = l
            .leaves
            .iter()
            .filter_map(|&i| match &l.nodes[i].kind {
                NodeKind::Genome { id, .. } => Some((*id, i)),
                _ => None,
            })
            .collect();
        let gsize = group_size(row);
        let right = bracket_x - 6.0;
        for (gi, band) in m.groups.iter().enumerate() {
            let mut tips: Vec<usize> = band
                .members
                .iter()
                .filter_map(|g| leaf_of.get(g).copied())
                .collect();
            if tips.is_empty() {
                continue;
            }
            tips.sort_by(|&a, &b| l.nodes[a].y.total_cmp(&l.nodes[b].y));
            // Contiguous runs of rows: a group the tree splits gets one band
            // per run, labelled once (on its tallest run).
            let mut runs: Vec<Vec<usize>> = Vec::new();
            for t in tips {
                match runs.last_mut() {
                    Some(r)
                        if (l.nodes[t].y - l.nodes[*r.last().unwrap()].y - 1.0).abs() < 1e-6 =>
                    {
                        r.push(t)
                    }
                    _ => runs.push(vec![t]),
                }
            }
            let tallest = (0..runs.len()).max_by_key(|&r| runs[r].len()).unwrap_or(0);
            for (ri, run) in runs.iter().enumerate() {
                let top = px(run[0]).y - row * 0.48;
                let bottom = px(*run.last().unwrap()).y + row * 0.48;
                if bottom < g.clip.top() || top > g.clip.bottom() {
                    continue;
                }
                let anchor = mrca(l, run);
                let left = if run.len() == 1 {
                    px(run[0]).x - 6.0
                } else {
                    px(anchor).x - 5.0
                };
                let r = Rect::from_min_max(Pos2::new(left, top), Pos2::new(right, bottom));
                s.prims.push(Prim::Rect {
                    r,
                    fill: ink.group_fill(band.number),
                    rounding: (row * 0.45).clamp(2.0, 7.0),
                });
                s.used.groups = true;
                // Every run is labelled, the tallest in full: a group the
                // tree splits must not leave an unexplained band behind.
                if r.height() >= gsize * 1.1 && show_labels {
                    let t = if ri == tallest {
                        group_label(band, m.group_threshold)
                    } else {
                        format!("G{} (continued)", band.number)
                    };
                    let w = text_width(&t, gsize);
                    let at = Pos2::new(group_x, (top + bottom) / 2.0);
                    s.prims.push(Prim::Text {
                        pos: at,
                        align: HAlign::Left,
                        text: t,
                        size: gsize,
                        color: ink.group_text(band.number),
                        bold: true,
                        tail: None,
                    });
                    s.hits.push((
                        Hit::Group(gi),
                        Rect::from_min_size(
                            Pos2::new(at.x - 2.0, at.y - gsize * 0.7),
                            Vec2::new(w + 4.0, gsize * 1.4),
                        ),
                    ));
                }
            }
        }
    }

    // --- branches --------------------------------------------------------------
    // Horizontal: from the parent's x to the node, at the node's height.
    for n in &l.nodes {
        let Some(p) = n.parent else { continue };
        let c = px(n.idx);
        if !g.row_visible(c.y, 2.0) {
            continue;
        }
        let a = Pos2::new(px(p).x, c.y);
        if (c.x - a.x).abs() < 0.05 {
            continue;
        }
        let exact = is_exact(n.edge_focus, n.edge_approximate);
        mark(&mut s.used, exact);
        s.prims.push(Prim::Line {
            a,
            b: c,
            width: BRANCH_W,
            color: ink.branch(exact),
        });
    }
    // Vertical: exactly from the first child to the last, plus half a line
    // width at each end so the elbows are square.
    for n in &l.nodes {
        if k.count[n.idx] < 2 {
            continue;
        }
        let x = px(n.idx).x;
        let top = g.screen(0.0, k.lo[n.idx]).y - BRANCH_W / 2.0;
        let bottom = g.screen(0.0, k.hi[n.idx]).y + BRANCH_W / 2.0;
        if bottom < g.clip.top() - 2.0 || top > g.clip.bottom() + 2.0 || bottom - top < 0.5 {
            continue;
        }
        let exact = is_exact(n.focus, n.approximate);
        mark(&mut s.used, exact);
        s.prims.push(Prim::Line {
            a: Pos2::new(x, top),
            b: Pos2::new(x, bottom),
            width: BRANCH_W,
            color: ink.branch(exact),
        });
    }

    // --- SNP counts on internal branches -----------------------------------------
    if o.branch_labels && row >= 10.0 {
        let fs = (row * 0.55).clamp(8.0, 10.5);
        for n in &l.nodes {
            let Some(p) = n.parent else { continue };
            if k.count[n.idx] == 0 || n.len < 0.5 || matches!(n.kind, NodeKind::ClusterRoot { .. })
            {
                continue;
            }
            let b = px(n.idx);
            let a = Pos2::new(px(p).x, b.y);
            if !g.row_visible(b.y, 2.0) {
                continue;
            }
            let t = format!(
                "{}{}",
                if n.edge_approximate { "≈" } else { "" },
                fmt::snps(n.len)
            );
            if b.x - a.x < text_width(&t, fs) + 6.0 {
                continue;
            }
            s.prims.push(text(
                Pos2::new((a.x + b.x) / 2.0, b.y - fs * 0.62 - 1.0),
                HAlign::Center,
                t,
                fs,
                ink.text_weak,
            ));
        }
    }

    // --- tips, wedges, cluster roots -------------------------------------------
    for &i in &l.leaves {
        let n = &l.nodes[i];
        let tip = px(i);
        if !g.row_visible(tip.y, row + 4.0) {
            continue;
        }
        match &n.kind {
            NodeKind::Genome { id, .. } => {
                let selected = m.selection.contains(id);
                let colour = m
                    .genome_cluster
                    .get(id)
                    .map(|&c| ink.cluster(c))
                    .unwrap_or(ink.text_weak);
                let mut hit = Rect::from_center_size(tip, Vec2::splat(8.0));
                if show_labels || selected {
                    s.prims.push(Prim::Circle {
                        c: tip,
                        r: if selected { 2.6 } else { 1.9 },
                        fill: colour,
                        stroke: selected.then_some((1.4, ink.select_ink)),
                    });
                }
                if Some(*id) == m.search_hit {
                    s.prims.push(Prim::Circle {
                        c: tip,
                        r: 5.0,
                        fill: Color32::TRANSPARENT,
                        stroke: Some((1.4, ink.search)),
                    });
                }
                if show_labels {
                    let x = label_x(l, g, o, tip.x, false);
                    if o.aligned_tips && x - tip.x > 12.0 {
                        s.prims.push(Prim::Line {
                            a: Pos2::new(tip.x + 4.0, tip.y),
                            b: Pos2::new(x - 4.0 - if o.color_strip { 10.0 } else { 0.0 }, tip.y),
                            width: 0.6,
                            color: ink.hairline,
                        });
                    }
                    if o.color_strip {
                        let h = (row * 0.8).max(2.0);
                        s.prims.push(Prim::Rect {
                            r: Rect::from_min_size(
                                Pos2::new(x - 10.0, tip.y - h / 2.0),
                                Vec2::new(6.0, h),
                            ),
                            fill: colour,
                            rounding: 1.0,
                        });
                    }
                    let (name, tail, w) = tip_texts(l, &k, m, o, i, size);
                    if tail.is_some() {
                        s.used.tip_snps = true;
                    }
                    s.prims.push(Prim::Text {
                        pos: Pos2::new(x, tip.y),
                        align: HAlign::Left,
                        text: name,
                        size,
                        color: ink.text,
                        bold: selected,
                        tail: tail.map(|t| Tail {
                            text: t,
                            size: size * 0.86,
                            color: ink.text_weak,
                        }),
                    });
                    hit = hit.union(Rect::from_min_size(
                        Pos2::new(x, tip.y - row / 2.0),
                        Vec2::new(w, row),
                    ));
                }
                s.hits.push((Hit::Genome(*id), hit));
            }
            NodeKind::ClusterCollapsed { id, members, .. } => {
                s.used.wedges = true;
                let half = (layout::collapsed_half_rows(*members) as f32 * row).max(2.0);
                let wpx = wedge_px(g);
                let pts = vec![
                    tip,
                    Pos2::new(tip.x + wpx, tip.y - half),
                    Pos2::new(tip.x + wpx, tip.y + half),
                ];
                let colour = ink.cluster(*id);
                let mut hit = Rect::from_points(&pts).expand(2.0);
                s.prims.push(Prim::Poly {
                    pts,
                    fill: crate::theme::with_alpha(colour, if ink.dark { 120 } else { 96 }),
                    stroke: (1.0, colour),
                });
                if show_labels {
                    let x = label_x(l, g, o, tip.x, true);
                    let t = wedge_label(*id, m.clusters.get(id), *members);
                    let w = text_width(&t, size);
                    s.prims
                        .push(text(Pos2::new(x, tip.y), HAlign::Left, t, size, ink.text));
                    hit = hit.union(Rect::from_min_size(
                        Pos2::new(x, tip.y - row / 2.0),
                        Vec2::new(w, row),
                    ));
                    if let Some(&k) = m.selected_in.get(id) {
                        let bx = x + w + 8.0;
                        let bs = size * 0.85;
                        let t = format!("{k} selected");
                        let bw = text_width(&t, bs) + 10.0;
                        let bh = (bs * 1.45).max(10.0);
                        let r =
                            Rect::from_min_size(Pos2::new(bx, tip.y - bh / 2.0), Vec2::new(bw, bh));
                        s.prims.push(Prim::Rect {
                            r,
                            fill: ink.select_ink,
                            rounding: bh / 2.0,
                        });
                        s.prims.push(Prim::Text {
                            pos: Pos2::new(bx + 5.0, tip.y),
                            align: HAlign::Left,
                            text: t,
                            size: bs,
                            color: ink.paper,
                            bold: true,
                            tail: None,
                        });
                        s.hits.push((Hit::Badge(*id), r));
                    }
                }
                s.hits.push((Hit::Cluster(*id), hit));
            }
            _ => {}
        }
    }
    // An opened cluster's root is a click target that closes it.
    for n in &l.nodes {
        if let NodeKind::ClusterRoot { id, .. } = &n.kind {
            let c = px(n.idx);
            if g.row_visible(c.y, 4.0) {
                s.hits.push((
                    Hit::Cluster(*id),
                    Rect::from_center_size(c, Vec2::splat(9.0)),
                ));
            }
        }
    }

    // --- cluster brackets ------------------------------------------------------------
    if o.cluster_brackets {
        let bs = bracket_size(row);
        for n in &l.nodes {
            let NodeKind::ClusterRoot { id, members, .. } = &n.kind else {
                continue;
            };
            let top = g.screen(0.0, n.row_min).y - row * 0.42;
            let bottom = g.screen(0.0, n.row_max).y + row * 0.42;
            if bottom < g.clip.top() || top > g.clip.bottom() {
                continue;
            }
            let colour = ink.cluster(*id);
            let x = bracket_x;
            s.prims.push(Prim::Line {
                a: Pos2::new(x, top),
                b: Pos2::new(x, bottom),
                width: 2.0,
                color: colour,
            });
            for yy in [top, bottom] {
                s.prims.push(Prim::Line {
                    a: Pos2::new(x - 4.0, yy),
                    b: Pos2::new(x + 1.0, yy),
                    width: 2.0,
                    color: colour,
                });
            }
            let brief = m.clusters.get(id);
            let label = format!("c{id} · {members}");
            let w = text_width(&label, bs);
            let ty = ((top + bottom) / 2.0).clamp(g.clip.top() + 8.0, g.clip.bottom() - 8.0);
            s.prims.push(Prim::Text {
                pos: Pos2::new(x + 7.0, ty),
                align: HAlign::Left,
                text: label,
                size: bs,
                color: colour,
                bold: true,
                tail: None,
            });
            if bottom - top > bs * 2.6 {
                s.prims.push(text(
                    Pos2::new(x + 7.0, ty + bs * 1.25),
                    HAlign::Left,
                    cluster_state(brief),
                    bs * 0.85,
                    ink.text_weak,
                ));
            }
            s.hits.push((
                Hit::Cluster(*id),
                Rect::from_min_max(Pos2::new(x - 5.0, top), Pos2::new(x + 10.0 + w, bottom)),
            ));
        }
    }
    s
}

fn mark(u: &mut Used, exact: bool) {
    if exact {
        u.exact = true;
    } else {
        u.estimate = true;
    }
}

// ---------------------------------------------------------------------------
// Footer: scale bars and legend
// ---------------------------------------------------------------------------

/// One item of the footer, laid out in a flowing row.
enum FooterItem {
    Bar { px: f32, exact: bool, label: String },
    Legend { mark: LegendMark, label: String },
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum LegendMark {
    Branch(bool),
    Group,
    Selected,
    Wedge,
    TipNumber,
}

const FOOTER_TEXT: f32 = 11.0;

impl FooterItem {
    fn width(&self) -> f32 {
        match self {
            FooterItem::Bar { px, label, .. } => px + 8.0 + text_width(label, FOOTER_TEXT),
            FooterItem::Legend { mark, label } => {
                let m = match mark {
                    LegendMark::TipNumber => text_width("+3", FOOTER_TEXT) + 6.0,
                    _ => 24.0,
                };
                m + text_width(label, FOOTER_TEXT)
            }
        }
    }
}

/// The scale bars: one for the exact clusters (they share one scale), one
/// for a cluster still estimated, one between clusters. Each about `target`
/// pixels wide, rounded to 1, 2 or 5 × 10ⁿ SNPs.
fn bars(l: &TreeLayout, g: &Geometry, target: f32) -> Vec<FooterItem> {
    let mut out = Vec::new();
    let mut add = |r: &layout::Region, what: &str| {
        let Some(bar) = layout::scale_bar(r, f64::from(target / g.px_per_unit.max(1e-6))) else {
            return;
        };
        let px = bar.layout_width as f32 * g.px_per_unit;
        if !px.is_finite() || px < 10.0 {
            return;
        }
        let exact = !(r.approximate || r.focus < fmt::FOCUS_EXACT);
        out.push(FooterItem::Bar {
            px,
            exact,
            label: format!("{}{what}", bar.label),
        });
    };
    let clusters: Vec<&layout::Region> = l.regions.iter().filter(|r| r.cluster.is_some()).collect();
    let several = clusters.len() > 1 || l.regions.first().is_some_and(|r| r.raw_max > 0.0);
    if let Some(r) = clusters
        .iter()
        .find(|r| !r.approximate && r.focus >= fmt::FOCUS_EXACT)
    {
        add(r, if several { " within clusters" } else { "" });
    }
    if let Some(r) = clusters
        .iter()
        .find(|r| r.approximate || r.focus < fmt::FOCUS_EXACT)
    {
        let what = format!(" within c{}", r.cluster.unwrap_or_default());
        add(r, &what);
    }
    if let Some(bb) = l.regions.first().filter(|r| r.raw_max > 0.0) {
        add(bb, " between clusters");
    }
    out
}

fn legend_items(m: &Marks, used: &Used) -> Vec<FooterItem> {
    let mut v = Vec::new();
    let mut push = |mark, label: String| v.push(FooterItem::Legend { mark, label });
    if used.exact {
        push(LegendMark::Branch(true), "exact SNPs".into());
    }
    if used.estimate {
        push(LegendMark::Branch(false), "≈ estimate".into());
    }
    if used.tip_snps {
        push(
            LegendMark::TipNumber,
            "SNPs on the genome's own branch".into(),
        );
    }
    if used.groups {
        push(
            LegendMark::Group,
            format!(
                "SNP group: linked by ≤{} SNP{}",
                m.group_threshold,
                if m.group_threshold == 1 { "" } else { "s" }
            ),
        );
    }
    if used.wedges {
        push(LegendMark::Wedge, "collapsed cluster".into());
    }
    if used.selection {
        push(
            LegendMark::Selected,
            format!("selected ({})", m.selection.len()),
        );
    }
    v
}

fn footer_flow(items: &[FooterItem], width: f32) -> (Vec<(f32, usize)>, usize) {
    let mut out = Vec::with_capacity(items.len());
    let (mut x, mut line) = (0.0f32, 0usize);
    for it in items {
        let w = it.width();
        if x > 0.0 && x + w > width {
            line += 1;
            x = 0.0;
        }
        out.push((x, line));
        x += w + 22.0;
    }
    (out, if items.is_empty() { 0 } else { line + 1 })
}

fn footer_items(
    l: &TreeLayout,
    m: &Marks,
    g: &Geometry,
    used: &Used,
    legend: bool,
) -> Vec<FooterItem> {
    let mut items = bars(l, g, 64.0);
    if legend {
        items.extend(legend_items(m, used));
    }
    items
}

/// Height of the footer for a picture `width` pixels wide.
pub fn footer_height(
    l: &TreeLayout,
    m: &Marks,
    g: &Geometry,
    used: &Used,
    legend: bool,
    width: f32,
) -> f32 {
    let items = footer_items(l, m, g, used, legend);
    let (_, lines) = footer_flow(&items, (width - 2.0 * FOOTER_PAD).max(50.0));
    if lines == 0 {
        0.0
    } else {
        lines as f32 * FOOTER_LINE + FOOTER_PAD
    }
}

/// The footer, in `area`: scale bars first, then the legend of what the
/// scene drew. `rule` draws a hairline along its top edge.
pub fn footer(
    l: &TreeLayout,
    m: &Marks,
    g: &Geometry,
    used: &Used,
    legend: bool,
    ink: &Ink,
    area: Rect,
    rule: bool,
) -> Vec<Prim> {
    let items = footer_items(l, m, g, used, legend);
    let mut out = Vec::new();
    if items.is_empty() {
        return out;
    }
    out.push(Prim::Rect {
        r: area,
        fill: ink.paper,
        rounding: 0.0,
    });
    if rule {
        out.push(Prim::Line {
            a: area.left_top(),
            b: area.right_top(),
            width: 1.0,
            color: ink.hairline,
        });
    }
    let (place, _) = footer_flow(&items, (area.width() - 2.0 * FOOTER_PAD).max(50.0));
    for (it, (x, line)) in items.iter().zip(place) {
        let x0 = area.left() + FOOTER_PAD + x;
        let y = area.top() + FOOTER_PAD / 2.0 + FOOTER_LINE * (line as f32 + 0.5);
        match it {
            FooterItem::Bar { px, exact, label } => {
                let c = ink.branch(*exact);
                out.push(Prim::Line {
                    a: Pos2::new(x0, y),
                    b: Pos2::new(x0 + px, y),
                    width: BRANCH_W,
                    color: c,
                });
                for xx in [x0, x0 + px] {
                    out.push(Prim::Line {
                        a: Pos2::new(xx, y - 3.5),
                        b: Pos2::new(xx, y + 3.5),
                        width: BRANCH_W,
                        color: c,
                    });
                }
                out.push(text(
                    Pos2::new(x0 + px + 8.0, y),
                    HAlign::Left,
                    label.clone(),
                    FOOTER_TEXT,
                    if *exact { ink.text } else { ink.text_weak },
                ));
            }
            FooterItem::Legend { mark, label } => {
                let tx = match mark {
                    LegendMark::Branch(exact) => {
                        out.push(Prim::Line {
                            a: Pos2::new(x0, y),
                            b: Pos2::new(x0 + 18.0, y),
                            width: BRANCH_W,
                            color: ink.branch(*exact),
                        });
                        x0 + 24.0
                    }
                    LegendMark::Group => {
                        out.push(Prim::Rect {
                            r: Rect::from_min_size(Pos2::new(x0, y - 6.0), Vec2::new(18.0, 12.0)),
                            fill: ink.group_fill(1),
                            rounding: 3.0,
                        });
                        x0 + 24.0
                    }
                    LegendMark::Selected => {
                        out.push(Prim::Rect {
                            r: Rect::from_min_size(Pos2::new(x0, y - 6.0), Vec2::new(18.0, 12.0)),
                            fill: ink.select_band,
                            rounding: 2.0,
                        });
                        x0 + 24.0
                    }
                    LegendMark::Wedge => {
                        let c = ink.text_weak;
                        out.push(Prim::Poly {
                            pts: vec![
                                Pos2::new(x0 + 1.0, y),
                                Pos2::new(x0 + 17.0, y - 6.0),
                                Pos2::new(x0 + 17.0, y + 6.0),
                            ],
                            fill: crate::theme::with_alpha(c, 70),
                            stroke: (1.0, c),
                        });
                        x0 + 24.0
                    }
                    LegendMark::TipNumber => {
                        out.push(text(
                            Pos2::new(x0, y),
                            HAlign::Left,
                            "+3",
                            FOOTER_TEXT,
                            ink.text_weak,
                        ));
                        x0 + text_width("+3", FOOTER_TEXT) + 6.0
                    }
                };
                out.push(text(
                    Pos2::new(tx, y),
                    HAlign::Left,
                    label.clone(),
                    FOOTER_TEXT,
                    ink.text,
                ));
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// The whole figure (SVG, tests)
// ---------------------------------------------------------------------------

const FIG_MARGIN: f32 = 18.0;

/// The whole tree at `row_px` per row and `tree_width` pixels for the tree
/// itself, with the footer under it: the picture an SVG export draws.
pub fn whole_scene(
    f: &Figure,
    o: &PubOptions,
    ink: &Ink,
    row_px: f32,
    tree_width: f32,
) -> (Vec<Prim>, Vec2) {
    let l = &f.layout;
    let px_per_unit = tree_width / l.x_max.max(1e-9) as f32;
    let origin = Pos2::new(FIG_MARGIN, FIG_MARGIN + row_px * 0.5);
    let tree_h = l.rows as f32 * row_px;
    let room = label_room(l, &f.marks, o, row_px);
    let width = (FIG_MARGIN * 2.0 + tree_width + room).ceil();
    let mut g = Geometry {
        origin,
        px_per_row: row_px,
        px_per_unit,
        clip: Rect::from_min_size(Pos2::ZERO, Vec2::new(width, FIG_MARGIN + tree_h + 4.0)),
    };
    let sc = scene(l, &f.marks, &g, o, ink, None);
    let fh = footer_height(l, &f.marks, &g, &sc.used, o.legend, width);
    let height = (FIG_MARGIN + tree_h + 10.0 + fh + FIG_MARGIN).ceil();
    g.clip = Rect::from_min_size(Pos2::ZERO, Vec2::new(width, height));
    let mut prims = sc.prims;
    let area = Rect::from_min_size(
        Pos2::new(FIG_MARGIN - FOOTER_PAD, FIG_MARGIN + tree_h + 10.0),
        Vec2::new(width - 2.0 * FIG_MARGIN + 2.0 * FOOTER_PAD, fh),
    );
    prims.extend(footer(
        l, &f.marks, &g, &sc.used, o.legend, ink, area, false,
    ));
    (prims, Vec2::new(width, height))
}

/// Cluster colours and states from the model: the part of [`Marks`] every
/// picture of the tree needs.
pub fn marks_from_model(model: &crate::state::UiModel) -> Marks {
    let mut m = Marks::default();
    for c in &model.clusters {
        m.clusters.insert(
            c.id,
            ClusterBrief {
                size: c.size(),
                focus: c.focus,
                approximate: c.approximate,
            },
        );
        for &g in &c.members {
            m.genome_cluster.insert(g, c.id);
        }
    }
    m
}

// ---------------------------------------------------------------------------
// Painting
// ---------------------------------------------------------------------------

/// Paint the primitives with egui. The SVG writer reads the same list.
/// Axis-aligned lines are snapped to pixel centres so the branches are crisp.
pub fn paint(painter: &egui::Painter, prims: &[Prim]) {
    for p in prims {
        match p {
            Prim::Line { a, b, width, color } => {
                let (mut a, mut b) = (*a, *b);
                if (a.y - b.y).abs() < 0.01 {
                    let y = painter.round_to_pixel_center(a.y);
                    a.y = y;
                    b.y = y;
                } else if (a.x - b.x).abs() < 0.01 {
                    let x = painter.round_to_pixel_center(a.x);
                    a.x = x;
                    b.x = x;
                }
                painter.line_segment([a, b], egui::Stroke::new(*width, *color));
            }
            Prim::Circle { c, r, fill, stroke } => {
                if fill.a() > 0 {
                    painter.circle_filled(*c, *r, *fill);
                }
                if let Some((w, col)) = stroke {
                    painter.circle_stroke(*c, *r + w * 0.5, egui::Stroke::new(*w, *col));
                }
            }
            Prim::Rect { r, fill, rounding } => {
                painter.rect_filled(*r, *rounding, *fill);
            }
            Prim::Poly { pts, fill, stroke } => {
                painter.add(egui::Shape::convex_polygon(
                    pts.clone(),
                    *fill,
                    egui::Stroke::new(stroke.0, stroke.1),
                ));
            }
            Prim::Text {
                pos,
                align,
                text,
                size,
                color,
                bold,
                tail,
            } => {
                let anchor = match align {
                    HAlign::Left => egui::Align2::LEFT_CENTER,
                    HAlign::Center => egui::Align2::CENTER_CENTER,
                };
                let font = egui::FontId::proportional(*size);
                let r = painter.text(*pos, anchor, text, font.clone(), *color);
                if *bold {
                    // No bold face is loaded: a second pass, half a pixel over.
                    painter.text(*pos + Vec2::new(0.55, 0.0), anchor, text, font, *color);
                }
                if let Some(t) = tail {
                    painter.text(
                        Pos2::new(r.right() + tail_gap(*size), pos.y),
                        egui::Align2::LEFT_CENTER,
                        &t.text,
                        egui::FontId::proportional(t.size),
                        t.color,
                    );
                }
            }
        }
    }
}

/// A small figure for tests here and in the tree window.
#[cfg(test)]
pub(crate) mod tests_support {
    use super::*;
    use crate::state::UiModel;
    use infinity_core::events::{ClusterFocusInfo, Event};

    /// Two clusters, the first exact and opened, the second collapsed; one
    /// genome of each selected; g1 and g2 in one SNP group.
    pub fn figure() -> Figure {
        let mut m = UiModel::default();
        for i in 0..5u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2]), (1, 3, &[3, 4])]),
        );
        m.apply(
            0.0,
            Event::BackboneTree {
                newick: "(g0:0.002,g3:0.004);".into(),
                pivots: vec!["g0".into(), "g3".into()],
            },
        );
        m.apply(
            0.0,
            Event::ClusterFocus(ClusterFocusInfo {
                cluster: 0,
                focus: 2,
                approximate: false,
                newick: "(g0:1,(g1:2,g2:4):1);".into(),
                members: vec![],
            }),
        );
        let open: HashSet<u32> = [0u32, 1].into_iter().collect();
        let layout = layout::layout(&m.display_tree(&open).unwrap());
        let mut marks = marks_from_model(&m);
        marks.selection = [1u32, 4].into_iter().collect();
        marks.selected_in = [(1u32, 1usize)].into_iter().collect();
        marks.groups = vec![GroupBand {
            number: 1,
            members: vec![1, 2],
        }];
        marks.group_threshold = 5;
        Figure { layout, marks }
    }
}

#[cfg(test)]
mod tests {
    use super::tests_support::figure;
    use super::*;

    fn texts(p: &[Prim]) -> Vec<String> {
        p.iter()
            .filter_map(|p| match p {
                Prim::Text { text, .. } => Some(text.clone()),
                _ => None,
            })
            .collect()
    }

    fn tails(p: &[Prim]) -> Vec<(String, String)> {
        p.iter()
            .filter_map(|p| match p {
                Prim::Text {
                    text,
                    tail: Some(t),
                    ..
                } => Some((text.clone(), t.text.clone())),
                _ => None,
            })
            .collect()
    }

    fn whole(o: &PubOptions) -> Vec<Prim> {
        whole_scene(&figure(), o, &Ink::PAPER, 16.0, 500.0).0
    }

    #[test]
    fn tips_carry_the_snps_of_their_own_branch() {
        let p = whole(&PubOptions::default());
        let t = tails(&p);
        // (g0:1,(g1:2,g2:4):1): g0 has 1 SNP of its own, g1 2, g2 4.
        assert!(t.contains(&("g0".into(), "+1".into())), "{t:?}");
        assert!(t.contains(&("g1".into(), "+2".into())));
        assert!(t.contains(&("g2".into(), "+4".into())));
        let off = whole(&PubOptions {
            tip_snps: false,
            ..Default::default()
        });
        assert!(tails(&off).is_empty());
        assert_eq!(tip_number(0.0, false), "+0");
        assert_eq!(tip_number(1200.0, true), "≈+1\u{2009}200");
    }

    #[test]
    fn every_branch_has_one_width_and_no_dashes_and_exact_ones_are_ink() {
        let p = whole(&PubOptions::default());
        let branches: Vec<(f32, Color32)> = p
            .iter()
            .filter_map(|p| match p {
                Prim::Line { width, color, .. }
                    if *color == Ink::PAPER.exact || *color == Ink::PAPER.estimate =>
                {
                    Some((*width, *color))
                }
                _ => None,
            })
            .collect();
        assert!(branches.iter().all(|(w, _)| *w == BRANCH_W));
        assert!(branches.iter().any(|(_, c)| *c == Ink::PAPER.exact));
        assert!(branches.iter().any(|(_, c)| *c == Ink::PAPER.estimate));
    }

    #[test]
    fn vertical_connectors_span_exactly_their_children() {
        let f = figure();
        let g = Geometry {
            origin: Pos2::new(10.0, 10.0),
            px_per_row: 20.0,
            px_per_unit: 300.0,
            clip: Rect::from_min_size(Pos2::ZERO, Vec2::new(2000.0, 2000.0)),
        };
        let s = scene(
            &f.layout,
            &f.marks,
            &g,
            &PubOptions::default(),
            &Ink::PAPER,
            None,
        );
        let l = &f.layout;
        for n in &l.nodes {
            let kids: Vec<f32> = l
                .nodes
                .iter()
                .filter(|c| c.parent == Some(n.idx))
                .map(|c| g.screen(c.x, c.y).y)
                .collect();
            if kids.len() < 2 {
                continue;
            }
            let x = g.screen(n.x, 0.0).x;
            let lo = kids.iter().copied().fold(f32::MAX, f32::min);
            let hi = kids.iter().copied().fold(f32::MIN, f32::max);
            let v = s
                .prims
                .iter()
                .find_map(|p| match p {
                    Prim::Line { a, b, .. } if (a.x - x).abs() < 1e-3 && (b.x - x).abs() < 1e-3 => {
                        Some((a.y, b.y))
                    }
                    _ => None,
                })
                .expect("a vertical connector");
            assert!((v.0 - (lo - BRANCH_W / 2.0)).abs() < 1e-3, "{v:?} vs {lo}");
            assert!((v.1 - (hi + BRANCH_W / 2.0)).abs() < 1e-3, "{v:?} vs {hi}");
        }
    }

    #[test]
    fn groups_are_shaded_and_labelled_and_the_legend_lists_only_what_is_drawn() {
        let p = whole(&PubOptions::default());
        let t = texts(&p);
        assert!(t.iter().any(|x| x == "G1 · 2 genomes ≤5 SNPs"), "{t:?}");
        assert!(p.iter().any(|p| matches!(p,
            Prim::Rect { fill, .. } if *fill == Ink::PAPER.group_fill(1))));
        assert!(t.iter().any(|x| x.starts_with("SNP group: linked by ≤5")));
        assert!(t.iter().any(|x| x == "exact SNPs"));
        assert!(t.iter().any(|x| x == "selected (2)"));
        // Cluster 1 has no tree yet, so it stays a wedge, and says so.
        assert!(t.iter().any(|x| x == "collapsed cluster"));
        assert!(
            !t.iter().any(|x| x.contains("Focus")),
            "no legend for what is not drawn"
        );
        let off = whole(&PubOptions {
            groups: false,
            ..Default::default()
        });
        assert!(!texts(&off).iter().any(|x| x.starts_with("G1")));
        assert!(!texts(&off).iter().any(|x| x.starts_with("SNP group")));
    }

    #[test]
    fn the_scale_bars_name_their_scope_and_estimates_say_so() {
        let t = texts(&whole(&PubOptions::default()));
        assert!(
            t.iter().any(|x| x.ends_with("SNPs within clusters")),
            "{t:?}"
        );
        let bb = t
            .iter()
            .find(|x| x.ends_with("between clusters"))
            .expect("a between-cluster bar");
        assert!(bb.starts_with('≈'), "{bb}");
    }

    #[test]
    fn brackets_and_branch_labels_follow_their_switches() {
        let with = whole(&PubOptions {
            branch_labels: true,
            ..Default::default()
        });
        let t = texts(&with);
        assert!(t.iter().any(|x| x == "c0 · 3"), "{t:?}");
        // The internal (g1,g2) branch is 1 SNP long; tips are not labelled
        // above their branch (their number is after the name).
        assert!(t.iter().any(|x| x == "1"), "{t:?}");
        assert!(!t.iter().any(|x| x == "4"));
        let without = whole(&PubOptions {
            cluster_brackets: false,
            legend: false,
            ..Default::default()
        });
        let t = texts(&without);
        assert!(!t.iter().any(|x| x.starts_with("c0 ·")));
        assert!(!t.iter().any(|x| x == "exact SNPs"));
    }

    #[test]
    fn a_collapsed_cluster_is_a_wedge_with_a_clean_label() {
        let mut f = figure();
        let mut m = crate::state::UiModel::default();
        for i in 0..5u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2]), (1, 3, &[3, 4])]),
        );
        m.apply(
            0.0,
            infinity_core::events::Event::BackboneTree {
                newick: "(g0:0.002,g3:0.004);".into(),
                pivots: vec!["g0".into(), "g3".into()],
            },
        );
        f.layout = layout::layout(&m.display_tree(&HashSet::new()).unwrap());
        let (p, _) = whole_scene(&f, &PubOptions::default(), &Ink::PAPER, 16.0, 400.0);
        let t = texts(&p);
        // c0 is exact in the marks, c1 not refined yet.
        assert!(t.iter().any(|x| x == "c0 · 3 genomes · exact"), "{t:?}");
        assert!(t.iter().any(|x| x == "c1 · 2 genomes · ≈ estimate"));
        assert!(t.iter().any(|x| x == "1 selected"));
        assert!(p.iter().any(|p| matches!(p, Prim::Poly { .. })));
        assert!(t.iter().any(|x| x == "collapsed cluster"));
    }

    #[test]
    fn clicks_find_genomes_clusters_and_groups() {
        let f = figure();
        let g = Geometry {
            origin: Pos2::new(10.0, 10.0),
            px_per_row: 16.0,
            px_per_unit: 300.0,
            clip: Rect::from_min_size(Pos2::ZERO, Vec2::new(2000.0, 2000.0)),
        };
        let s = scene(
            &f.layout,
            &f.marks,
            &g,
            &PubOptions::default(),
            &Ink::PAPER,
            None,
        );
        assert!(s.hits.iter().any(|(h, _)| *h == Hit::Genome(1)));
        assert!(s.hits.iter().any(|(h, _)| *h == Hit::Cluster(0)));
        assert!(s.hits.iter().any(|(h, _)| *h == Hit::Group(0)));
    }

    #[test]
    fn the_mrca_of_two_sisters_is_their_parent() {
        let f = figure();
        let l = &f.layout;
        let find = |name: &str| {
            l.nodes
                .iter()
                .position(|n| n.kind.name() == Some(name))
                .unwrap()
        };
        let (a, b, c) = (find("g1"), find("g2"), find("g0"));
        assert_eq!(mrca(l, &[a, b]), l.nodes[a].parent.unwrap());
        assert_eq!(mrca(l, &[a]), a);
        let top = mrca(l, &[a, c]);
        assert!(matches!(l.nodes[top].kind, NodeKind::ClusterRoot { .. }));
    }

    #[test]
    fn the_room_for_labels_grows_with_names_and_decorations() {
        let f = figure();
        let o = PubOptions::default();
        let bare = PubOptions {
            tip_snps: false,
            groups: false,
            cluster_brackets: false,
            ..o
        };
        assert!(
            label_room(&f.layout, &f.marks, &o, 16.0)
                > label_room(&f.layout, &f.marks, &bare, 16.0)
        );
        assert!(
            label_room(&f.layout, &f.marks, &bare, 3.0) < 20.0,
            "no labels at 3 px rows"
        );
    }
}
