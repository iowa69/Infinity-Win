//! The tree in publication style: white background, near-black branches,
//! leaves coloured by cluster, cluster brackets, SNP scale bars, a legend and
//! the selection highlight, whatever the theme of the window.
//!
//! The picture is built once as a list of drawing primitives in pixels
//! ([`Prim`]), and the same list is painted in the tree window and written as
//! SVG by [`crate::svg::from_scene`], so the saved figure is exactly what the
//! window shows. Pure data, no `egui::Ui`: testable without a window.
//!
//! Honesty rules carry over from the main tree (section 0 rule 8): exact SNP
//! branches are solid black, estimated ones grey (dashed for sketch-only), and
//! every estimated number, including the between-cluster scale bar, carries
//! "≈" and the word "estimate".

use std::collections::{BTreeMap, HashMap, HashSet};

use egui::{Color32, Pos2, Rect, Vec2};

use crate::colors::{self, FocusStyle};
use crate::fmt;
use crate::layout::{self, NodeKind, TreeLayout};

/// Paper.
pub const PAPER: Color32 = Color32::WHITE;
/// Exact branches and main text.
pub const INK: Color32 = Color32::from_rgb(22, 22, 22);
/// Focus 1 branches and secondary text.
pub const INK_MID: Color32 = Color32::from_rgb(96, 96, 96);
/// Focus 0 branches and hairlines.
pub const INK_SOFT: Color32 = Color32::from_rgb(160, 160, 160);
/// The selection highlight: a highlighter-pen band behind the row.
pub const HIGHLIGHT: Color32 = Color32::from_rgb(255, 214, 90);
/// The selection marks drawn on top: ring, tick, badge.
pub const HIGHLIGHT_INK: Color32 = Color32::from_rgb(196, 98, 0);

/// Branch style on paper: black and solid when exact, grey when estimated,
/// dashed when only the sketch is behind it.
pub fn branch_style(focus: u8) -> FocusStyle {
    match focus {
        fmt::FOCUS_EXACT => FocusStyle {
            color: INK,
            width: 1.4,
            dashed: false,
            dash: 0.0,
        },
        fmt::FOCUS_BLOCKS => FocusStyle {
            color: INK_MID,
            width: 1.1,
            dashed: false,
            dash: 0.0,
        },
        _ => FocusStyle {
            color: INK_SOFT,
            width: 1.0,
            dashed: true,
            dash: 3.0,
        },
    }
}

/// A cluster's colour on paper: the light-theme colour, which is dark enough
/// to read as text on white.
pub fn cluster_ink(id: u32) -> Color32 {
    colors::cluster_color(id, false)
}

/// One thing to draw, in pixels.
#[derive(Debug, Clone, PartialEq)]
pub enum Prim {
    Line {
        a: Pos2,
        b: Pos2,
        width: f32,
        color: Color32,
        /// Dash length, for a dashed line.
        dash: Option<f32>,
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
    /// Text whose left edge (or centre, or right edge) is at `pos.x`, centred
    /// vertically on `pos.y`.
    Text {
        pos: Pos2,
        align: HAlign,
        text: String,
        size: f32,
        color: Color32,
        bold: bool,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HAlign {
    Left,
    Center,
}

/// What a click on the picture hit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hit {
    Genome(u32),
    Cluster(u32),
}

/// One cluster, as far as the picture needs it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterBrief {
    pub size: usize,
    pub focus: u8,
    pub approximate: bool,
}

/// What to draw. Everything is owned, so the tree window can keep it between
/// frames without holding on to the model.
#[derive(Debug, Clone, Default)]
pub struct Figure {
    pub layout: TreeLayout,
    pub genome_cluster: HashMap<u32, u32>,
    pub clusters: BTreeMap<u32, ClusterBrief>,
    pub selection: HashSet<u32>,
    /// Selected genomes per cluster, for the badge on a collapsed cluster.
    pub selected_in: HashMap<u32, usize>,
}

/// Decorations the user can switch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PubOptions {
    pub branch_labels: bool,
    pub legend: bool,
    pub cluster_brackets: bool,
}

impl Default for PubOptions {
    fn default() -> Self {
        Self {
            branch_labels: false,
            legend: true,
            cluster_brackets: true,
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

/// Estimated width of a text, the same in the window and in the SVG so the
/// two lay out alike.
pub fn text_width(text: &str, size: f32) -> f32 {
    text.chars().count() as f32 * size * 0.56
}

/// Label font size for a row height.
pub fn label_size(px_per_row: f32) -> f32 {
    (px_per_row * 0.72).clamp(7.5, 14.0)
}

/// Below this row height the tip labels are dropped (they would overlap).
pub const MIN_ROW_FOR_LABELS: f32 = 6.0;

fn genome_label(name: &str, approximate: bool) -> String {
    if approximate {
        format!("≈ {name}")
    } else {
        name.to_string()
    }
}

fn collapsed_label(id: u32, c: Option<&ClusterBrief>, members: usize) -> String {
    format!(
        "c{id} ({members} genomes){}",
        if c.is_some_and(|c| c.approximate) {
            " ≈"
        } else {
            ""
        }
    )
}

/// The right end of the widest tip label, in pixels: where the cluster
/// brackets line up. With the labels hidden (rows too thin), the right end of
/// the tips themselves.
pub fn labels_right(f: &Figure, g: &Geometry) -> f32 {
    let size = label_size(g.px_per_row);
    let labels = g.px_per_row >= MIN_ROW_FOR_LABELS;
    let mut right = f32::MIN;
    for &i in &f.layout.leaves {
        let n = &f.layout.nodes[i];
        let s = g.screen(n.x, n.y);
        if !labels {
            let tip = match &n.kind {
                NodeKind::ClusterCollapsed { .. } => {
                    layout::COLLAPSED_SPREAD as f32 * g.px_per_unit
                }
                _ => 0.0,
            };
            right = right.max(s.x + 4.0 + tip);
            continue;
        }
        let w = match &n.kind {
            NodeKind::Genome { name, .. } => text_width(&genome_label(name, n.approximate), size),
            NodeKind::ClusterCollapsed { id, members, .. } => {
                text_width(&collapsed_label(*id, f.clusters.get(id), *members), size)
                    + (layout::COLLAPSED_SPREAD as f32 * g.px_per_unit)
                    + if f.selected_in.contains_key(id) {
                        70.0
                    } else {
                        0.0
                    }
            }
            _ => 0.0,
        };
        right = right.max(s.x + 7.0 + w);
    }
    if right == f32::MIN {
        g.origin.x
    } else {
        right
    }
}

/// The tree itself: highlight bands, branches, tips, labels, branch lengths
/// and cluster brackets. Returns the primitives and the click targets.
pub fn scene(f: &Figure, g: &Geometry, o: &PubOptions) -> (Vec<Prim>, Vec<(Hit, Rect)>) {
    let l = &f.layout;
    let mut out: Vec<Prim> = Vec::with_capacity(l.nodes.len() * 4);
    let mut hits: Vec<(Hit, Rect)> = Vec::new();
    let row = g.px_per_row;
    let show_labels = row >= MIN_ROW_FOR_LABELS;
    let size = label_size(row);

    // --- selection bands, behind everything -----------------------------------
    for &i in &l.leaves {
        let n = &l.nodes[i];
        let y = g.screen(n.x, n.y).y;
        if !g.row_visible(y, row) {
            continue;
        }
        let strength = match &n.kind {
            NodeKind::Genome { id, .. } if f.selection.contains(id) => Some(150u8),
            NodeKind::ClusterCollapsed { id, .. } if f.selected_in.contains_key(id) => Some(80),
            _ => None,
        };
        if let Some(a) = strength {
            let h = row.max(2.0);
            out.push(Prim::Rect {
                r: Rect::from_min_size(
                    Pos2::new(g.clip.left(), y - h / 2.0),
                    Vec2::new(g.clip.width(), h),
                ),
                fill: colors::wash(HIGHLIGHT, a),
                rounding: 0.0,
            });
        }
    }

    // --- branches -------------------------------------------------------------
    for n in &l.nodes {
        let Some(p) = n.parent.map(|p| &l.nodes[p]) else {
            continue;
        };
        let child = g.screen(n.x, n.y);
        let parent = g.screen(p.x, n.y);
        if !g.row_visible(child.y, 2.0) {
            continue;
        }
        push_line(&mut out, parent, child, &branch_style(n.edge_focus));
    }
    for n in l.nodes.iter().filter(|n| n.leaf_count >= 2) {
        let top = g.screen(n.x, n.row_min);
        let bottom = g.screen(n.x, n.row_max);
        if bottom.y < g.clip.top() - 2.0 || top.y > g.clip.bottom() + 2.0 {
            continue;
        }
        push_line(&mut out, top, bottom, &branch_style(n.focus));
    }

    // --- branch lengths ---------------------------------------------------------
    if o.branch_labels && row >= 10.0 {
        let fs = (row * 0.55).clamp(8.0, 11.0);
        for n in &l.nodes {
            let Some(p) = n.parent.map(|p| &l.nodes[p]) else {
                continue;
            };
            if n.len <= 0.0 || matches!(n.kind, NodeKind::Root) {
                continue;
            }
            let a = g.screen(p.x, n.y);
            let b = g.screen(n.x, n.y);
            if !g.row_visible(b.y, 2.0) {
                continue;
            }
            let text = format!(
                "{}{}",
                if n.edge_approximate { "≈" } else { "" },
                fmt::snps(n.len)
            );
            if b.x - a.x < text_width(&text, fs) + 6.0 {
                continue;
            }
            out.push(Prim::Text {
                pos: Pos2::new((a.x + b.x) / 2.0, b.y - fs * 0.6 - 1.0),
                align: HAlign::Center,
                text,
                size: fs,
                color: if n.edge_approximate {
                    INK_SOFT
                } else {
                    INK_MID
                },
                bold: false,
            });
        }
    }

    // --- tips -----------------------------------------------------------------
    for n in &l.nodes {
        let s = g.screen(n.x, n.y);
        if !g.row_visible(s.y, row + 4.0) {
            continue;
        }
        match &n.kind {
            NodeKind::Genome { id, name } => {
                let selected = f.selection.contains(id);
                let colour = f
                    .genome_cluster
                    .get(id)
                    .map(|&c| cluster_ink(c))
                    .unwrap_or(INK_SOFT);
                out.push(Prim::Circle {
                    c: s,
                    r: if selected { 3.0 } else { 2.2 },
                    fill: colour,
                    stroke: selected.then_some((1.6, HIGHLIGHT_INK)),
                });
                let mut hit = Rect::from_center_size(s, Vec2::splat(8.0));
                if show_labels {
                    let text = genome_label(name, n.approximate);
                    let x = s.x + 7.0;
                    let w = text_width(&text, size);
                    out.push(Prim::Text {
                        pos: Pos2::new(x, s.y),
                        align: HAlign::Left,
                        text,
                        size,
                        color: colour,
                        bold: selected,
                    });
                    hit = hit.union(Rect::from_min_size(
                        Pos2::new(x, s.y - row / 2.0),
                        Vec2::new(w, row),
                    ));
                }
                hits.push((Hit::Genome(*id), hit));
            }
            NodeKind::ClusterCollapsed { id, members, .. } => {
                let half = layout::collapsed_half_rows(*members);
                let t = layout::triangle(n.x, n.y, half);
                let pts: Vec<Pos2> = t.iter().map(|&(x, y)| g.screen(x, y)).collect();
                let colour = cluster_ink(*id);
                let base_x = pts[1].x;
                out.push(Prim::Poly {
                    pts: pts.clone(),
                    fill: colors::wash(colour, 90),
                    stroke: (1.2, colour),
                });
                let mut hit = Rect::from_points(&pts);
                if show_labels {
                    let text = collapsed_label(*id, f.clusters.get(id), *members);
                    let w = text_width(&text, size);
                    out.push(Prim::Text {
                        pos: Pos2::new(base_x + 5.0, s.y),
                        align: HAlign::Left,
                        text,
                        size,
                        color: colour,
                        bold: false,
                    });
                    hit = hit.union(Rect::from_min_size(
                        Pos2::new(base_x + 5.0, s.y - row / 2.0),
                        Vec2::new(w, row),
                    ));
                    if let Some(&k) = f.selected_in.get(id) {
                        let bx = base_x + 12.0 + w;
                        let text = format!("{k} selected");
                        let bw = text_width(&text, size * 0.85) + 10.0;
                        let bh = size * 1.3;
                        out.push(Prim::Rect {
                            r: Rect::from_min_size(
                                Pos2::new(bx, s.y - bh / 2.0),
                                Vec2::new(bw, bh),
                            ),
                            fill: HIGHLIGHT_INK,
                            rounding: bh / 2.0,
                        });
                        out.push(Prim::Text {
                            pos: Pos2::new(bx + 5.0, s.y),
                            align: HAlign::Left,
                            text,
                            size: size * 0.85,
                            color: PAPER,
                            bold: true,
                        });
                    }
                }
                hits.push((Hit::Cluster(*id), hit));
            }
            NodeKind::ClusterRoot { id, .. } => {
                out.push(Prim::Circle {
                    c: s,
                    r: 2.6,
                    fill: cluster_ink(*id),
                    stroke: None,
                });
            }
            NodeKind::Internal | NodeKind::Root => {}
        }
    }

    // --- cluster brackets -------------------------------------------------------
    if o.cluster_brackets {
        let x = labels_right(f, g) + 8.0;
        for n in &l.nodes {
            let NodeKind::ClusterRoot { id, members, .. } = &n.kind else {
                continue;
            };
            let top = g.screen(n.x, n.row_min - 0.38).y;
            let bottom = g.screen(n.x, n.row_max + 0.38).y;
            if bottom < g.clip.top() || top > g.clip.bottom() {
                continue;
            }
            let colour = cluster_ink(*id);
            out.push(Prim::Line {
                a: Pos2::new(x, top),
                b: Pos2::new(x, bottom),
                width: 2.0,
                color: colour,
                dash: None,
            });
            for yy in [top, bottom] {
                out.push(Prim::Line {
                    a: Pos2::new(x - 4.0, yy),
                    b: Pos2::new(x, yy),
                    width: 2.0,
                    color: colour,
                    dash: None,
                });
            }
            let approx = f.clusters.get(id).is_some_and(|c| c.approximate);
            let text = format!("c{id} · {members}{}", if approx { " ≈" } else { "" });
            let ty = ((top + bottom) / 2.0).clamp(g.clip.top() + 8.0, g.clip.bottom() - 8.0);
            out.push(Prim::Text {
                pos: Pos2::new(x + 6.0, ty),
                align: HAlign::Left,
                text,
                size: (size * 1.05).clamp(9.0, 15.0),
                color: colour,
                bold: true,
            });
            hits.push((
                Hit::Cluster(*id),
                Rect::from_min_max(Pos2::new(x - 4.0, top), Pos2::new(x + 60.0, bottom)),
            ));
        }
    }
    (out, hits)
}

fn push_line(out: &mut Vec<Prim>, a: Pos2, b: Pos2, style: &FocusStyle) {
    if (b - a).length() < 0.05 {
        return;
    }
    out.push(Prim::Line {
        a,
        b,
        width: style.width,
        color: style.color,
        dash: style.dashed.then_some(style.dash.max(1.0)),
    });
}

/// The SNP scale bars, stacked upwards from `bottom_left`: one for the exact
/// clusters (they share one scale), one for a cluster still estimated, and one
/// for the distances between clusters, marked as an estimate unless the
/// backbone is exact.
pub fn scale_bars(f: &Figure, g: &Geometry, bottom_left: Pos2) -> Vec<Prim> {
    let l = &f.layout;
    let mut bars: Vec<(&layout::Region, String)> = Vec::new();
    if let Some(r) = l
        .regions
        .iter()
        .find(|r| r.cluster.is_some() && !r.approximate && r.focus >= fmt::FOCUS_EXACT)
    {
        bars.push((r, "within clusters, exact".to_string()));
    }
    if let Some(r) = l
        .regions
        .iter()
        .find(|r| r.cluster.is_some() && (r.approximate || r.focus < fmt::FOCUS_EXACT))
    {
        bars.push((
            r,
            format!("within c{}, estimate", r.cluster.unwrap_or_default()),
        ));
    }
    if let Some(bb) = l.regions.first().filter(|r| r.raw_max > 0.0) {
        let what = if bb.approximate || bb.focus < fmt::FOCUS_EXACT {
            "between clusters, estimate"
        } else {
            "between clusters"
        };
        bars.push((bb, what.to_string()));
    }
    let mut out = Vec::new();
    let mut y = bottom_left.y;
    for (region, what) in bars {
        let Some(bar) = layout::scale_bar(region, 0.18) else {
            continue;
        };
        let w = bar.layout_width as f32 * g.px_per_unit;
        if !w.is_finite() || w <= 4.0 {
            continue;
        }
        let x0 = bottom_left.x;
        let estimate = region.approximate || region.focus < fmt::FOCUS_EXACT;
        let style = FocusStyle {
            dashed: false,
            ..branch_style(if estimate {
                fmt::FOCUS_BLOCKS
            } else {
                fmt::FOCUS_EXACT
            })
        };
        let ink = if estimate { INK_MID } else { INK };
        push_line(&mut out, Pos2::new(x0, y), Pos2::new(x0 + w, y), &style);
        for x in [x0, x0 + w] {
            out.push(Prim::Line {
                a: Pos2::new(x, y - 4.0),
                b: Pos2::new(x, y + 4.0),
                width: 1.4,
                color: ink,
                dash: None,
            });
        }
        out.push(Prim::Text {
            pos: Pos2::new(x0 + w + 7.0, y),
            align: HAlign::Left,
            text: format!("{}  ({what})", bar.label),
            size: 11.0,
            color: ink,
            bold: false,
        });
        y -= 18.0;
    }
    out
}

/// Height the scale bars need, for reserving room under the tree.
pub fn scale_bars_height(f: &Figure) -> f32 {
    let n = f
        .layout
        .regions
        .iter()
        .filter(|r| r.raw_max > 0.0)
        .count()
        .min(3);
    18.0 * n as f32 + 10.0
}

/// Most clusters a legend lists one by one.
const LEGEND_CLUSTERS: usize = 18;

/// The legend, with its top-left corner at `at`: branch styles, cluster
/// colours, and the selection when there is one. Returns the primitives and
/// the box size.
pub fn legend(f: &Figure, at: Pos2) -> (Vec<Prim>, Vec2) {
    let mut items: Vec<(LegendMark, String)> = vec![
        (
            LegendMark::Branch(fmt::FOCUS_EXACT),
            "exact SNPs (Focus 2)".into(),
        ),
        (
            LegendMark::Branch(fmt::FOCUS_BLOCKS),
            "≈ estimate, blocks compared (Focus 1)".into(),
        ),
        (
            LegendMark::Branch(fmt::FOCUS_SKETCH),
            "≈ estimate, sketch only (Focus 0)".into(),
        ),
    ];
    if !f.selection.is_empty() {
        items.push((
            LegendMark::Selected,
            format!("selected ({})", f.selection.len()),
        ));
    }
    let mut cluster_items: Vec<(LegendMark, String)> = f
        .clusters
        .iter()
        .take(LEGEND_CLUSTERS)
        .map(|(&id, c)| {
            (
                LegendMark::Cluster(id),
                format!(
                    "c{id} · {}{}",
                    c.size,
                    if c.approximate { " ≈" } else { "" }
                ),
            )
        })
        .collect();
    if f.clusters.len() > LEGEND_CLUSTERS {
        cluster_items.push((
            LegendMark::None,
            format!("+ {} more clusters", f.clusters.len() - LEGEND_CLUSTERS),
        ));
    }
    let line_h = 15.0;
    let pad = 8.0;
    let size = 11.0;
    let col1_w = items
        .iter()
        .map(|(_, t)| text_width(t, size))
        .fold(0.0, f32::max)
        + 40.0;
    // Clusters in up to three columns beside the branch styles.
    let per_col = items.len().max(cluster_items.len().div_ceil(3)).max(1);
    let cols = cluster_items.len().div_ceil(per_col);
    let col_w = cluster_items
        .iter()
        .map(|(_, t)| text_width(t, size))
        .fold(0.0, f32::max)
        + 24.0;
    let rows = items.len().max(per_col.min(cluster_items.len()));
    let size_box = Vec2::new(
        pad * 2.0 + col1_w + cols as f32 * col_w,
        pad * 2.0 + rows as f32 * line_h,
    );
    let mut out = vec![Prim::Rect {
        r: Rect::from_min_size(at, size_box),
        fill: Color32::from_rgba_unmultiplied(255, 255, 255, 235),
        rounding: 5.0,
    }];
    // A hairline frame.
    let r = Rect::from_min_size(at, size_box);
    for (a, b) in [
        (r.left_top(), r.right_top()),
        (r.right_top(), r.right_bottom()),
        (r.right_bottom(), r.left_bottom()),
        (r.left_bottom(), r.left_top()),
    ] {
        out.push(Prim::Line {
            a,
            b,
            width: 0.8,
            color: Color32::from_rgb(210, 210, 210),
            dash: None,
        });
    }
    let draw = |out: &mut Vec<Prim>, x: f32, i: usize, mark: &LegendMark, text: &str| {
        let y = at.y + pad + line_h * (i as f32 + 0.5);
        match mark {
            LegendMark::Branch(focus) => push_line(
                out,
                Pos2::new(x, y),
                Pos2::new(x + 26.0, y),
                &branch_style(*focus),
            ),
            LegendMark::Cluster(id) => out.push(Prim::Circle {
                c: Pos2::new(x + 5.0, y),
                r: 4.0,
                fill: cluster_ink(*id),
                stroke: None,
            }),
            LegendMark::Selected => out.push(Prim::Rect {
                r: Rect::from_center_size(Pos2::new(x + 13.0, y), Vec2::new(26.0, 10.0)),
                fill: colors::wash(HIGHLIGHT, 150),
                rounding: 2.0,
            }),
            LegendMark::None => {}
        }
        let tx = match mark {
            LegendMark::Cluster(_) | LegendMark::None => x + 14.0,
            _ => x + 32.0,
        };
        out.push(Prim::Text {
            pos: Pos2::new(tx, y),
            align: HAlign::Left,
            text: text.to_string(),
            size,
            color: match mark {
                LegendMark::Cluster(id) => cluster_ink(*id),
                LegendMark::Branch(f) => branch_style(*f).color,
                _ => INK,
            },
            bold: false,
        });
    };
    for (i, (m, t)) in items.iter().enumerate() {
        draw(&mut out, at.x + pad, i, m, t);
    }
    for (k, (m, t)) in cluster_items.iter().enumerate() {
        let col = k / per_col;
        let i = k % per_col;
        draw(&mut out, at.x + pad + col1_w + col as f32 * col_w, i, m, t);
    }
    (out, size_box)
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum LegendMark {
    Branch(u8),
    Cluster(u32),
    Selected,
    None,
}

/// A geometry that shows the whole figure at `row_px` per row and
/// `tree_width` pixels for the tree, with room above for the legend and below
/// for the scale bars: the picture an SVG export draws. Returns the geometry
/// and the full picture size.
pub fn whole_figure(f: &Figure, o: &PubOptions, row_px: f32, tree_width: f32) -> (Geometry, Vec2) {
    let margin = 18.0;
    let legend_h = if o.legend {
        legend(f, Pos2::ZERO).1.y + 12.0
    } else {
        0.0
    };
    let px_per_unit = tree_width / f.layout.x_max.max(1e-9) as f32;
    let origin = Pos2::new(margin, margin + legend_h + row_px * 0.5);
    let height = origin.y + f.layout.rows as f32 * row_px + scale_bars_height(f) + margin;
    let mut g = Geometry {
        origin,
        px_per_row: row_px,
        px_per_unit,
        clip: Rect::from_min_size(Pos2::ZERO, Vec2::new(10_000.0, height)),
    };
    let mut width = labels_right(f, &g) + 24.0;
    if o.cluster_brackets {
        width += 90.0;
    }
    if o.legend {
        width = width.max(margin * 2.0 + legend(f, Pos2::ZERO).1.x);
    }
    g.clip = Rect::from_min_size(Pos2::ZERO, Vec2::new(width.ceil(), height.ceil()));
    (g, Vec2::new(width.ceil(), height.ceil()))
}

/// Everything of the whole figure, ready for [`crate::svg::from_scene`].
pub fn whole_scene(f: &Figure, o: &PubOptions, row_px: f32, tree_width: f32) -> (Vec<Prim>, Vec2) {
    let (g, size) = whole_figure(f, o, row_px, tree_width);
    let (mut prims, _) = scene(f, &g, o);
    // The bars stack upwards from the bottom margin, inside the room
    // `whole_figure` reserved under the last row.
    prims.extend(scale_bars(f, &g, Pos2::new(18.0, size.y - 22.0)));
    if o.legend {
        prims.extend(legend(f, Pos2::new(18.0, 18.0)).0);
    }
    (prims, size)
}

/// A small figure for tests here and in the tree window.
#[cfg(test)]
pub(crate) mod tests_support {
    use super::*;
    use crate::state::UiModel;
    use infinity_core::events::{ClusterFocusInfo, Event};

    /// Two clusters, the first exact and opened, the second collapsed; one
    /// genome of each selected.
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
        let mut f = Figure {
            layout,
            ..Default::default()
        };
        for c in &m.clusters {
            f.clusters.insert(
                c.id,
                ClusterBrief {
                    size: c.size(),
                    focus: c.focus,
                    approximate: c.approximate,
                },
            );
            for &g in &c.members {
                f.genome_cluster.insert(g, c.id);
            }
        }
        f.selection = [1u32, 4].into_iter().collect();
        f.selected_in = [(1u32, 1usize)].into_iter().collect();
        f
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

    #[test]
    fn the_figure_is_on_paper_with_black_exact_branches_and_cluster_coloured_tips() {
        let f = figure();
        let (prims, size) = whole_scene(&f, &PubOptions::default(), 16.0, 500.0);
        assert!(size.x > 500.0 && size.y > 5.0 * 16.0);
        let t = texts(&prims);
        for name in ["g0", "g1", "g2"] {
            assert!(t.iter().any(|x| x == name), "{name} missing from {t:?}");
        }
        // Exact branches are near-black and solid.
        assert!(prims.iter().any(|p| matches!(p,
            Prim::Line { color, dash: None, .. } if *color == INK)));
        // Tip labels take their cluster's colour.
        assert!(prims.iter().any(|p| matches!(p,
            Prim::Text { text, color, .. } if text == "g1" && *color == cluster_ink(0))));
        // Brackets name the clusters; the legend explains the styles.
        assert!(t.iter().any(|x| x.starts_with("c0 · 3")), "{t:?}");
        assert!(t.iter().any(|x| x == "exact SNPs (Focus 2)"));
        assert!(t.iter().any(|x| x == "selected (2)"));
    }

    #[test]
    fn selected_tips_are_bold_on_a_highlight_band() {
        let f = figure();
        let (prims, _) = whole_scene(&f, &PubOptions::default(), 16.0, 500.0);
        assert!(prims.iter().any(|p| matches!(p,
            Prim::Text { text, bold: true, .. } if text == "g1")));
        assert!(prims.iter().any(|p| matches!(p,
            Prim::Text { text, bold: false, .. } if text == "g2")));
        let bands = prims
            .iter()
            .filter(|p| {
                matches!(p, Prim::Rect { fill, .. }
                if fill.to_srgba_unmultiplied()[..3] == HIGHLIGHT.to_srgba_unmultiplied()[..3])
            })
            .count();
        assert!(bands >= 2, "one band per selected genome, got {bands}");
    }

    #[test]
    fn the_between_cluster_bar_says_it_is_an_estimate() {
        let f = figure();
        let (prims, _) = whole_scene(&f, &PubOptions::default(), 16.0, 500.0);
        let t = texts(&prims);
        let bar = t
            .iter()
            .find(|x| x.contains("between clusters"))
            .expect("a between-cluster bar");
        assert!(bar.contains("estimate") && bar.contains('≈'), "{bar}");
        assert!(
            t.iter().any(|x| x.contains("within clusters, exact")),
            "{t:?}"
        );
    }

    #[test]
    fn switching_decorations_off_removes_them() {
        let f = figure();
        let o = PubOptions {
            legend: false,
            cluster_brackets: false,
            branch_labels: false,
        };
        let (prims, _) = whole_scene(&f, &o, 16.0, 500.0);
        let t = texts(&prims);
        assert!(!t.iter().any(|x| x.contains("Focus 2")));
        assert!(!t.iter().any(|x| x.starts_with("c0 ·")));
        let with = whole_scene(
            &f,
            &PubOptions {
                branch_labels: true,
                ..o
            },
            20.0,
            500.0,
        );
        assert!(texts(&with.0).iter().any(|x| x == "4"), "the 4-SNP branch");
    }

    #[test]
    fn clicks_find_genomes() {
        let f = figure();
        let (g, _) = whole_figure(&f, &PubOptions::default(), 16.0, 500.0);
        let (_, hits) = scene(&f, &g, &PubOptions::default());
        assert!(hits.iter().any(|(h, _)| *h == Hit::Genome(1)));
        assert!(hits.iter().any(|(h, _)| *h == Hit::Cluster(0)));
    }
}
