//! The tree view as a vector picture (SVG), drawn with the same vocabulary as
//! the panel: focus styles (dashed and pale at Focus 0, thin orange at Focus 1,
//! solid and strong at Focus 2), cluster colours, collapsed-cluster triangles,
//! SNP scale bars, the legend, and the optional branch labels, aligned tips and
//! colour strip.
//!
//! Where the PNG export is a screenshot of what is on screen, this draws the
//! *whole* laid-out tree at a readable row height, so a tree of three hundred
//! genomes comes out as a figure rather than as a crop. Pure string building
//! over [`TreeLayout`]: testable without a window.

use std::fmt::Write as _;

use egui::Color32;

use crate::layout::{self, NodeKind, TreeLayout};
use crate::state::UiModel;
use crate::{colors, fmt, theme};

/// What the picture should look like.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SvgOptions {
    pub dark: bool,
    pub legend: bool,
    pub branch_labels: bool,
    pub aligned_tips: bool,
    pub color_strip: bool,
    /// Width of the tree drawing itself, without the tip labels, in px.
    pub tree_width: f32,
    /// Height of one leaf row, in px.
    pub row_px: f32,
    /// The cluster whose scale bar is drawn, when several are open.
    pub selected_cluster: Option<u32>,
}

impl Default for SvgOptions {
    fn default() -> Self {
        Self {
            dark: false,
            legend: true,
            branch_labels: false,
            aligned_tips: false,
            color_strip: true,
            tree_width: 760.0,
            row_px: 16.0,
            selected_cluster: None,
        }
    }
}

const MARGIN: f32 = 16.0;
const LEGEND_H: f32 = 70.0;
const FONT: &str = "Helvetica, Arial, sans-serif";

/// Render a laid-out tree as a complete SVG document.
pub fn render(l: &TreeLayout, model: &UiModel, o: &SvgOptions) -> String {
    let p = theme::palette_for(o.dark);
    let row = o.row_px.clamp(4.0, 60.0);
    let label_px = (row * 0.72).clamp(8.5, 14.0);
    let top = MARGIN + if o.legend { LEGEND_H } else { 0.0 };
    let px_per_unit = o.tree_width / (l.x_max.max(1e-9) as f32);
    let sx = |x: f64| MARGIN + x as f32 * px_per_unit;
    let sy = |y: f64| top + (y as f32 + 0.5) * row;

    // The widest tip label decides the picture's width.
    let longest = l
        .nodes
        .iter()
        .map(|n| match &n.kind {
            NodeKind::Genome { name, .. } => name.chars().count() + 2,
            NodeKind::ClusterCollapsed { .. } => 28,
            _ => 0,
        })
        .max()
        .unwrap_or(0) as f32;
    let label_w = longest * label_px * 0.6 + if o.color_strip { 14.0 } else { 0.0 } + 16.0;
    let width = (MARGIN * 2.0 + o.tree_width + label_w).ceil();
    let height = (top + l.rows as f32 * row + 56.0).ceil();
    let tip_x = sx(l.x_max) + 10.0;

    let mut s = String::with_capacity(64 * l.nodes.len() + 2048);
    let _ = writeln!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="{FONT}">"#
    );
    let _ = writeln!(
        s,
        r#"<rect x="0" y="0" width="{width}" height="{height}" fill="{}"/>"#,
        rgb(p.canvas)
    );

    // --- branches -----------------------------------------------------------
    s.push_str("<g fill=\"none\" stroke-linecap=\"butt\">\n");
    for n in &l.nodes {
        let Some(parent) = n.parent.map(|i| &l.nodes[i]) else {
            continue;
        };
        let style = colors::focus_style(n.edge_focus, o.dark);
        line(&mut s, (sx(parent.x), sy(n.y)), (sx(n.x), sy(n.y)), &style);
    }
    for n in l.nodes.iter().filter(|n| n.leaf_count >= 2) {
        let style = colors::focus_style(n.focus, o.dark);
        line(
            &mut s,
            (sx(n.x), sy(n.row_min)),
            (sx(n.x), sy(n.row_max)),
            &style,
        );
    }
    s.push_str("</g>\n");

    // --- branch labels ------------------------------------------------------
    if o.branch_labels && row >= 10.0 {
        let fs = (row * 0.55).clamp(8.0, 11.0);
        let _ = writeln!(s, r#"<g font-size="{fs:.1}" text-anchor="middle">"#);
        for n in &l.nodes {
            let Some(parent) = n.parent.map(|i| &l.nodes[i]) else {
                continue;
            };
            if n.len <= 0.0 {
                continue;
            }
            let text = format!(
                "{}{}",
                if n.edge_approximate { "≈" } else { "" },
                fmt::snps(n.len)
            );
            let (x0, x1) = (sx(parent.x), sx(n.x));
            if x1 - x0 < text.chars().count() as f32 * fs * 0.56 + 6.0 {
                continue;
            }
            let _ = writeln!(
                s,
                r#"<text x="{:.1}" y="{:.1}" fill="{}">{}</text>"#,
                (x0 + x1) / 2.0,
                sy(n.y) - 2.0,
                rgb(colors::focus_color(n.edge_focus, o.dark)),
                esc(&text)
            );
        }
        s.push_str("</g>\n");
    }

    // --- nodes and tip labels -----------------------------------------------
    let _ = writeln!(s, r#"<g font-size="{label_px:.1}">"#);
    for n in &l.nodes {
        let (x, y) = (sx(n.x), sy(n.y));
        match &n.kind {
            NodeKind::Genome { id, name } => {
                let colour = model
                    .genome(*id)
                    .and_then(|g| g.cluster)
                    .map(|c| colors::cluster_color(c, o.dark))
                    .unwrap_or(p.text_weak);
                let _ = writeln!(
                    s,
                    r#"<circle cx="{x:.1}" cy="{y:.1}" r="2" fill="{}"/>"#,
                    rgb(colour)
                );
                let mut lx = if o.aligned_tips {
                    if tip_x - x > 8.0 {
                        let _ = writeln!(
                            s,
                            r#"<line x1="{:.1}" y1="{y:.1}" x2="{:.1}" y2="{y:.1}" stroke="{}" stroke-width="0.8" stroke-dasharray="1.2 2.8"/>"#,
                            x + 4.0,
                            tip_x - 3.0,
                            rgb(p.text_weak)
                        );
                    }
                    tip_x
                } else {
                    x + 5.0
                };
                if o.color_strip {
                    let h = (row * 0.86).max(2.0);
                    let _ = writeln!(
                        s,
                        r#"<rect x="{lx:.1}" y="{:.1}" width="7" height="{h:.1}" rx="1" fill="{}"/>"#,
                        y - h / 2.0,
                        rgb(colour)
                    );
                    lx += 11.0;
                }
                let text = if n.approximate {
                    format!("≈ {name}")
                } else {
                    name.clone()
                };
                let _ = writeln!(
                    s,
                    r#"<text x="{lx:.1}" y="{y:.1}" dominant-baseline="central" fill="{}">{}</text>"#,
                    rgb(p.text),
                    esc(&text)
                );
            }
            NodeKind::ClusterCollapsed { id, members, focus } => {
                let half = layout::collapsed_half_rows(*members);
                let t = layout::triangle(n.x, n.y, half);
                let colour = colors::cluster_color(*id, o.dark);
                let pts: Vec<String> = t
                    .iter()
                    .map(|&(tx, ty)| format!("{:.1},{:.1}", sx(tx), sy(ty)))
                    .collect();
                let _ = writeln!(
                    s,
                    r#"<polygon points="{}" fill="{}" fill-opacity="0.43" stroke="{}" stroke-width="1.2"/>"#,
                    pts.join(" "),
                    rgb(colour),
                    rgb(colour)
                );
                let approx = model.cluster(*id).map(|c| c.approximate).unwrap_or(true);
                let label = format!(
                    "c{id} ({members} genomes) {}{}",
                    fmt::focus_word(*focus),
                    if approx { " ≈" } else { "" }
                );
                let _ = writeln!(
                    s,
                    r#"<text x="{:.1}" y="{y:.1}" dominant-baseline="central" fill="{}">{}</text>"#,
                    sx(t[1].0) + 5.0,
                    rgb(colour),
                    esc(&label)
                );
            }
            NodeKind::ClusterRoot { id, members, focus } => {
                let colour = colors::cluster_color(*id, o.dark);
                let _ = writeln!(
                    s,
                    r#"<circle cx="{x:.1}" cy="{y:.1}" r="3" fill="{}"/>"#,
                    rgb(colour)
                );
                let label = format!(
                    "c{id} · {members} genomes · {}{}",
                    fmt::focus_word(*focus),
                    if n.approximate { " ≈" } else { "" }
                );
                let _ = writeln!(
                    s,
                    r#"<text x="{:.1}" y="{:.1}" font-size="{:.1}" fill="{}">{}</text>"#,
                    x + 7.0,
                    y - row * 0.6,
                    (label_px * 0.86).max(8.0),
                    rgb(colour),
                    esc(&label)
                );
            }
            NodeKind::Internal | NodeKind::Root => {}
        }
    }
    s.push_str("</g>\n");

    // --- scale bars ---------------------------------------------------------
    let mut y = top + l.rows as f32 * row + 20.0;
    let cluster_region = l
        .regions
        .iter()
        .find(|r| r.cluster.is_some() && r.cluster == o.selected_cluster)
        .or_else(|| l.regions.iter().find(|r| r.cluster.is_some()));
    let mut bars: Vec<(&layout::Region, String)> = Vec::new();
    if let Some(r) = cluster_region {
        bars.push((r, format!("in cluster c{}", r.cluster.unwrap_or(0))));
    }
    if let Some(bb) = l.regions.first().filter(|r| r.raw_max > 0.0) {
        bars.push((bb, "between clusters".to_string()));
    }
    for (region, what) in bars {
        let Some(bar) = layout::scale_bar(region, 0.18) else {
            continue;
        };
        let w = bar.layout_width as f32 * px_per_unit;
        if !w.is_finite() || w <= 4.0 {
            continue;
        }
        let style = colors::focus_style(region.focus, o.dark);
        let x0 = MARGIN;
        line(&mut s, (x0, y), (x0 + w, y), &style);
        let solid = colors::FocusStyle {
            dashed: false,
            ..style
        };
        line(&mut s, (x0, y - 4.0), (x0, y + 4.0), &solid);
        line(&mut s, (x0 + w, y - 4.0), (x0 + w, y + 4.0), &solid);
        let _ = writeln!(
            s,
            r#"<text x="{:.1}" y="{y:.1}" font-size="11" dominant-baseline="central" fill="{}">{}  ({what})</text>"#,
            x0 + w + 6.0,
            rgb(style.color),
            esc(&bar.label)
        );
        y += 18.0;
    }

    // --- legend -------------------------------------------------------------
    if o.legend {
        let fill = theme::mix(p.canvas, p.panel, 0.6);
        let _ = writeln!(
            s,
            r#"<rect x="{MARGIN}" y="{MARGIN}" width="214" height="62" rx="6" fill="{}" stroke="{}"/>"#,
            rgb(fill),
            rgb(p.border)
        );
        for (i, text) in [
            "Focus 0 — sketch only, ≈ SNPs",
            "Focus 1 — blocks compared, ≈ SNPs",
            "Focus 2 — exact SNPs",
        ]
        .iter()
        .enumerate()
        {
            let ly = MARGIN + 14.0 + i as f32 * 17.0;
            let style = colors::focus_style(i as u8, o.dark);
            line(&mut s, (MARGIN + 10.0, ly), (MARGIN + 40.0, ly), &style);
            let _ = writeln!(
                s,
                r#"<text x="{:.1}" y="{ly:.1}" font-size="11" dominant-baseline="central" fill="{}">{}</text>"#,
                MARGIN + 46.0,
                rgb(style.color),
                esc(text)
            );
        }
    }
    s.push_str("</svg>\n");
    s
}

/// One branch segment in a focus style.
fn line(s: &mut String, a: (f32, f32), b: (f32, f32), style: &colors::FocusStyle) {
    let dash = if style.dashed {
        format!(r#" stroke-dasharray="{0:.1} {0:.1}""#, style.dash.max(1.0))
    } else {
        String::new()
    };
    let _ = writeln!(
        s,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="{:.2}"{dash}/>"#,
        a.0,
        a.1,
        b.0,
        b.1,
        rgb(style.color),
        style.width
    );
}

fn rgb(c: Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}

/// Escape text for an SVG text node or attribute.
pub fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use infinity_core::events::{ClusterFocusInfo, Event};
    use std::collections::HashSet;

    fn model() -> UiModel {
        let mut m = UiModel::default();
        for i in 0..4u32 {
            m.apply(0.0, crate::state::scan_event(i, &format!("g{i}"), 1000, 1));
        }
        m.apply(
            0.0,
            crate::state::clustered_event(&[(0, 0, &[0, 1, 2]), (1, 3, &[3])]),
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
        m
    }

    #[test]
    fn the_svg_is_a_complete_document_with_every_tip_and_the_legend() {
        let m = model();
        let open: HashSet<u32> = [0u32].into_iter().collect();
        let l = layout::layout(&m.display_tree(&open).unwrap());
        let svg = render(&l, &m, &SvgOptions::default());
        assert!(svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""));
        assert!(svg.trim_end().ends_with("</svg>"));
        for name in ["g0", "g1", "g2"] {
            assert!(svg.contains(&format!(">{name}</text>")), "{name} missing");
        }
        // The singleton cluster is still collapsed: a triangle with its label.
        assert!(svg.contains("<polygon"));
        assert!(svg.contains("c1 (1 genomes)"));
        // Focus 0 backbone branches are dashed, the exact cluster's are not.
        assert!(svg.contains("stroke-dasharray"));
        let exact = rgb(colors::focus_color(2, false));
        assert!(svg.contains(&format!("stroke=\"{exact}\" stroke-width=\"1.70\"/>")));
        // Legend and a SNP scale bar.
        assert!(svg.contains("Focus 2 — exact SNPs"));
        assert!(svg.contains("SNPs  (in cluster c0)"), "{svg}");
        // Tags are balanced: every opened group is closed.
        assert_eq!(svg.matches("<g").count(), svg.matches("</g>").count());
    }

    #[test]
    fn options_add_branch_labels_strips_and_leaders() {
        let m = model();
        let open: HashSet<u32> = [0u32].into_iter().collect();
        let l = layout::layout(&m.display_tree(&open).unwrap());
        let plain = render(
            &l,
            &m,
            &SvgOptions {
                color_strip: false,
                legend: false,
                ..Default::default()
            },
        );
        assert!(!plain.contains("rx=\"1\""), "no colour strip");
        assert!(!plain.contains("Focus 0 —"), "no legend");
        let rich = render(
            &l,
            &m,
            &SvgOptions {
                branch_labels: true,
                aligned_tips: true,
                row_px: 20.0,
                ..Default::default()
            },
        );
        assert!(rich.contains("rx=\"1\""), "colour strip");
        assert!(rich.contains("stroke-dasharray=\"1.2 2.8\""), "leaders");
        assert!(rich.contains(">4</text>"), "the 4-SNP branch is labelled");
    }

    #[test]
    fn text_is_escaped() {
        assert_eq!(esc("a<b & \"c\"'"), "a&lt;b &amp; &quot;c&quot;&apos;");
    }
}
