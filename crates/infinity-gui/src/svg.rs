//! The tree as a vector picture (SVG): the publication-style figure of
//! [`crate::pubtree`], written primitive by primitive, on white.
//!
//! Where the PNG export is a screenshot of what is on screen, this draws the
//! *whole* laid-out tree at a readable row height, so a tree of three hundred
//! genomes comes out as a figure rather than as a crop. The tree window paints
//! the very same primitives, so its "Save SVG" is what it shows. Pure string
//! building: testable without a window.

use std::fmt::Write as _;

use egui::Color32;

use crate::pubtree::{Figure, HAlign, Ink, Prim, PubOptions};

/// What the picture should look like.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SvgOptions {
    pub opts: PubOptions,
    /// Width of the tree drawing itself, without the tip labels, in px.
    pub tree_width: f32,
    /// Height of one leaf row, in px.
    pub row_px: f32,
}

impl Default for SvgOptions {
    fn default() -> Self {
        Self {
            opts: PubOptions::default(),
            tree_width: 720.0,
            row_px: 16.0,
        }
    }
}

const FONT: &str = "Helvetica, Arial, sans-serif";

/// The whole figure as a complete SVG document, on white.
pub fn render(f: &Figure, o: &SvgOptions) -> String {
    let (prims, size) = crate::pubtree::whole_scene(
        f,
        &o.opts,
        &Ink::PAPER,
        o.row_px.clamp(4.0, 60.0),
        o.tree_width.max(100.0),
    );
    from_scene(&prims, size)
}

/// Write a list of drawing primitives as a complete SVG document of `size`
/// pixels on a white page. The windows paint the same list, so the file is
/// what they show.
pub fn from_scene(prims: &[Prim], size: egui::Vec2) -> String {
    let (w, h) = (size.x.ceil(), size.y.ceil());
    let mut s = String::with_capacity(96 * prims.len() + 512);
    let _ = writeln!(
        s,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="{FONT}">"#
    );
    let _ = writeln!(
        s,
        r##"<rect x="0" y="0" width="{w}" height="{h}" fill="#ffffff"/>"##
    );
    for p in prims {
        match p {
            Prim::Line { a, b, width, color } => {
                let _ = writeln!(
                    s,
                    r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}"{} stroke-width="{width:.2}"/>"#,
                    a.x,
                    a.y,
                    b.x,
                    b.y,
                    rgbu(*color),
                    opacity("stroke-opacity", *color)
                );
            }
            Prim::Circle { c, r, fill, stroke } => {
                let st = stroke
                    .map(|(w, col)| format!(r#" stroke="{}" stroke-width="{w:.2}""#, rgbu(col)))
                    .unwrap_or_default();
                let fill_attr = if fill.a() == 0 {
                    r#"fill="none""#.to_string()
                } else {
                    format!(
                        r#"fill="{}"{}"#,
                        rgbu(*fill),
                        opacity("fill-opacity", *fill)
                    )
                };
                let _ = writeln!(
                    s,
                    r#"<circle cx="{:.2}" cy="{:.2}" r="{r:.2}" {fill_attr}{st}/>"#,
                    c.x, c.y,
                );
            }
            Prim::Rect { r, fill, rounding } => {
                let _ = writeln!(
                    s,
                    r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" rx="{rounding:.1}" fill="{}"{}/>"#,
                    r.left(),
                    r.top(),
                    r.width(),
                    r.height(),
                    rgbu(*fill),
                    opacity("fill-opacity", *fill)
                );
            }
            Prim::Poly { pts, fill, stroke } => {
                let pts: Vec<String> = pts
                    .iter()
                    .map(|p| format!("{:.2},{:.2}", p.x, p.y))
                    .collect();
                let _ = writeln!(
                    s,
                    r#"<polygon points="{}" fill="{}"{} stroke="{}" stroke-width="{:.2}" stroke-linejoin="round"/>"#,
                    pts.join(" "),
                    rgbu(*fill),
                    opacity("fill-opacity", *fill),
                    rgbu(stroke.1),
                    stroke.0
                );
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
                    HAlign::Left => "start",
                    HAlign::Center => "middle",
                };
                let weight = if *bold { r#" font-weight="bold""# } else { "" };
                let tail = tail
                    .as_ref()
                    .map(|t| {
                        format!(
                            r#"<tspan dx="{:.1}" font-size="{:.1}" font-weight="normal" fill="{}">{}</tspan>"#,
                            crate::pubtree::tail_gap(*size),
                            t.size,
                            rgbu(t.color),
                            esc(&t.text)
                        )
                    })
                    .unwrap_or_default();
                let _ = writeln!(
                    s,
                    r#"<text x="{:.2}" y="{:.2}" font-size="{size:.1}" text-anchor="{anchor}" dominant-baseline="central" fill="{}"{weight}>{}{tail}</text>"#,
                    pos.x,
                    pos.y,
                    rgbu(*color),
                    esc(text)
                );
            }
        }
    }
    s.push_str("</svg>\n");
    s
}

/// An SVG opacity attribute for a colour that is not fully opaque.
fn opacity(attr: &str, c: Color32) -> String {
    let a = c.to_srgba_unmultiplied()[3];
    if a == 255 {
        String::new()
    } else {
        format!(r#" {attr}="{:.3}""#, f32::from(a) / 255.0)
    }
}

/// The colour itself, without the premultiplied alpha `Color32` stores: a
/// translucent wash keeps its hue and gets its alpha as an opacity attribute.
fn rgbu(c: Color32) -> String {
    let [r, g, b, _] = c.to_srgba_unmultiplied();
    format!("#{r:02x}{g:02x}{b:02x}")
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

    #[test]
    fn the_svg_is_a_complete_white_document_with_every_tip_and_its_snps() {
        let f = crate::pubtree::tests_support::figure();
        let svg = render(&f, &SvgOptions::default());
        assert!(svg.starts_with("<svg xmlns=\"http://www.w3.org/2000/svg\""));
        assert!(svg.trim_end().ends_with("</svg>"));
        assert!(svg.contains(r##"fill="#ffffff""##), "white page");
        for name in ["g0", "g1", "g2"] {
            assert!(svg.contains(&format!(">{name}<tspan")), "{name} missing");
        }
        // The tip numbers ride along as a lighter tspan.
        assert!(svg.contains(">+4</tspan>"), "{svg}");
        // No dashes anywhere: estimates are grey, not dotted.
        assert!(!svg.contains("stroke-dasharray"));
        // The SNP group band and its label.
        assert!(svg.contains("G1 · 2 genomes ≤5 SNPs"));
        // A scale bar in SNPs.
        assert!(svg.contains("SNPs within clusters"));
        // The selected genome is bold.
        assert!(svg.contains(r#"font-weight="bold">g1<tspan"#));
    }

    #[test]
    fn options_switch_decorations_off() {
        let f = crate::pubtree::tests_support::figure();
        let plain = render(
            &f,
            &SvgOptions {
                opts: PubOptions {
                    tip_snps: false,
                    groups: false,
                    legend: false,
                    cluster_brackets: false,
                    ..Default::default()
                },
                ..Default::default()
            },
        );
        assert!(!plain.contains("<tspan"), "no tip numbers");
        assert!(!plain.contains("G1 ·"));
        assert!(!plain.contains("exact SNPs<"), "no legend");
        let rich = render(
            &f,
            &SvgOptions {
                opts: PubOptions {
                    branch_labels: true,
                    aligned_tips: true,
                    color_strip: true,
                    ..Default::default()
                },
                row_px: 20.0,
                ..Default::default()
            },
        );
        assert!(rich.contains(">1</text>"), "the 1-SNP internal branch");
    }

    #[test]
    fn text_is_escaped() {
        assert_eq!(esc("a<b & \"c\"'"), "a&lt;b &amp; &quot;c&quot;&apos;");
    }
}
