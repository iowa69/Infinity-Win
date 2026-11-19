//! Colours and line styles: the visual vocabulary of focus levels, clusters and
//! variant kinds (section 8).
//!
//! Focus must always be distinguishable, because a Focus 0 number and a
//! Focus 2 number are different kinds of number and rule 8 of section 0 forbids
//! letting one pass for the other. The three styles differ in colour, weight
//! and dashing; the genome list and the pairwise readout take their colours
//! from here. The tree itself is drawn as a publication figure
//! ([`crate::pubtree`]): exact branches in ink, estimates in grey, every
//! estimated number with `≈`.
//!
//! `Color32` is plain data, so this module is testable without a window.

use egui::Color32;
use infinity_core::variant::VariantKind;

use crate::fmt;

/// How one focus level is drawn.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FocusStyle {
    pub color: Color32,
    pub width: f32,
    /// Focus 0 branches are dashed (section 8).
    pub dashed: bool,
    /// Dash length in screen pixels.
    pub dash: f32,
}

/// Focus 0: sketch only. Dashed and pale.
pub fn focus0(dark: bool) -> FocusStyle {
    FocusStyle {
        color: if dark {
            Color32::from_rgb(120, 132, 152)
        } else {
            Color32::from_rgb(168, 180, 200)
        },
        width: 0.9,
        dashed: true,
        dash: 3.0,
    }
}

/// Focus 1: blocks compared, differences known but not resolved. Solid, thin,
/// orange, and always accompanied by the "≈" marker.
pub fn focus1(dark: bool) -> FocusStyle {
    FocusStyle {
        color: if dark {
            Color32::from_rgb(228, 158, 62)
        } else {
            Color32::from_rgb(206, 132, 30)
        },
        width: 1.1,
        dashed: false,
        dash: 0.0,
    }
}

/// Focus 2: exact, base level. Solid and strong. "Dark" means high contrast
/// against the background, which is near-white on the dark theme.
pub fn focus2(dark: bool) -> FocusStyle {
    FocusStyle {
        color: if dark {
            Color32::from_rgb(228, 232, 238)
        } else {
            Color32::from_rgb(34, 38, 46)
        },
        width: 1.7,
        dashed: false,
        dash: 0.0,
    }
}

/// The style of a branch at a given focus level.
pub fn focus_style(focus: u8, dark: bool) -> FocusStyle {
    match focus {
        fmt::FOCUS_SKETCH => focus0(dark),
        fmt::FOCUS_BLOCKS => focus1(dark),
        fmt::FOCUS_EXACT => focus2(dark),
        _ => focus0(dark),
    }
}

/// Just the colour of a focus level, for text and swatches.
pub fn focus_color(focus: u8, dark: bool) -> Color32 {
    focus_style(focus, dark).color
}

/// Convert HSV to an opaque RGB colour. `h` is wrapped into `[0, 1)`.
pub fn hsv(h: f32, s: f32, v: f32) -> Color32 {
    let h = (h - h.floor()).clamp(0.0, 1.0);
    let s = s.clamp(0.0, 1.0);
    let v = v.clamp(0.0, 1.0);
    let i = (h * 6.0).floor() as i32;
    let f = h * 6.0 - i as f32;
    let p = v * (1.0 - s);
    let q = v * (1.0 - f * s);
    let t = v * (1.0 - (1.0 - f) * s);
    let (r, g, b) = match i.rem_euclid(6) {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    Color32::from_rgb(
        (r * 255.0).round() as u8,
        (g * 255.0).round() as u8,
        (b * 255.0).round() as u8,
    )
}

/// The golden-ratio hue step, so that neighbouring cluster ids never look alike.
const GOLDEN: f32 = 0.618_034;

/// A stable colour for a cluster. The same id always gives the same colour
/// within a session, which is what makes the genome list swatch, the tree
/// triangle and the map row recognisable as the same cluster.
pub fn cluster_color(id: u32, dark: bool) -> Color32 {
    let h = (id as f32) * GOLDEN;
    if dark {
        hsv(h.fract(), 0.52, 0.92)
    } else {
        hsv(h.fract(), 0.62, 0.72)
    }
}

/// A pale wash of a cluster colour, for selected and hovered backgrounds.
pub fn wash(c: Color32, alpha: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), alpha)
}

/// Identical stretches: a thin grey line (section 8, map panel).
pub const IDENTICAL: Color32 = Color32::from_rgb(150, 150, 150);

/// A difference piece that has not been resolved yet: an orange box.
pub const PIECE_UNRESOLVED: Color32 = Color32::from_rgb(235, 140, 40);

/// A difference piece that has been resolved: a paler, narrower orange box.
/// Drawing both widths is the only part of the section 8 "zoom made visible"
/// requirement the event stream can support (see `map_panel`).
pub const PIECE_RESOLVED: Color32 = Color32::from_rgb(240, 190, 120);

/// A breakpoint: a triangle.
pub const BREAKPOINT: Color32 = Color32::from_rgb(150, 90, 30);

/// The colour of each variant kind on the map.
pub fn variant_color(kind: VariantKind) -> Color32 {
    match kind {
        VariantKind::Snp => Color32::from_rgb(214, 40, 40),
        VariantKind::Insertion | VariantKind::Deletion => Color32::from_rgb(45, 110, 220),
        VariantKind::Complex => Color32::from_rgb(150, 60, 190),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_three_focus_styles_differ_in_every_channel() {
        let light = (
            focus_style(0, false),
            focus_style(1, false),
            focus_style(2, false),
        );
        assert_ne!(light.0.color, light.1.color);
        assert_ne!(light.1.color, light.2.color);
        assert_ne!(light.0.width, light.1.width);
        assert_ne!(light.1.width, light.2.width);
        assert!(light.0.dashed && !light.1.dashed && !light.2.dashed);
        // Focus 0 is the palest, Focus 2 the strongest, on the light theme.
        let lum = |c: Color32| c.r() as i32 + c.g() as i32 + c.b() as i32;
        assert!(lum(light.0.color) > lum(light.2.color));
    }

    #[test]
    fn focus_styles_adapt_to_the_theme() {
        // On a dark background Focus 2 must be light, or it disappears.
        let dark2 = focus_style(2, true).color;
        assert!(dark2.r() > 150 && dark2.g() > 150 && dark2.b() > 150);
        let light2 = focus_style(2, false).color;
        assert!(light2.r() < 100);
    }

    #[test]
    fn hsv_primary_colours_are_exact() {
        assert_eq!(hsv(0.0, 1.0, 1.0), Color32::from_rgb(255, 0, 0));
        assert_eq!(hsv(1.0 / 3.0, 1.0, 1.0), Color32::from_rgb(0, 255, 0));
        assert_eq!(hsv(2.0 / 3.0, 1.0, 1.0), Color32::from_rgb(0, 0, 255));
        assert_eq!(hsv(0.5, 0.0, 1.0), Color32::from_rgb(255, 255, 255));
        // Hue wraps rather than panicking.
        assert_eq!(hsv(1.0, 1.0, 1.0), hsv(0.0, 1.0, 1.0));
        assert_eq!(hsv(-0.25, 1.0, 1.0), hsv(0.75, 1.0, 1.0));
    }

    #[test]
    fn cluster_colours_are_stable_and_distinct() {
        assert_eq!(cluster_color(7, false), cluster_color(7, false));
        assert_ne!(cluster_color(0, false), cluster_color(1, false));
        assert_ne!(cluster_color(1, false), cluster_color(2, false));
        // No two of the first 24 clusters may collide: the swatches are the
        // only thing tying a genome list row to a tree triangle.
        let mut seen = Vec::new();
        for id in 0..24u32 {
            let c = cluster_color(id, false);
            assert!(!seen.contains(&c), "cluster {id} reuses a colour");
            seen.push(c);
        }
    }

    #[test]
    fn variant_kinds_have_the_colours_the_spec_names() {
        assert_eq!(variant_color(VariantKind::Snp).r(), 214);
        assert_eq!(
            variant_color(VariantKind::Insertion),
            variant_color(VariantKind::Deletion)
        );
        assert_ne!(
            variant_color(VariantKind::Snp),
            variant_color(VariantKind::Complex)
        );
    }
}
