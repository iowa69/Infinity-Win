//! One visual language for the whole window: palette, spacing, typography and
//! the few building blocks every panel uses (panel headers, section titles,
//! primary buttons, labelled rows with plain-word tooltips).
//!
//! Both the dark and the light theme are installed up front, so switching
//! between them (View menu, or following the system) never leaves a panel in
//! the other theme's colours. Panels ask [`palette`] for colours instead of
//! hard-coding them, which is what keeps the two themes consistent.
//!
//! The accent is a cool cyan: it marks what is interactive or selected, the
//! panel titles, and the "live" parts of the window (the compute view uses it
//! as its glow). Focus levels and cluster colours keep their own vocabulary in
//! [`crate::colors`], because they carry meaning rather than decoration.

use eframe::egui;
use egui::{Color32, CornerRadius, FontId, Margin, RichText, Stroke, TextStyle, Theme, Vec2};

/// The colours of one theme.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Palette {
    pub dark: bool,
    /// Behind drawing canvases (tree, map, compute view) and text fields.
    pub canvas: Color32,
    /// Panel backgrounds.
    pub panel: Color32,
    /// Resting widget background.
    pub surface: Color32,
    /// Hovered widget background.
    pub surface_hover: Color32,
    /// Hairlines between panels and around boxes.
    pub border: Color32,
    pub text: Color32,
    pub text_strong: Color32,
    pub text_weak: Color32,
    /// Interactive, selected, live.
    pub accent: Color32,
    /// A deeper accent, for filled primary buttons.
    pub accent_fill: Color32,
    /// Text on an accent fill.
    pub on_accent: Color32,
    pub warn: Color32,
    pub error: Color32,
    pub ok: Color32,
}

pub const DARK: Palette = Palette {
    dark: true,
    canvas: Color32::from_rgb(12, 16, 22),
    panel: Color32::from_rgb(19, 24, 32),
    surface: Color32::from_rgb(30, 38, 50),
    surface_hover: Color32::from_rgb(40, 51, 67),
    border: Color32::from_rgb(44, 55, 71),
    text: Color32::from_rgb(208, 217, 228),
    text_strong: Color32::from_rgb(240, 245, 250),
    text_weak: Color32::from_rgb(128, 142, 160),
    accent: Color32::from_rgb(64, 216, 222),
    accent_fill: Color32::from_rgb(22, 140, 150),
    on_accent: Color32::from_rgb(240, 252, 252),
    warn: Color32::from_rgb(232, 166, 64),
    error: Color32::from_rgb(240, 108, 108),
    ok: Color32::from_rgb(112, 210, 140),
};

pub const LIGHT: Palette = Palette {
    dark: false,
    canvas: Color32::from_rgb(252, 253, 254),
    panel: Color32::from_rgb(243, 245, 248),
    surface: Color32::from_rgb(228, 233, 239),
    surface_hover: Color32::from_rgb(214, 222, 231),
    border: Color32::from_rgb(206, 213, 223),
    text: Color32::from_rgb(36, 44, 56),
    text_strong: Color32::from_rgb(12, 18, 26),
    text_weak: Color32::from_rgb(100, 112, 126),
    accent: Color32::from_rgb(10, 132, 140),
    accent_fill: Color32::from_rgb(12, 128, 136),
    on_accent: Color32::from_rgb(255, 255, 255),
    warn: Color32::from_rgb(178, 104, 0),
    error: Color32::from_rgb(196, 56, 56),
    ok: Color32::from_rgb(30, 140, 70),
};

/// The palette of the theme a `Ui` is drawn in.
pub fn palette(ui: &egui::Ui) -> Palette {
    palette_for(ui.visuals().dark_mode)
}

pub fn palette_for(dark: bool) -> Palette {
    if dark {
        DARK
    } else {
        LIGHT
    }
}

/// Install both themes' styles. Called once at start-up; harmless to repeat.
pub fn install(ctx: &egui::Context) {
    ctx.style_mut_of(Theme::Dark, |s| apply(s, &DARK));
    ctx.style_mut_of(Theme::Light, |s| apply(s, &LIGHT));
}

/// Which theme to follow, as remembered in the preferences.
pub fn preference(name: &str) -> egui::ThemePreference {
    match name {
        "light" => egui::ThemePreference::Light,
        "system" => egui::ThemePreference::System,
        _ => egui::ThemePreference::Dark,
    }
}

fn apply(style: &mut egui::Style, p: &Palette) {
    // --- typography: one scale, used everywhere -----------------------------
    style.text_styles = [
        (TextStyle::Heading, FontId::proportional(16.0)),
        (TextStyle::Body, FontId::proportional(13.0)),
        (TextStyle::Button, FontId::proportional(13.0)),
        (TextStyle::Small, FontId::proportional(10.5)),
        (TextStyle::Monospace, FontId::monospace(12.0)),
    ]
    .into();

    // --- spacing ------------------------------------------------------------
    let sp = &mut style.spacing;
    sp.item_spacing = Vec2::new(7.0, 5.0);
    sp.button_padding = Vec2::new(9.0, 3.0);
    sp.interact_size = Vec2::new(36.0, 22.0);
    sp.menu_margin = Margin::same(6);
    sp.window_margin = Margin::same(12);
    sp.indent = 16.0;
    sp.slider_width = 120.0;
    sp.combo_width = 120.0;
    sp.text_edit_width = 200.0;
    sp.tooltip_width = 360.0;

    // --- colours ------------------------------------------------------------
    let v = &mut style.visuals;
    *v = if p.dark {
        egui::Visuals::dark()
    } else {
        egui::Visuals::light()
    };
    v.panel_fill = p.panel;
    v.window_fill = p.panel;
    v.extreme_bg_color = p.canvas;
    v.text_edit_bg_color = Some(p.canvas);
    v.faint_bg_color = if p.dark {
        Color32::from_rgb(24, 30, 40)
    } else {
        Color32::from_rgb(236, 240, 245)
    };
    v.code_bg_color = p.surface;
    v.warn_fg_color = p.warn;
    v.error_fg_color = p.error;
    v.hyperlink_color = p.accent;
    v.window_stroke = Stroke::new(1.0, p.border);
    v.window_corner_radius = CornerRadius::same(8);
    v.menu_corner_radius = CornerRadius::same(6);
    v.selection.bg_fill = with_alpha(p.accent, if p.dark { 70 } else { 60 });
    v.selection.stroke = Stroke::new(1.0, p.accent);
    v.slider_trailing_fill = true;
    v.striped = true;
    v.indent_has_left_vline = false;

    let r = CornerRadius::same(5);
    let w = &mut v.widgets;
    w.noninteractive.bg_fill = p.panel;
    w.noninteractive.weak_bg_fill = p.panel;
    w.noninteractive.bg_stroke = Stroke::new(1.0, p.border);
    w.noninteractive.fg_stroke = Stroke::new(1.0, p.text);
    w.noninteractive.corner_radius = r;

    w.inactive.bg_fill = p.surface;
    w.inactive.weak_bg_fill = p.surface;
    w.inactive.bg_stroke = Stroke::NONE;
    w.inactive.fg_stroke = Stroke::new(1.0, p.text);
    w.inactive.corner_radius = r;

    w.hovered.bg_fill = p.surface_hover;
    w.hovered.weak_bg_fill = p.surface_hover;
    w.hovered.bg_stroke = Stroke::new(1.0, with_alpha(p.accent, 140));
    w.hovered.fg_stroke = Stroke::new(1.2, p.text_strong);
    w.hovered.corner_radius = r;
    w.hovered.expansion = 0.0;

    w.active.bg_fill = p.accent_fill;
    w.active.weak_bg_fill = p.accent_fill;
    w.active.bg_stroke = Stroke::new(1.0, p.accent);
    w.active.fg_stroke = Stroke::new(1.5, p.on_accent);
    w.active.corner_radius = r;
    w.active.expansion = 0.0;

    w.open.bg_fill = p.surface_hover;
    w.open.weak_bg_fill = p.surface_hover;
    w.open.bg_stroke = Stroke::new(1.0, p.accent);
    w.open.fg_stroke = Stroke::new(1.0, p.text_strong);
    w.open.corner_radius = r;
}

/// A colour with a new alpha (unmultiplied).
pub fn with_alpha(c: Color32, a: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(c.r(), c.g(), c.b(), a)
}

/// Linear blend of two opaque colours; `t = 0` is `a`, `t = 1` is `b`.
pub fn mix(a: Color32, b: Color32, t: f32) -> Color32 {
    let t = if t.is_finite() {
        t.clamp(0.0, 1.0)
    } else {
        0.0
    };
    let l = |x: u8, y: u8| (f32::from(x) + (f32::from(y) - f32::from(x)) * t).round() as u8;
    Color32::from_rgba_unmultiplied(
        l(a.r(), b.r()),
        l(a.g(), b.g()),
        l(a.b(), b.b()),
        l(a.a(), b.a()),
    )
}

// ---------------------------------------------------------------------------
// Building blocks
// ---------------------------------------------------------------------------

/// The title of a panel: small, spaced capitals in the accent colour.
pub fn title(text: &str, p: &Palette) -> RichText {
    RichText::new(text.to_uppercase())
        .size(11.5)
        .strong()
        .extra_letter_spacing(1.2)
        .color(p.accent)
}

/// A panel header: the title, then whatever the caller adds on the same
/// (wrapping) line, then a hairline. Every panel starts with one of these.
pub fn panel_header<R>(ui: &mut egui::Ui, text: &str, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let p = palette(ui);
    let r = ui
        .horizontal_wrapped(|ui| {
            ui.label(title(text, &p));
            add(ui)
        })
        .inner;
    hairline(ui);
    r
}

/// A thin full-width rule in the border colour.
pub fn hairline(ui: &mut egui::Ui) {
    let p = palette(ui);
    let (rect, _) =
        ui.allocate_exact_size(Vec2::new(ui.available_width(), 5.0), egui::Sense::hover());
    ui.painter().line_segment(
        [rect.left_center(), rect.right_center()],
        Stroke::new(1.0, p.border),
    );
}

/// A section title inside a panel or the settings window, with an accent tick.
pub fn section(ui: &mut egui::Ui, text: &str) {
    let p = palette(ui);
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        let (rect, _) = ui.allocate_exact_size(Vec2::new(3.0, 13.0), egui::Sense::hover());
        ui.painter().rect_filled(rect, 1, p.accent);
        ui.label(RichText::new(text).size(12.5).strong().color(p.text_strong));
    });
}

/// Small explanatory text under a control.
pub fn hint(ui: &mut egui::Ui, text: &str) {
    let p = palette(ui);
    ui.label(RichText::new(text).small().color(p.text_weak));
}

/// The one filled button per area: the thing the user most likely wants next.
pub fn primary_button(text: &str, p: &Palette) -> egui::Button<'static> {
    egui::Button::new(RichText::new(text.to_string()).strong().color(p.on_accent))
        .fill(p.accent_fill)
        .stroke(Stroke::new(1.0, p.accent))
}

/// A two-column row of a settings grid: a label with a plain-word tooltip on
/// the left, the control on the right. Must be called inside an `egui::Grid`.
pub fn grid_row<R>(
    ui: &mut egui::Ui,
    label: &str,
    tip: &str,
    add: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    ui.label(label).on_hover_text(tip);
    let r = ui.horizontal(|ui| add(ui)).inner;
    ui.end_row();
    r
}

/// A label and a value on one line, the value right-aligned: the counters.
///
/// When the two do not fit side by side (a narrow panel, a long value), the
/// value moves to its own right-aligned line instead of overlapping the label.
pub fn stat_row(ui: &mut egui::Ui, label: &str, value: &str) {
    let p = palette(ui);
    let font = TextStyle::Body.resolve(ui.style());
    let need = ui.ctx().fonts_mut(|f| {
        f.layout_no_wrap(label.to_string(), font.clone(), p.text)
            .size()
            .x
            + f.layout_no_wrap(value.to_string(), font.clone(), p.text)
                .size()
                .x
    }) + 16.0;
    if need <= ui.available_width() {
        ui.horizontal(|ui| {
            ui.label(RichText::new(label).color(p.text_weak));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(RichText::new(value).color(p.text));
            });
        });
    } else {
        ui.label(RichText::new(label).color(p.text_weak));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(RichText::new(value).color(p.text));
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mixing_hits_both_ends_and_the_middle() {
        let a = Color32::from_rgb(0, 100, 200);
        let b = Color32::from_rgb(200, 100, 0);
        assert_eq!(mix(a, b, 0.0), a);
        assert_eq!(mix(a, b, 1.0), b);
        assert_eq!(mix(a, b, 0.5), Color32::from_rgb(100, 100, 100));
        // Out of range and NaN never produce garbage.
        assert_eq!(mix(a, b, 7.0), b);
        assert_eq!(mix(a, b, f32::NAN), a);
    }

    #[test]
    fn the_two_palettes_are_each_readable() {
        // Text must contrast with the panel in both themes.
        let lum =
            |c: Color32| 0.2126 * c.r() as f32 + 0.7152 * c.g() as f32 + 0.0722 * c.b() as f32;
        assert!(lum(DARK.text) - lum(DARK.panel) > 120.0);
        assert!(lum(LIGHT.panel) - lum(LIGHT.text) > 120.0);
        // Weak text is weaker than text, but still readable.
        assert!(lum(DARK.text_weak) < lum(DARK.text));
        assert!(lum(DARK.text_weak) - lum(DARK.panel) > 60.0);
        assert!(lum(LIGHT.text_weak) > lum(LIGHT.text));
        assert!(lum(LIGHT.panel) - lum(LIGHT.text_weak) > 60.0);
    }

    #[test]
    fn theme_names_map_to_preferences() {
        assert_eq!(preference("light"), egui::ThemePreference::Light);
        assert_eq!(preference("system"), egui::ThemePreference::System);
        assert_eq!(preference("dark"), egui::ThemePreference::Dark);
        assert_eq!(preference("anything"), egui::ThemePreference::Dark);
    }

    #[test]
    fn installing_sets_both_themes() {
        let ctx = egui::Context::default();
        install(&ctx);
        ctx.set_theme(egui::ThemePreference::Dark);
        assert_eq!(ctx.global_style().visuals.panel_fill, DARK.panel);
        ctx.set_theme(egui::ThemePreference::Light);
        assert_eq!(ctx.global_style().visuals.panel_fill, LIGHT.panel);
    }
}
