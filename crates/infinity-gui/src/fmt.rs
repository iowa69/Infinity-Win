//! Number and label formatting shared by every panel.
//!
//! Section 0 rule 8 ("never present an approximate number as exact") is
//! enforced here rather than at each call site: every formatter that can be
//! handed an estimated value takes an `approximate` flag and prefixes the
//! result with "≈", so a panel cannot forget the marker.

use infinity_core::variant::VariantKind;

/// Focus levels of section 5.
pub const FOCUS_SKETCH: u8 = 0;
pub const FOCUS_BLOCKS: u8 = 1;
pub const FOCUS_EXACT: u8 = 2;

/// The "≈" marker required by section 0 rule 8, or the empty string.
pub fn approx_mark(approximate: bool) -> &'static str {
    if approximate {
        "≈ "
    } else {
        ""
    }
}

/// Convert a 0-based half-open start into the 1-based position a user reads
/// (section 5: the conversion happens only where a position is displayed).
///
/// An insertion sits *between* two bases, so its position is the base before
/// the insertion site; adding one would point at the base after it. This
/// mirrors `output.rs`, so a position on screen and the same position in
/// `cluster_<id>.variants.tsv` always agree.
pub fn one_based(start: u64, kind: VariantKind) -> u64 {
    match kind {
        VariantKind::Insertion => start,
        VariantKind::Snp | VariantKind::Deletion | VariantKind::Complex => start.saturating_add(1),
    }
}

/// A 0-based half-open interval as an inclusive 1-based range, for display.
pub fn one_based_range(start: u64, end: u64, kind: VariantKind) -> String {
    let s = one_based(start, kind);
    let e = one_based(end.saturating_sub(1), kind);
    if e <= s {
        format!("{s}")
    } else {
        format!("{s}-{e}")
    }
}

/// An integer with thin-space thousands separators, for counters.
pub fn count(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push('\u{2009}');
        }
        out.push(c);
    }
    out
}

/// A base-pair length, scaled to the unit that keeps it readable.
pub fn bp(n: u64) -> String {
    if n < 1_000 {
        format!("{n} bp")
    } else if n < 1_000_000 {
        format!("{:.1} kb", n as f64 / 1e3)
    } else if n < 1_000_000_000 {
        format!("{:.2} Mb", n as f64 / 1e6)
    } else {
        format!("{:.2} Gb", n as f64 / 1e9)
    }
}

/// Milliseconds as a human duration.
pub fn secs(ms: u64) -> String {
    if ms < 1_000 {
        format!("{ms} ms")
    } else if ms < 60_000 {
        format!("{:.1} s", ms as f64 / 1000.0)
    } else {
        let s = ms / 1000;
        format!("{} m {:02} s", s / 60, s % 60)
    }
}

/// A fraction in `[0, 1]` as a percentage. `None` stays "n/a": a ratio that
/// was never measured is not zero (section 0 rule 8).
pub fn pct(x: Option<f64>) -> String {
    match x {
        Some(v) if v.is_finite() => format!("{:.2} %", v * 100.0),
        _ => "n/a".to_string(),
    }
}

/// A distance with its focus-dependent unit, prefixed with "≈" when the number
/// behind it is an estimate rather than a measurement.
pub fn distance(value: f64, approximate: bool) -> String {
    if !value.is_finite() {
        return "n/a".to_string();
    }
    format!("{}{}", approx_mark(approximate), sci(value))
}

/// Scientific notation with three significant digits, or a plain integer when
/// the value is a small whole number (SNP counts read better that way).
pub fn sci(v: f64) -> String {
    if !v.is_finite() {
        return "n/a".to_string();
    }
    if v == 0.0 {
        return "0".to_string();
    }
    if v.abs() >= 1e-4 && v.abs() < 1e6 && (v - v.round()).abs() < 1e-9 {
        return format!("{}", v.round() as i64);
    }
    format!("{v:.3e}")
}

/// The word the mock-up of section 8 uses for a focus level.
pub fn focus_word(focus: u8) -> &'static str {
    match focus {
        FOCUS_SKETCH => "blurry",
        FOCUS_BLOCKS => "medium",
        FOCUS_EXACT => "sharp",
        _ => "unknown",
    }
}

/// "Focus 2 (sharp)", the label used in the genome list and the tooltips.
pub fn focus_label(focus: u8, approximate: bool) -> String {
    let mark = if approximate { "≈" } else { "" };
    format!("Focus {focus}{mark} ({})", focus_word(focus))
}

/// The unit every tree is *drawn* in. Inside a cluster the engine's Focus 2
/// tree is neighbour joining on SNP counts, so a phylogram in SNPs is exact
/// there. Everything else is converted into the same unit and carries the "≈"
/// marker (see [`source_unit`] for where each estimate comes from), so one
/// scale-bar vocabulary serves the whole picture.
pub const TREE_UNIT: &str = "SNPs";

/// The unit of the branch lengths a tree is drawn in. A function of the focus
/// level so that a future split stays a one-line change; today every level is
/// drawn in SNPs and only the "≈" marker differs. The marker is added by the
/// caller through [`approx_mark`], so this stays a plain unit name.
pub fn unit_label(_focus: u8) -> &'static str {
    TREE_UNIT
}

/// Where the SNP numbers of a tree at this focus level come from, for the
/// legend and the tooltips. Focus 0 trees are sketch distances (substitutions
/// per site) times the median genome length, Focus 1 trees count differing
/// difference pieces (each is at least one difference), Focus 2 trees are SNP
/// counts (6.4, 6.9).
pub fn source_unit(focus: u8) -> &'static str {
    match focus {
        FOCUS_SKETCH => "sketch distance × median genome length",
        FOCUS_BLOCKS => "difference pieces",
        FOCUS_EXACT => "SNP counts",
        _ => "unknown",
    }
}

/// True when a tree at this focus level carries estimated rather than measured
/// distances, and so must be labelled with the "≈" marker.
pub fn is_approximate_focus(focus: u8) -> bool {
    focus < FOCUS_EXACT
}

/// A unit name with its marker: "≈ SNPs" below Focus 2, "SNPs" at Focus 2.
pub fn marked_unit(focus: u8) -> String {
    format!(
        "{}{}",
        approx_mark(is_approximate_focus(focus)),
        unit_label(focus)
    )
}

/// A branch length in SNPs the way a reader wants it: whole numbers with
/// thousands separators, one decimal under ten when neighbour joining split a
/// difference in half, and plain "0" for anything under a hundredth.
pub fn snps(v: f64) -> String {
    if !v.is_finite() || v < 0.0 {
        return "n/a".to_string();
    }
    if v < 0.005 {
        return "0".to_string();
    }
    let r = v.round();
    if (v - r).abs() < 0.05 || v >= 10.0 {
        return count(r as u64);
    }
    format!("{v:.1}")
}

/// An empty allele is written as a dash, as the output files do.
pub fn or_dash(s: &str) -> &str {
    if s.is_empty() {
        "-"
    } else {
        s
    }
}

/// Byte counts, for the memory budget parameter.
pub fn bytes(n: u64) -> String {
    if n == 0 {
        return "0 B".to_string();
    }
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut v = n as f64;
    let mut i = 0;
    while v >= 1024.0 && i + 1 < UNITS.len() {
        v /= 1024.0;
        i += 1;
    }
    if i == 0 {
        format!("{n} B")
    } else {
        format!("{v:.2} {}", UNITS[i])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insertion_positions_are_not_shifted_by_one() {
        assert_eq!(one_based(10, VariantKind::Insertion), 10);
        assert_eq!(one_based(10, VariantKind::Snp), 11);
        assert_eq!(one_based(10, VariantKind::Deletion), 11);
        assert_eq!(one_based(10, VariantKind::Complex), 11);
        assert_eq!(one_based(0, VariantKind::Snp), 1);
    }

    #[test]
    fn ranges_are_inclusive_and_one_based() {
        // [10, 13) on the pivot is bases 11..13 for a deletion.
        assert_eq!(
            one_based_range(10, 13, VariantKind::Deletion),
            "11-13".to_string()
        );
        // A single-base SNP is one position, not a range.
        assert_eq!(one_based_range(10, 11, VariantKind::Snp), "11".to_string());
        // An empty interval (an insertion) does not produce an inverted range.
        assert_eq!(
            one_based_range(10, 10, VariantKind::Insertion),
            "10".to_string()
        );
    }

    #[test]
    fn counters_get_thousands_separators() {
        assert_eq!(count(0), "0");
        assert_eq!(count(999), "999");
        assert_eq!(count(1000), "1\u{2009}000");
        assert_eq!(count(1234567), "1\u{2009}234\u{2009}567");
    }

    #[test]
    fn lengths_pick_their_own_unit() {
        assert_eq!(bp(999), "999 bp");
        assert_eq!(bp(1500), "1.5 kb");
        assert_eq!(bp(5_000_000), "5.00 Mb");
        assert_eq!(bp(3_000_000_000), "3.00 Gb");
    }

    #[test]
    fn durations_stay_readable() {
        assert_eq!(secs(57), "57 ms");
        assert_eq!(secs(2100), "2.1 s");
        assert_eq!(secs(65_000), "1 m 05 s");
    }

    #[test]
    fn an_unmeasured_ratio_is_never_zero() {
        assert_eq!(pct(None), "n/a");
        assert_eq!(pct(Some(0.997)), "99.70 %");
        assert_eq!(pct(Some(f64::NAN)), "n/a");
    }

    #[test]
    fn approximate_numbers_carry_their_marker() {
        assert_eq!(approx_mark(true), "≈ ");
        assert_eq!(approx_mark(false), "");
        assert_eq!(distance(0.005, true), "≈ 5.000e-3");
        assert_eq!(distance(42.0, false), "42");
        assert_eq!(distance(f64::NAN, false), "n/a");
    }

    #[test]
    fn focus_labels_and_units_match_the_vocabulary() {
        assert_eq!(focus_word(0), "blurry");
        assert_eq!(focus_word(1), "medium");
        assert_eq!(focus_word(2), "sharp");
        assert_eq!(focus_label(1, true), "Focus 1≈ (medium)");
        assert_eq!(focus_label(2, false), "Focus 2 (sharp)");
        // Every tree is drawn in SNPs; only the marker and the source differ.
        assert_eq!(unit_label(0), "SNPs");
        assert_eq!(unit_label(1), "SNPs");
        assert_eq!(unit_label(2), "SNPs");
        assert!(source_unit(0).contains("sketch"));
        assert!(source_unit(1).contains("pieces"));
        assert_eq!(source_unit(2), "SNP counts");
    }

    #[test]
    fn snp_lengths_read_as_counts() {
        assert_eq!(snps(0.0), "0");
        assert_eq!(snps(0.001), "0");
        assert_eq!(snps(3.0), "3");
        assert_eq!(snps(2.5), "2.5");
        assert_eq!(snps(12.4), "12");
        assert_eq!(snps(12_345.0), "12\u{2009}345");
        assert_eq!(snps(f64::NAN), "n/a");
        assert_eq!(snps(-1.0), "n/a");
    }

    #[test]
    fn only_an_exact_focus_reads_without_a_marker() {
        assert!(!is_approximate_focus(2));
        assert!(is_approximate_focus(0) && is_approximate_focus(1));
        assert_eq!(marked_unit(2), "SNPs");
        assert_eq!(marked_unit(1), "≈ SNPs");
        assert_eq!(marked_unit(0), "≈ SNPs");
    }

    #[test]
    fn byte_counts_scale() {
        assert_eq!(bytes(0), "0 B");
        assert_eq!(bytes(512), "512 B");
        assert_eq!(bytes(1 << 30), "1.00 GiB");
    }

    #[test]
    fn an_empty_allele_is_a_dash() {
        assert_eq!(or_dash(""), "-");
        assert_eq!(or_dash("ACGT"), "ACGT");
    }
}
