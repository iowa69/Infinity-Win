//! Divergent assemblies: genomes that do not belong with the rest of the input.
//!
//! Infinity-Win assumes its input is one species, assembled reasonably well.
//! A genome of another species, a mislabelled one, a contaminated assembly or
//! a badly incomplete one breaks that assumption quietly: it forms a cluster of
//! its own, stretches the backbone, and can shrink the core genome of a cluster
//! it lands in. This module finds such genomes from what the scan already
//! produced (sketches, lengths, GC), so the run can warn about them before any
//! refinement and offer to leave them out.
//!
//! Three independent signals, each with a documented threshold in
//! [`OutlierThresholds`]:
//!
//! 1. **No close relative.** The sketch distance to the nearest other genome is
//!    above 0.05, below about 95% identity: the conventional species boundary,
//!    so likely another species or a mislabelled genome.
//! 2. **Unusual length.** Outside 0.75-1.30 x the median length: extra sequence
//!    from contamination, or an incomplete assembly.
//! 3. **Unusual GC.** More than 2 percentage points from the median GC: a
//!    mixture with another organism, or another species.
//!
//! Every signal compares a genome with the *rest* of the input, so nothing is
//! flagged with fewer than three genomes.

use rayon::prelude::*;

use crate::events::OutlierInfo;
use crate::genome::GenomeRecord;
use crate::sketch;

/// When a genome counts as divergent. The defaults are deliberately loose: a
/// warning that fires on good genomes would teach users to ignore it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OutlierThresholds {
    /// Below this many genomes there is no "rest" to compare with: nothing is
    /// flagged. Default 3.
    pub min_genomes: usize,
    /// Largest sketch distance (about 1 - identity) to the nearest other
    /// genome. Default 0.05, the ~95% identity species boundary.
    pub max_nearest_dist: f64,
    /// Smallest length as a fraction of the median length. Default 0.75.
    pub min_length_ratio: f64,
    /// Largest length as a fraction of the median length. Default 1.30.
    pub max_length_ratio: f64,
    /// Largest distance of the GC fraction from the median GC fraction.
    /// Default 0.02, two percentage points.
    pub max_gc_diff: f64,
}

impl Default for OutlierThresholds {
    fn default() -> Self {
        Self {
            min_genomes: 3,
            max_nearest_dist: 0.05,
            min_length_ratio: 0.75,
            max_length_ratio: 1.30,
            max_gc_diff: 0.02,
        }
    }
}

/// For every genome, its nearest other genome by sketch distance: `(distance,
/// index)`, or `(1.0, None)` when it is the only one. Computed all-versus-all
/// in parallel; ties go to the lower index, so the answer does not depend on
/// the thread count.
///
/// `same[i] == same[j]` marks two entries as the same genome given twice (the
/// same file name in two input folders): a copy of itself is not a relative,
/// and counting it would hide exactly the genome this is meant to find.
pub fn nearest_neighbours(
    sketches: &[&[u64]],
    same: &[&str],
    k: usize,
) -> Vec<(f64, Option<usize>)> {
    let n = sketches.len();
    (0..n)
        .into_par_iter()
        .map(|i| {
            let mut best = (f64::INFINITY, None);
            for j in 0..n {
                if j == i || same.get(i).is_some_and(|s| same.get(j) == Some(s)) {
                    continue;
                }
                let d = sketch::distance(sketches[i], sketches[j], k);
                if d < best.0 {
                    best = (d, Some(j));
                }
            }
            if best.1.is_none() {
                (1.0, None)
            } else {
                best
            }
        })
        .collect()
}

/// The middle value, averaging the two middle ones of an even count. `None`
/// for an empty slice.
pub fn median(values: &[f64]) -> Option<f64> {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.total_cmp(b));
    let mid = v.len() / 2;
    Some(if v.len() % 2 == 1 {
        v[mid]
    } else {
        (v[mid - 1] + v[mid]) / 2.0
    })
}

/// Find the divergent genomes among `records`, in record order.
///
/// `names` are the names to report (the record names, which are what
/// `--exclude` matches). Must be called inside the run's thread pool to use
/// its thread count.
pub fn find(records: &[GenomeRecord], k: usize, t: &OutlierThresholds) -> Vec<OutlierInfo> {
    let n = records.len();
    if n < t.min_genomes.max(2) {
        return Vec::new();
    }
    let sketches: Vec<&[u64]> = records.iter().map(|r| r.sketch.as_slice()).collect();
    let names: Vec<&str> = records.iter().map(|r| r.name.as_str()).collect();
    let nearest = nearest_neighbours(&sketches, &names, k);
    let lengths: Vec<f64> = records.iter().map(|r| r.stats.length as f64).collect();
    let gcs: Vec<f64> = records.iter().map(|r| r.stats.gc).collect();
    let med_len = median(&lengths).unwrap_or(0.0);
    let med_gc = median(&gcs).unwrap_or(0.0);

    let mut out: Vec<OutlierInfo> = Vec::new();
    for (i, r) in records.iter().enumerate() {
        if out.iter().any(|o| o.name == r.name) {
            // The same genome given twice is reported once: `--exclude` drops
            // every file of that name anyway.
            continue;
        }
        let (nd, nj) = nearest[i];
        // At distance 1 nothing is shared, so "the nearest" is just the first
        // genome looked at and naming it would mislead.
        let nearest_name = nj
            .filter(|_| nd < 1.0)
            .map(|j| records[j].name.clone())
            .unwrap_or_default();
        let mut reasons = Vec::new();
        if nd > t.max_nearest_dist {
            reasons.push(if nd >= 1.0 {
                "no close relative: shares no sketch k-mer with any other genome \
                 (far below 95% identity); likely another species or a mislabelled genome"
                    .to_string()
            } else {
                format!(
                    "no close relative: the nearest genome ({nearest_name}) is {nd:.3} away \
                     (about {:.1}% identity, below the ~95% species boundary); likely another \
                     species or a mislabelled genome",
                    (1.0 - nd) * 100.0
                )
            });
        }
        if med_len > 0.0 {
            let ratio = r.stats.length as f64 / med_len;
            if ratio > t.max_length_ratio {
                reasons.push(format!(
                    "length {:.2} Mb is {ratio:.2}x the median ({:.2} Mb); likely contamination",
                    r.stats.length as f64 / 1e6,
                    med_len / 1e6
                ));
            } else if ratio < t.min_length_ratio {
                reasons.push(format!(
                    "length {:.2} Mb is {ratio:.2}x the median ({:.2} Mb); likely an incomplete \
                     assembly",
                    r.stats.length as f64 / 1e6,
                    med_len / 1e6
                ));
            }
        }
        let gc_diff = r.stats.gc - med_gc;
        if gc_diff.abs() > t.max_gc_diff {
            reasons.push(format!(
                "GC {:.1}% against a median of {:.1}% ({:+.1} points); likely contamination or \
                 another species",
                r.stats.gc * 100.0,
                med_gc * 100.0,
                gc_diff * 100.0
            ));
        }
        if !reasons.is_empty() {
            out.push(OutlierInfo {
                id: r.id,
                name: r.name.clone(),
                length: r.stats.length,
                gc: r.stats.gc,
                nearest_dist: nd,
                nearest_name,
                reasons,
            });
        }
    }
    out
}

/// The note written to `run.json` for a list of flagged genomes: one line per
/// genome, then how to leave them all out.
pub fn notes_for(list: &[OutlierInfo]) -> Vec<String> {
    if list.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<String> = list
        .iter()
        .map(|o| format!("divergent assembly {}: {}", o.name, o.reasons.join("; ")))
        .collect();
    out.push(format!(
        "{} divergent assembl{} may compromise the results; to leave them out run again with {}",
        list.len(),
        if list.len() == 1 { "y" } else { "ies" },
        exclude_option(list.iter().map(|o| o.name.as_str()))
    ));
    out
}

/// The command-line option that excludes these names: `--exclude a,b,c`, with
/// quotes when a name needs them. A name with a comma in it cannot be written
/// that way, so the option then points at `--exclude-file` instead.
pub fn exclude_option<'a>(names: impl IntoIterator<Item = &'a str>) -> String {
    let names: Vec<&str> = names.into_iter().collect();
    if names.iter().any(|n| n.contains(',')) {
        return "--exclude-file <file with one name per line>".to_string();
    }
    let joined = names.join(",");
    let needs_quotes = joined
        .chars()
        .any(|c| c.is_whitespace() || "'\"$&|;<>()*?[]{}!`\\".contains(c));
    if needs_quotes {
        format!("--exclude \"{}\"", joined.replace('"', "\\\""))
    } else {
        format!("--exclude {joined}")
    }
}

/// A list of genome names as typed in a file: one per line, surrounding space
/// trimmed, blank lines and `#` comments ignored.
pub fn parse_name_list(text: &str) -> Vec<String> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

/// Split input files into the ones to scan and the excluded names that matched
/// something, plus the excluded names that matched nothing (worth a note: a
/// typo would otherwise silently keep the genome in).
pub fn apply_exclusions(
    files: Vec<std::path::PathBuf>,
    exclude: &[String],
) -> (Vec<std::path::PathBuf>, Vec<String>, Vec<String>) {
    if exclude.is_empty() {
        return (files, Vec::new(), Vec::new());
    }
    let wanted: std::collections::BTreeSet<&str> = exclude
        .iter()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    let mut kept = Vec::with_capacity(files.len());
    let mut matched: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    let mut dropped = Vec::new();
    for f in files {
        let name = crate::fasta::genome_name(&f);
        if wanted.contains(name.as_str()) {
            matched.insert(name.clone());
            dropped.push(name);
        } else {
            kept.push(f);
        }
    }
    let unmatched = wanted
        .iter()
        .filter(|w| !matched.contains(**w))
        .map(|w| w.to_string())
        .collect();
    (kept, dropped, unmatched)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::{GenomeRecord, GenomeStats};
    use std::path::PathBuf;

    fn rec(id: u32, name: &str, length: u64, gc: f64, sketch: Vec<u64>) -> GenomeRecord {
        let stats = GenomeStats {
            length,
            gc,
            contigs: 1,
            n50: length,
            longest: length,
            ..Default::default()
        };
        GenomeRecord::new(id, name.into(), PathBuf::from(name), stats, sketch)
    }

    /// A sketch sharing `shared` of its 1000 hashes with the base sketch.
    fn sketch_sharing(shared: u64, salt: u64) -> Vec<u64> {
        let mut v: Vec<u64> = (0..shared)
            .chain((0..1000 - shared).map(|i| 1_000_000 * (salt + 1) + i))
            .collect();
        v.sort_unstable();
        v
    }

    fn base(n: u32) -> Vec<GenomeRecord> {
        (0..n)
            .map(|i| {
                rec(
                    i,
                    &format!("g{i}"),
                    2_800_000 + u64::from(i) * 1000,
                    0.328,
                    sketch_sharing(990, u64::from(i)),
                )
            })
            .collect()
    }

    #[test]
    fn a_homogeneous_set_flags_nothing() {
        let recs = base(6);
        assert!(find(&recs, 31, &OutlierThresholds::default()).is_empty());
    }

    #[test]
    fn another_species_is_flagged_for_its_distance_and_gc() {
        let mut recs = base(6);
        // Shares nothing, different length and GC: a different species.
        recs.push(rec(
            6,
            "odd",
            2_000_000,
            0.65,
            (900_000_000..900_001_000).collect(),
        ));
        let out = find(&recs, 31, &OutlierThresholds::default());
        assert_eq!(out.len(), 1, "{out:?}");
        let o = &out[0];
        assert_eq!(o.name, "odd");
        assert_eq!(o.nearest_dist, 1.0);
        assert_eq!(o.reasons.len(), 3, "{:?}", o.reasons);
        assert!(o.reasons[0].contains("no close relative"));
        assert!(o.reasons[1].contains("incomplete"));
        assert!(o.reasons[2].contains("GC 65.0%"));
    }

    #[test]
    fn contamination_shows_as_length_and_gc_even_with_a_close_relative() {
        let mut recs = base(6);
        let mut big = rec(6, "contaminated", 3_900_000, 0.36, sketch_sharing(990, 99));
        big.id = 6;
        recs.push(big);
        let out = find(&recs, 31, &OutlierThresholds::default());
        assert_eq!(out.len(), 1);
        assert!(out[0].reasons.iter().any(|r| r.contains("contamination")));
        assert!(out[0].nearest_dist < 0.05, "its sketch still matches");
        assert!(!out[0].nearest_name.is_empty());
    }

    #[test]
    fn two_genomes_are_never_compared_against_a_rest() {
        let mut recs = base(1);
        recs.push(rec(1, "odd", 1_000_000, 0.65, vec![1, 2, 3]));
        assert!(find(&recs, 31, &OutlierThresholds::default()).is_empty());
    }

    #[test]
    fn nearest_neighbours_break_ties_by_index_and_handle_a_single_genome() {
        let a = vec![1u64, 2, 3];
        let b = vec![1u64, 2, 3];
        let c = vec![1u64, 2, 3];
        let nn = nearest_neighbours(&[&a, &b, &c], &["a", "b", "c"], 31);
        assert_eq!(nn[0], (0.0, Some(1)));
        assert_eq!(nn[1], (0.0, Some(0)));
        assert_eq!(nn[2], (0.0, Some(0)));
        assert_eq!(nearest_neighbours(&[&a], &["a"], 31), vec![(1.0, None)]);
        // A copy of itself under the same name is not a relative.
        let nn = nearest_neighbours(&[&a, &b, &c], &["a", "a", "c"], 31);
        assert_eq!(nn[0], (0.0, Some(2)));
    }

    #[test]
    fn a_genome_given_twice_is_still_found_and_reported_once() {
        let mut recs = base(5);
        let odd: Vec<u64> = (900_000_000..900_001_000).collect();
        recs.push(rec(5, "odd", 2_800_000, 0.328, odd.clone()));
        recs.push(rec(6, "odd", 2_800_000, 0.328, odd));
        let out = find(&recs, 31, &OutlierThresholds::default());
        assert_eq!(out.len(), 1, "{out:?}");
        assert_eq!(out[0].nearest_dist, 1.0);
        assert!(out[0].nearest_name.is_empty(), "nothing shared, no name");
    }

    #[test]
    fn medians_average_the_middle_pair() {
        assert_eq!(median(&[]), None);
        assert_eq!(median(&[3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(median(&[4.0, 1.0, 2.0, 3.0]), Some(2.5));
    }

    #[test]
    fn the_exclude_option_is_ready_to_paste() {
        assert_eq!(exclude_option(["a", "b"]), "--exclude a,b");
        assert_eq!(exclude_option(["my genome"]), "--exclude \"my genome\"");
        assert!(exclude_option(["a,b"]).starts_with("--exclude-file"));
    }

    #[test]
    fn name_lists_skip_blanks_and_comments() {
        let names = parse_name_list("a\n\n  b  \n# a comment\nc\r\n");
        assert_eq!(names, vec!["a", "b", "c"]);
    }

    #[test]
    fn exclusions_drop_files_by_name_and_report_typos() {
        let files = vec![
            PathBuf::from("/x/a.fasta"),
            PathBuf::from("/x/b.fna.gz"),
            PathBuf::from("/x/c.fa"),
        ];
        let (kept, dropped, unmatched) =
            apply_exclusions(files.clone(), &["b".into(), "zz".into(), " ".into()]);
        assert_eq!(kept, vec![files[0].clone(), files[2].clone()]);
        assert_eq!(dropped, vec!["b".to_string()]);
        assert_eq!(unmatched, vec!["zz".to_string()]);
        let (kept, dropped, unmatched) = apply_exclusions(files.clone(), &[]);
        assert_eq!(kept, files);
        assert!(dropped.is_empty() && unmatched.is_empty());
    }

    #[test]
    fn notes_end_with_the_way_to_exclude() {
        let o = OutlierInfo {
            id: 0,
            name: "x".into(),
            length: 1,
            gc: 0.5,
            nearest_dist: 1.0,
            nearest_name: String::new(),
            reasons: vec!["r".into()],
        };
        let n = notes_for(&[o]);
        assert_eq!(n.len(), 2);
        assert_eq!(n[0], "divergent assembly x: r");
        assert!(n[1].ends_with("--exclude x"), "{}", n[1]);
        assert!(notes_for(&[]).is_empty());
    }
}
