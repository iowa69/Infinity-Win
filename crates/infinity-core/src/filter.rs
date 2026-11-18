//! SNP rules: which differences count towards a distance.
//!
//! The comparison engine reports every difference between a member and its
//! pivot. Not every reported difference is a mutation: collapsed repeats,
//! contig ends, short high-copy contigs, recombination and misassemblies all
//! produce runs of "SNPs" that a core-genome pipeline (SKA2, Snippy-core,
//! Parsnp, Gubbins) would never count. This module applies the same kind of
//! rules, using the one thing this engine has that a k-mer sketch does not:
//! the exact pivot position of every difference and of every uncovered base.
//!
//! A rule never just deletes a SNP. It **masks** a stretch of the pivot for
//! one member (or for the whole cluster), so the masked bases stop counting in
//! the comparable denominator as well as in the difference count. Distances
//! therefore stay "SNPs per comparable site" with an honest denominator.
//!
//! Rules, in the order they are applied:
//!
//! | rule | scope | masks |
//! |---|---|---|
//! | `short_contig` | cluster | pivot contigs shorter than `min_contig_len` |
//! | `contig_end` | cluster | `end_margin` bases at both ends of every pivot contig |
//! | `high_copy` | cluster | pivot contigs whose coverage (from the assembler's contig name) exceeds `high_copy_ratio` x the genome's typical coverage |
//! | `not_core` | cluster | every pivot base that any member does not cover |
//! | `flag:<names>` | member | variants the engine itself flagged (repeat, end extension, low identity, ...) |
//! | `snp_cluster` | member | windows holding more SNPs than a Poisson process allows (recombination, misassembly) |
//!
//! The SNP-cluster test: under clonal descent mutations fall independently along
//! the genome, so the number of SNPs a member carries in a window of `w` bases
//! is Poisson with mean `lambda = w * n / L` (n SNPs over L core bases). A
//! window holding `k` SNPs is rejected when `P(X >= k) * (L / w) < alpha`, a
//! Bonferroni correction over the `L / w` windows of the genome. This is the
//! idea behind the spatial scan statistic of Gubbins, applied per member.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::interval::{contains, merge, total_len, Iv};
use crate::profile::Profile;
use crate::variant::{flags, VariantKind};

/// User-facing settings of the SNP rules.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SnpFilter {
    /// Master switch. Off reproduces the unfiltered whole-assembly counts.
    pub enabled: bool,
    /// Count only pivot bases in the cluster's core genome.
    pub core: bool,
    /// A pivot base is core when at least this fraction of the cluster's
    /// genomes cover it (1.0 = every genome). Below 1.0, a base that one
    /// fragmented assembly misses at a contig break stays core for everyone
    /// else, as in SKA2's `--min-freq`.
    pub core_min_fraction: f64,
    /// Pivot contigs shorter than this are masked.
    pub min_contig_len: u64,
    /// Bases masked at each end of every pivot contig.
    pub end_margin: u64,
    /// Mask a pivot contig whose coverage exceeds this multiple of the typical
    /// coverage. Needs coverage in the contig names (SPAdes, Shovill, SKESA,
    /// Unicycler style); `0` turns the rule off.
    pub high_copy_ratio: f64,
    /// Exclude differences the engine flagged as uncertain.
    pub drop_flagged: bool,
    /// Count SNPs found by extending across a breakpoint (flag `break_ext`)
    /// even though they carry a flag. Off by default: on real assemblies these
    /// sit next to contig junctions and produced a false SNP in the Salmonella
    /// acceptance trio, while the coverage the extension adds is kept either way.
    pub keep_break_ext: bool,
    /// Window of the SNP-cluster scan, in bases; `0` turns the scan off.
    pub scan_window: u64,
    /// Family-wise false-rejection rate of the SNP-cluster scan.
    pub scan_alpha: f64,
    /// Window of the four-gamete test, in bases; `0` turns it off.
    pub incompat_window: u64,
    /// A site is rejected when it is incompatible with at least this many
    /// other sites inside the window.
    pub incompat_min: u32,
}

impl Default for SnpFilter {
    fn default() -> Self {
        Self {
            enabled: true,
            core: true,
            core_min_fraction: 1.0,
            min_contig_len: 5_000,
            end_margin: 150,
            high_copy_ratio: 1.8,
            drop_flagged: true,
            keep_break_ext: false,
            scan_window: 1_000,
            scan_alpha: 0.01,
            incompat_window: 1_000,
            incompat_min: 2,
        }
    }
}

impl SnpFilter {
    /// The rules switched off: every reported difference counts.
    pub fn off() -> Self {
        Self {
            enabled: false,
            ..Self::default()
        }
    }

    pub fn validate(&self) -> Result<(), String> {
        if !(self.high_copy_ratio == 0.0 || self.high_copy_ratio > 1.0) {
            return Err(format!(
                "high-copy ratio must be 0 (off) or greater than 1, got {}",
                self.high_copy_ratio
            ));
        }
        if !(self.core_min_fraction > 0.0 && self.core_min_fraction <= 1.0) {
            return Err(format!(
                "core fraction must be in (0, 1], got {}",
                self.core_min_fraction
            ));
        }
        if !(self.scan_alpha > 0.0 && self.scan_alpha < 1.0) {
            return Err(format!(
                "SNP-cluster alpha must be between 0 and 1, got {}",
                self.scan_alpha
            ));
        }
        Ok(())
    }
}

/// Why a variant does not count. `None` in [`FilterOutcome::reasons`] is PASS.
pub const PASS: &str = "PASS";

/// The result of applying the rules to one cluster.
#[derive(Debug, Clone)]
pub struct FilterOutcome {
    /// One filtered profile per input profile, in the same order: rejected
    /// SNPs and indels removed, masked stretches added to `uncovered`, and only
    /// passing variants kept in `variants`.
    pub profiles: Vec<Profile>,
    /// Per input profile, per variant (in the profile's `variants` order):
    /// `PASS` or the `;`-separated reasons it does not count.
    pub reasons: Vec<Vec<String>>,
    /// Pivot bases left after the cluster-wide masks: the core genome.
    pub core_bp: u64,
    /// Bases removed by each cluster-wide rule (overlaps counted once, in rule
    /// order), for `run.json` and the report.
    pub masked_by_rule: Vec<(String, u64)>,
    /// Passing SNPs per input profile.
    pub pass_snps: Vec<u64>,
}

/// Coverage written into a contig name by common assemblers:
/// `NODE_7_length_774551_cov_139.9` (SPAdes/Shovill), `..._cov42.0`,
/// `Contig_1_42.5` (SKESA writes it last), `... depth=1.23x` (Unicycler).
pub fn coverage_from_name(name: &str) -> Option<f64> {
    let lower = name.to_ascii_lowercase();
    for key in ["_cov_", "_cov", "depth=", "cov="] {
        if let Some(i) = lower.find(key) {
            let rest = &lower[i + key.len()..];
            let num: String = rest
                .chars()
                .take_while(|c| c.is_ascii_digit() || *c == '.')
                .collect();
            if let Ok(v) = num.parse::<f64>() {
                if v.is_finite() && v > 0.0 {
                    return Some(v);
                }
            }
        }
    }
    None
}

/// Typical coverage of an assembly: the length-weighted median over contigs of
/// at least 50 kb (the chromosome), falling back to all contigs with coverage.
fn typical_coverage(names: &[String], lens: &[u64]) -> Option<f64> {
    let pick = |min_len: u64| -> Vec<(f64, u64)> {
        names
            .iter()
            .zip(lens)
            .filter(|(_, &l)| l >= min_len)
            .filter_map(|(n, &l)| coverage_from_name(n).map(|c| (c, l)))
            .collect()
    };
    let mut v = pick(50_000);
    if v.is_empty() {
        v = pick(0);
    }
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: u64 = v.iter().map(|x| x.1).sum();
    let mut acc = 0u64;
    for (c, l) in &v {
        acc += l;
        if acc * 2 >= total {
            return Some(*c);
        }
    }
    v.last().map(|x| x.0)
}

/// `P(X >= k)` for `X ~ Poisson(lambda)`, summed from the tail so that tiny
/// probabilities do not cancel to zero.
pub fn poisson_tail(k: u64, lambda: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if lambda <= 0.0 {
        return 0.0;
    }
    // log of the first term, e^-l l^k / k!
    let mut log_term = -lambda + (k as f64) * lambda.ln() - ln_factorial(k);
    let mut sum = 0.0f64;
    let mut i = k;
    loop {
        let t = log_term.exp();
        sum += t;
        i += 1;
        log_term += lambda.ln() - (i as f64).ln();
        if t < sum * 1e-15 || i > k + 10_000 {
            break;
        }
    }
    sum.min(1.0)
}

fn ln_factorial(n: u64) -> f64 {
    (1..=n).map(|i| (i as f64).ln()).sum()
}

/// Positions (flat pivot coordinates, sorted) of the SNPs that form
/// significant clusters, given the member's candidate SNPs.
fn snp_cluster_windows(snps: &[u64], core_bp: u64, window: u64, alpha: f64) -> Vec<Iv> {
    if window == 0 || snps.len() < 2 || core_bp == 0 {
        return Vec::new();
    }
    let n = snps.len() as f64;
    let lambda = n * window as f64 / core_bp as f64;
    let windows = (core_bp as f64 / window as f64).max(1.0);
    let half = window / 2;
    let mut flagged = vec![false; snps.len()];
    let (mut lo, mut hi) = (0usize, 0usize);
    for (i, &p) in snps.iter().enumerate() {
        while snps[lo] + half < p {
            lo += 1;
        }
        while hi + 1 < snps.len() && snps[hi + 1] <= p + half {
            hi += 1;
        }
        let k = (hi - lo + 1) as u64;
        if k >= 2 && poisson_tail(k, lambda) * windows < alpha {
            for f in flagged.iter_mut().take(hi + 1).skip(lo) {
                *f = true;
            }
        }
        let _ = i;
    }
    // Mask from the first to the last flagged SNP of each run (runs are broken
    // by a gap wider than half a window).
    let mut out = Vec::new();
    let mut run: Option<(u64, u64)> = None;
    for (i, &p) in snps.iter().enumerate() {
        if !flagged[i] {
            continue;
        }
        run = match run {
            Some((s, e)) if p <= e + half => Some((s, p + 1)),
            Some(r) => {
                out.push(r);
                Some((p, p + 1))
            }
            None => Some((p, p + 1)),
        };
    }
    if let Some(r) = run {
        out.push(r);
    }
    out
}

/// Pivot intervals that more than `(1 - fraction) * n` of the `n` genomes do
/// not cover. With `fraction == 1.0` this is the union of all uncovered lists.
pub fn soft_core_mask(profiles: &[Profile], fraction: f64) -> Vec<Iv> {
    let n = profiles.len();
    let allowed = (((1.0 - fraction) * n as f64) + 1e-9).floor() as i64;
    let mut ev: Vec<(u64, i64)> = Vec::new();
    for p in profiles {
        for &(s, e) in &p.uncovered {
            if e > s {
                ev.push((s, 1));
                ev.push((e, -1));
            }
        }
    }
    // ends before starts at the same coordinate: half-open intervals
    ev.sort_unstable_by_key(|&(x, d)| (x, d));
    let mut out = Vec::new();
    let mut depth = 0i64;
    let mut open: Option<u64> = None;
    for (x, d) in ev {
        depth += d;
        if depth > allowed && open.is_none() {
            open = Some(x);
        } else if depth <= allowed {
            if let Some(s) = open.take() {
                if x > s {
                    out.push((s, x));
                }
            }
        }
    }
    merge(out)
}

/// Sites that fail the four-gamete test against at least `min_conflicts`
/// other sites within `window` bases.
///
/// Two substitutions that each happened once on a clonal tree split the genomes
/// in a nested or disjoint way, so at most three of the four combinations of
/// their alleles occur (`00` is always present: the pivot). When all four occur,
/// at least one of the two sites is not a single clonal mutation: recombination,
/// or a collapsed repeat whose copies the assemblies resolved differently. Real
/// outbreak SNPs essentially never conflict with their neighbours; collapsed
/// repeats conflict with all of them. `survivors[i]` is genome `i`'s sorted
/// (position, base) list after the other rules.
pub fn incompatible_sites(
    survivors: &[Vec<(u64, u8)>],
    window: u64,
    min_conflicts: u32,
) -> Vec<u64> {
    use std::collections::BTreeMap;
    let n = survivors.len();
    let words = n.div_ceil(64).max(1);
    // position -> (base, carriers bitset, multi-allelic)
    let mut sites: BTreeMap<u64, (u8, Vec<u64>, bool)> = BTreeMap::new();
    for (gi, list) in survivors.iter().enumerate() {
        for &(p, b) in list {
            let e = sites
                .entry(p)
                .or_insert_with(|| (b, vec![0u64; words], false));
            if e.0 != b {
                e.2 = true;
            }
            e.1[gi / 64] |= 1u64 << (gi % 64);
        }
    }
    let pos: Vec<u64> = sites.keys().copied().collect();
    let info: Vec<&(u8, Vec<u64>, bool)> = sites.values().collect();
    let mut valid = vec![!0u64; words];
    if n % 64 != 0 {
        valid[words - 1] = (1u64 << (n % 64)) - 1;
    }
    let mut conflicts = vec![0u32; pos.len()];
    for i in 0..pos.len() {
        if info[i].2 {
            continue;
        }
        let mut j = i + 1;
        while j < pos.len() && pos[j] - pos[i] <= window {
            if !info[j].2 {
                let (a, b) = (&info[i].1, &info[j].1);
                let (mut h11, mut h10, mut h01, mut h00) = (false, false, false, false);
                for w in 0..words {
                    h11 |= a[w] & b[w] != 0;
                    h10 |= a[w] & !b[w] & valid[w] != 0;
                    h01 |= !a[w] & b[w] & valid[w] != 0;
                    h00 |= !a[w] & !b[w] & valid[w] != 0;
                }
                if h11 && h10 && h01 && h00 {
                    conflicts[i] += 1;
                    conflicts[j] += 1;
                }
            }
            j += 1;
        }
    }
    pos.iter()
        .zip(&conflicts)
        .filter(|(_, &c)| c >= min_conflicts)
        .map(|(&p, _)| p)
        .collect()
}

/// Apply the SNP rules to the profiles of one cluster.
///
/// `profiles[0]` is conventionally the pivot (an empty profile), but nothing
/// here depends on that. Coordinates are flat pivot coordinates; `offsets[c]`
/// is where pivot contig `c` starts.
pub fn apply(
    profiles: &[Profile],
    contig_names: &[String],
    contig_lens: &[u64],
    offsets: &[u64],
    pivot_len: u64,
    rules: &SnpFilter,
) -> FilterOutcome {
    if !rules.enabled {
        return FilterOutcome {
            profiles: profiles.to_vec(),
            reasons: profiles
                .iter()
                .map(|p| vec![PASS.to_string(); p.variants.len()])
                .collect(),
            core_bp: pivot_len,
            masked_by_rule: Vec::new(),
            pass_snps: profiles.iter().map(|p| p.snp_map.len() as u64).collect(),
        };
    }

    // ---- cluster-wide masks, one list per rule ----
    let mut short = Vec::new();
    let mut ends = Vec::new();
    let mut high = Vec::new();
    let typical = if rules.high_copy_ratio > 0.0 {
        typical_coverage(contig_names, contig_lens)
    } else {
        None
    };
    for (c, &len) in contig_lens.iter().enumerate() {
        let off = offsets.get(c).copied().unwrap_or(0);
        if len < rules.min_contig_len {
            short.push((off, off + len));
            continue;
        }
        if let (Some(t), Some(cov)) = (
            typical,
            contig_names.get(c).and_then(|n| coverage_from_name(n)),
        ) {
            if cov > rules.high_copy_ratio * t {
                high.push((off, off + len));
                continue;
            }
        }
        if rules.end_margin > 0 {
            let m = rules.end_margin.min(len);
            ends.push((off, off + m));
            ends.push((off + len - m, off + len));
        }
    }
    let short = merge(short);
    let ends = merge(ends);
    let high = merge(high);
    let not_core = if rules.core {
        soft_core_mask(profiles, rules.core_min_fraction)
    } else {
        Vec::new()
    };
    let rule_lists: [(&str, &Vec<Iv>); 4] = [
        ("short_contig", &short),
        ("contig_end", &ends),
        ("high_copy", &high),
        ("not_core", &not_core),
    ];
    let mut cluster_mask: Vec<Iv> = Vec::new();
    let mut masked_by_rule = Vec::new();
    for (name, list) in rule_lists.iter() {
        let before = total_len(&cluster_mask);
        cluster_mask = merge(cluster_mask.iter().chain(list.iter()).copied().collect());
        masked_by_rule.push((name.to_string(), total_len(&cluster_mask) - before));
    }
    let core_bp = pivot_len.saturating_sub(total_len(&cluster_mask));

    // ---- per-member masks (pass 1) ----
    let flat = |c: u32, p: u64| offsets.get(c as usize).copied().unwrap_or(0) + p;
    let counted_flags = |f: u32| -> u32 {
        if rules.keep_break_ext {
            f & !flags::BREAK_EXT
        } else {
            f
        }
    };
    let mut member_masks: Vec<(Vec<Iv>, Vec<Iv>)> = Vec::with_capacity(profiles.len());
    let mut survivors: Vec<Vec<(u64, u8)>> = Vec::with_capacity(profiles.len());
    for prof in profiles {
        let mut flagged_mask: Vec<Iv> = Vec::new();
        if rules.drop_flagged {
            for v in &prof.variants {
                if counted_flags(v.flags) != 0 {
                    let s = flat(v.pivot_contig, v.pivot_start);
                    flagged_mask.push((s, flat(v.pivot_contig, v.pivot_end).max(s + 1)));
                }
            }
        }
        let flagged_mask = merge(flagged_mask);
        let mut candidates: Vec<u64> = prof
            .snp_map
            .keys()
            .copied()
            .filter(|&p| !contains(&cluster_mask, p) && !contains(&flagged_mask, p))
            .collect();
        candidates.sort_unstable();
        // Two scales: the configured window, and a window ten times finer that
        // catches a tight run of three or more SNPs (a collapsed repeat or a
        // misassembly) which a 1 kb window dilutes. Each scale is Bonferroni
        // corrected over its own number of windows, so a real double mutation
        // (two adjacent SNPs) is still kept.
        let mut scan: Vec<Iv> =
            snp_cluster_windows(&candidates, core_bp, rules.scan_window, rules.scan_alpha);
        if rules.scan_window >= 500 {
            scan.extend(snp_cluster_windows(
                &candidates,
                core_bp,
                rules.scan_window / 10,
                rules.scan_alpha,
            ));
        }
        let scan_mask = merge(scan);
        let member_mask = merge(
            flagged_mask
                .iter()
                .chain(scan_mask.iter())
                .copied()
                .collect(),
        );
        let mut alive: Vec<(u64, u8)> = candidates
            .iter()
            .filter(|&&p| !contains(&member_mask, p))
            .map(|&p| (p, prof.snp_map[&p]))
            .collect();
        alive.sort_unstable();
        survivors.push(alive);
        member_masks.push((scan_mask, member_mask));
    }

    // ---- cluster-wide: the four-gamete test ----
    let incompatible: Vec<Iv> = if rules.incompat_window > 0 {
        incompatible_sites(&survivors, rules.incompat_window, rules.incompat_min)
            .into_iter()
            .map(|p| (p, p + 1))
            .collect()
    } else {
        Vec::new()
    };

    // ---- outputs (pass 2) ----
    let mut out_profiles = Vec::with_capacity(profiles.len());
    let mut reasons = Vec::with_capacity(profiles.len());
    let mut pass_snps = Vec::with_capacity(profiles.len());
    for (prof, (scan_mask, member_mask)) in profiles.iter().zip(member_masks) {
        let member_mask = merge(
            member_mask
                .iter()
                .chain(incompatible.iter())
                .copied()
                .collect(),
        );
        let masked = |s: u64, e: u64| -> bool {
            let e = e.max(s + 1);
            (s..e).any(|p| contains(&cluster_mask, p) || contains(&member_mask, p))
        };

        // reasons, per variant
        let mut why = Vec::with_capacity(prof.variants.len());
        let mut keep = Vec::with_capacity(prof.variants.len());
        for v in &prof.variants {
            let s = flat(v.pivot_contig, v.pivot_start);
            let e = flat(v.pivot_contig, v.pivot_end).max(s + 1);
            let mut r: Vec<String> = Vec::new();
            let hit = |list: &[Iv]| (s..e).any(|p| contains(list, p));
            for (name, list) in rule_lists.iter() {
                if hit(list) {
                    r.push(name.to_string());
                }
            }
            if rules.drop_flagged && counted_flags(v.flags) != 0 {
                r.push(format!(
                    "flag:{}",
                    flags::names(counted_flags(v.flags)).join(",")
                ));
            }
            if v.kind == VariantKind::Snp && hit(&scan_mask) {
                r.push("snp_cluster".to_string());
            }
            if v.kind == VariantKind::Snp && hit(&incompatible) {
                r.push("incompatible".to_string());
            }
            keep.push(r.is_empty());
            why.push(if r.is_empty() {
                PASS.to_string()
            } else {
                r.join(";")
            });
        }

        let mut f = prof.clone();
        f.snp_map
            .retain(|&p, _| !contains(&cluster_mask, p) && !contains(&member_mask, p));
        f.indels.retain(|k| !masked(k.0, k.1));
        f.variants = prof
            .variants
            .iter()
            .zip(&keep)
            .filter(|(_, &k)| k)
            .map(|(v, _)| v.clone())
            .collect();
        f.uncovered = merge(
            prof.uncovered
                .iter()
                .chain(cluster_mask.iter())
                .chain(member_mask.iter())
                .copied()
                .collect(),
        );
        pass_snps.push(f.snp_map.len() as u64);
        out_profiles.push(f);
        reasons.push(why);
    }

    FilterOutcome {
        profiles: out_profiles,
        reasons,
        core_bp,
        masked_by_rule,
        pass_snps,
    }
}

/// Number of SNP clusters at each threshold: single-linkage groups of genomes
/// whose pairwise SNP distance is at most `t`. Returns, per threshold, a group
/// index per genome (groups numbered from 0 in order of first appearance).
pub fn snp_groups(snps: &[Vec<u64>], thresholds: &[u64]) -> Vec<Vec<u32>> {
    let n = snps.len();
    thresholds
        .iter()
        .map(|&t| {
            let mut parent: Vec<usize> = (0..n).collect();
            fn find(p: &mut [usize], x: usize) -> usize {
                let mut r = x;
                while p[r] != r {
                    r = p[r];
                }
                let mut y = x;
                while p[y] != r {
                    let nx = p[y];
                    p[y] = r;
                    y = nx;
                }
                r
            }
            for i in 0..n {
                for j in (i + 1)..n {
                    if snps[i][j] <= t {
                        let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                        if a != b {
                            parent[a.max(b)] = a.min(b);
                        }
                    }
                }
            }
            let mut label: HashMap<usize, u32> = HashMap::new();
            (0..n)
                .map(|i| {
                    let r = find(&mut parent, i);
                    let next = label.len() as u32;
                    *label.entry(r).or_insert(next)
                })
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::Variant;

    fn snp(contig: u32, pos: u64) -> Variant {
        Variant {
            kind: VariantKind::Snp,
            pivot_contig: contig,
            pivot_start: pos,
            pivot_end: pos + 1,
            ref_allele: "A".into(),
            alt_allele: "C".into(),
            genome_contig: 0,
            genome_start: pos,
            genome_end: pos + 1,
            strand: '+',
            flags: 0,
        }
    }

    fn profile(id: u32, vars: Vec<Variant>, offsets: &[u64], uncovered: Vec<Iv>) -> Profile {
        let mut p = Profile::pivot(id);
        for v in &vars {
            p.snp_map
                .insert(offsets[v.pivot_contig as usize] + v.pivot_start, b'C');
        }
        p.variants = vars;
        p.uncovered = merge(uncovered);
        p
    }

    #[test]
    fn coverage_is_read_from_common_assembler_names() {
        assert_eq!(
            coverage_from_name("NODE_7_length_774551_cov_139.9142_unk"),
            Some(139.9142)
        );
        assert_eq!(
            coverage_from_name("isolate-B_contig37_len1448_cov139.5"),
            Some(139.5)
        );
        assert_eq!(
            coverage_from_name("1 length=4800000 depth=1.00x"),
            Some(1.0)
        );
        assert_eq!(coverage_from_name("chromosome"), None);
    }

    #[test]
    fn poisson_tail_matches_closed_forms() {
        let l: f64 = 0.5;
        assert!((poisson_tail(1, l) - (1.0 - (-l).exp())).abs() < 1e-12);
        let p2 = 1.0 - (-l).exp() * (1.0 + l);
        assert!((poisson_tail(2, l) - p2).abs() < 1e-12);
        // tiny lambda: no cancellation to zero
        let t = poisson_tail(3, 1e-4);
        assert!(t > 1.6e-13 && t < 1.7e-13, "{t}");
        assert_eq!(poisson_tail(0, 3.0), 1.0);
    }

    #[test]
    fn short_and_high_copy_contigs_and_ends_are_masked_for_everyone() {
        // contig 0: 100 kb chromosome at 40x, contig 1: 1.4 kb at 140x,
        // contig 2: 20 kb plasmid at 160x
        let names = vec![
            "c1_len100000_cov40.0".to_string(),
            "c2_len1400_cov140.0".to_string(),
            "c3_len20000_cov160.0".to_string(),
        ];
        let lens = vec![100_000, 1_400, 20_000];
        let offsets = vec![0, 100_000, 101_400];
        let total = 121_400;
        let m = profile(
            1,
            vec![snp(0, 50), snp(0, 50_000), snp(1, 700), snp(2, 5_000)],
            &offsets,
            vec![],
        );
        let piv = Profile::pivot(0);
        let out = apply(
            &[piv, m],
            &names,
            &lens,
            &offsets,
            total,
            &SnpFilter::default(),
        );
        assert_eq!(out.reasons[1][0], "contig_end");
        assert_eq!(out.reasons[1][1], PASS);
        assert_eq!(out.reasons[1][2], "short_contig");
        assert_eq!(out.reasons[1][3], "high_copy");
        assert_eq!(out.pass_snps[1], 1);
        assert_eq!(out.core_bp, 100_000 - 300);
        // the pivot gets the same cluster-wide masks
        assert_eq!(total_len(&out.profiles[0].uncovered), total - out.core_bp);
    }

    #[test]
    fn a_base_one_member_lacks_is_not_core_for_anyone() {
        let names = vec!["chr".to_string()];
        let lens = vec![1_000_000];
        let offsets = vec![0];
        let a = profile(1, vec![snp(0, 10_000)], &offsets, vec![]);
        let b = profile(2, vec![], &offsets, vec![(9_000, 11_000)]);
        let out = apply(
            &[Profile::pivot(0), a, b],
            &names,
            &lens,
            &offsets,
            1_000_000,
            &SnpFilter::default(),
        );
        assert_eq!(out.reasons[1][0], "not_core");
        assert_eq!(out.pass_snps[1], 0);
    }

    #[test]
    fn a_dense_snp_window_is_rejected_and_isolated_snps_are_kept() {
        let names = vec!["chr".to_string()];
        let lens = vec![2_000_000];
        let offsets = vec![0];
        let mut vars: Vec<Variant> = (0..10).map(|i| snp(0, 200_000 * (i + 1))).collect();
        for p in [
            1_000_000 + 7,
            1_000_000 + 19,
            1_000_000 + 40,
            1_000_000 + 300,
        ] {
            vars.push(snp(0, p));
        }
        let m = profile(1, vars, &offsets, vec![]);
        let out = apply(
            &[Profile::pivot(0), m],
            &names,
            &lens,
            &offsets,
            2_000_000,
            &SnpFilter::default(),
        );
        let rejected = out.reasons[1]
            .iter()
            .filter(|r| r.contains("snp_cluster"))
            .count();
        // the four clustered SNPs plus the isolated SNP at exactly 1,000,000
        // that sits inside the same window
        assert!(rejected >= 4, "{:?}", out.reasons[1]);
        assert!(out.pass_snps[1] >= 9);
    }

    #[test]
    fn a_soft_core_tolerates_a_few_genomes_missing_a_base() {
        let names = vec!["chr".to_string()];
        let lens = vec![1_000_000];
        let offsets = vec![0];
        let mut profs = vec![Profile::pivot(0)];
        for i in 1..20u32 {
            // genome 1 misses 9-11 kb; everyone else covers it
            let unc = if i == 1 {
                vec![(9_000, 11_000)]
            } else {
                vec![]
            };
            profs.push(profile(i, vec![snp(0, 10_000)], &offsets, unc));
        }
        let strict = SnpFilter {
            core_min_fraction: 1.0,
            ..SnpFilter::default()
        };
        let out = apply(&profs, &names, &lens, &offsets, 1_000_000, &strict);
        assert_eq!(
            out.pass_snps[2], 0,
            "strict core drops the SNP for everyone"
        );
        let soft = SnpFilter {
            core_min_fraction: 0.9,
            ..SnpFilter::default()
        };
        let out = apply(&profs, &names, &lens, &offsets, 1_000_000, &soft);
        assert_eq!(out.pass_snps[2], 1, "the soft core keeps it");
        assert_eq!(out.reasons[1][0], PASS, "genome 1 has its own call there");
        assert_eq!(soft_core_mask(&profs, 0.9), vec![]);
        assert_eq!(soft_core_mask(&profs, 1.0), vec![(9_000, 11_000)]);
    }

    #[test]
    fn the_four_gamete_test_rejects_scrambled_sites_and_keeps_clonal_ones() {
        // genomes 0 (pivot) .. 5. Clonal: site 100 in {1,2}, site 200 in {1}.
        // Scrambled region: sites 5000.. carried by unrelated subsets.
        let g = |sites: &[(u64, u8)]| sites.to_vec();
        let survivors = vec![
            g(&[]),
            g(&[(100, b'A'), (200, b'C'), (5000, b'T'), (5100, b'T')]),
            g(&[(100, b'A'), (5100, b'T'), (5200, b'G')]),
            g(&[(5000, b'T'), (5200, b'G')]),
            g(&[(5000, b'T'), (5100, b'T'), (5200, b'G')]),
            g(&[]),
        ];
        let bad = incompatible_sites(&survivors, 1_000, 2);
        assert_eq!(bad, vec![5000, 5100, 5200]);
        assert!(incompatible_sites(&survivors, 1_000, 3).is_empty());
    }

    #[test]
    fn switched_off_nothing_changes() {
        let names = vec!["chr".to_string()];
        let m = profile(1, vec![snp(0, 10)], &[0], vec![]);
        let out = apply(&[m.clone()], &names, &[100], &[0], 100, &SnpFilter::off());
        assert_eq!(out.profiles[0], m);
        assert_eq!(out.reasons[0], vec![PASS.to_string()]);
    }

    #[test]
    fn snp_groups_are_single_linkage() {
        let d = vec![
            vec![0, 3, 30, 31],
            vec![3, 0, 28, 29],
            vec![30, 28, 0, 4],
            vec![31, 29, 4, 0],
        ];
        let g = snp_groups(&d, &[2, 5, 30]);
        assert_eq!(g[0], vec![0, 1, 2, 3]);
        assert_eq!(g[1], vec![0, 0, 1, 1]);
        assert_eq!(g[2], vec![0, 0, 0, 0]);
    }
}
