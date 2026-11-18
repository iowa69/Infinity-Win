//! Parameters of section 6.1, with validation and the small-genome rule.

use serde::{Deserialize, Serialize};

use crate::seq::LEVELS;

/// Everything that changes the result of a run.
///
/// The five `bool` switches at the bottom are the experiment switches of
/// section 11: each one removes a part of the idea while leaving the reported
/// variants, distances and trees identical. They are stored as positive flags
/// ("use this optimisation") so that `Parameters::default()` is the full idea.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameters {
    /// Anchor k-mer length. Must be odd and at most 31 (6.1).
    pub k: usize,
    /// One anchor per this many k-mers, at levels 1, 2 and 3. Descending, and
    /// each must divide the previous.
    pub densities: [u64; LEVELS],
    /// Level whose anchor hashes form the sketch (1-based).
    pub sketch_level: u8,
    /// Sketch distance threshold for joining a cluster.
    pub cluster_dist: f64,
    /// Collinearity limit: span ratio (6.6).
    pub max_span_ratio: f64,
    /// Collinearity limit: slack added to the ratio test, in bases (6.6).
    pub max_span_extra: u64,
    /// Longest trimmed piece aligned with dynamic programming (6.7).
    pub dp_max_len: usize,
    /// Below this identity an aligned piece is one COMPLEX variant (6.7).
    pub min_identity: f64,
    /// Stop rule for ungapped extension at contig ends (6.6 step D).
    pub xdrop: i32,
    /// Worker threads. 0 means "all cores".
    pub threads: usize,

    // --- experiment switches (section 11) ---
    /// Hash blocks and skip the identical ones (6.6 step A).
    pub block_skip: bool,
    /// Use the nested levels top-down instead of leaf blocks only (6.10).
    pub hierarchical: bool,
    /// Resolve an identical difference piece once (6.10).
    pub piece_cache: bool,
    /// Handle members with identical profiles once (6.9, 6.10).
    pub profile_grouping: bool,
    /// Compare every genome to a single pivot instead of clustering first.
    pub single_cluster: bool,
    /// Sequences of non-pivot genomes are kept while they fit in this budget
    /// (6.10). 0 means "never keep".
    ///
    /// Not yet wired up: the 6.10 eviction it describes is unimplemented, so
    /// this value is parsed, defaulted, validated and written to `run.json`,
    /// but nothing reads it. The only memory control that actually operates is
    /// the piece cache — a member whose difference pieces are all cached is not
    /// re-read for Focus 2 (`engine::resolve_focus1`). See DECISIONS.md.
    pub memory_budget_bytes: u64,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            k: 31,
            densities: [16384, 2048, 256],
            sketch_level: 2,
            cluster_dist: 0.005,
            max_span_ratio: 10.0,
            max_span_extra: 50_000,
            dp_max_len: 5000,
            min_identity: 0.70,
            xdrop: 15,
            threads: 0,
            block_skip: true,
            hierarchical: true,
            piece_cache: true,
            profile_grouping: true,
            single_cluster: false,
            memory_budget_bytes: 0,
        }
    }
}

/// Densities substituted when the median genome is under 500 kb (6.1).
pub const SMALL_GENOME_DENSITIES: [u64; LEVELS] = [1024, 128, 16];

/// Median genome length below which the default densities are too sparse.
pub const SMALL_GENOME_LENGTH: u64 = 500_000;

impl Parameters {
    /// A quarter of physical RAM, or 1 GiB when that cannot be determined
    /// (6.10 memory budget default).
    pub fn default_memory_budget() -> u64 {
        quarter_of_physical_ram().unwrap_or(1 << 30)
    }

    /// Reject anything that would make the algorithm ill-defined. The messages
    /// are written for a user typing them on the command line.
    pub fn validate(&self) -> Result<(), String> {
        if self.k == 0 || self.k % 2 == 0 {
            return Err(format!("k must be odd, got {}", self.k));
        }
        if self.k > 31 {
            return Err(format!(
                "k must be at most 31 to fit in 64 bits, got {}",
                self.k
            ));
        }
        if self.densities.contains(&0) {
            return Err("densities must all be at least 1".to_string());
        }
        for w in self.densities.windows(2) {
            if w[1] >= w[0] {
                return Err(format!(
                    "densities must descend (coarse to fine), got {:?}",
                    self.densities
                ));
            }
            if w[0] % w[1] != 0 {
                return Err(format!(
                    "each density must divide the previous one, got {:?}",
                    self.densities
                ));
            }
        }
        if !(1..=LEVELS as u8).contains(&self.sketch_level) {
            return Err(format!(
                "sketch_level must be between 1 and {LEVELS}, got {}",
                self.sketch_level
            ));
        }
        if !(self.cluster_dist > 0.0 && self.cluster_dist < 1.0) {
            return Err(format!(
                "cluster_dist must be in (0, 1), got {}",
                self.cluster_dist
            ));
        }
        if self.max_span_ratio < 1.0 {
            return Err(format!(
                "max_span_ratio must be at least 1, got {}",
                self.max_span_ratio
            ));
        }
        if !(self.min_identity > 0.0 && self.min_identity <= 1.0) {
            return Err(format!(
                "min_identity must be in (0, 1], got {}",
                self.min_identity
            ));
        }
        if self.xdrop < 1 {
            return Err(format!("xdrop must be at least 1, got {}", self.xdrop));
        }
        Ok(())
    }

    /// Apply the small-genome rule of 6.1: if the median genome length is under
    /// 500 kb the defaults are too sparse, so switch to 1024/128/16. Returns
    /// `Some(message)` when the substitution happened, so the caller can log it.
    pub fn apply_small_genome_rule(&mut self, median_length: u64) -> Option<String> {
        if median_length < SMALL_GENOME_LENGTH && self.densities != SMALL_GENOME_DENSITIES {
            let old = self.densities;
            self.densities = SMALL_GENOME_DENSITIES;
            return Some(format!(
                "median genome length {median_length} bp is under {SMALL_GENOME_LENGTH} bp: \
                 densities changed from {old:?} to {:?}",
                SMALL_GENOME_DENSITIES
            ));
        }
        None
    }

    /// Number of worker threads to actually use.
    pub fn resolved_threads(&self) -> usize {
        if self.threads == 0 {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        } else {
            self.threads
        }
    }
}

/// Total physical memory, from `/proc/meminfo` on Linux.
///
/// Other platforms get `None` and the caller falls back to a fixed, deliberately
/// conservative budget. Reading memory on Windows needs an FFI call, and this
/// crate is `#![forbid(unsafe_code)]`; trading a crate-wide safety guarantee for
/// one tuning number is not a good deal. The fallback only affects how many
/// non-pivot sequences are cached before being re-read from disk (6.10), so it
/// changes memory use and never a result.
fn quarter_of_physical_ram() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let txt = std::fs::read_to_string("/proc/meminfo").ok()?;
        for line in txt.lines() {
            if let Some(rest) = line.strip_prefix("MemTotal:") {
                let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
                return Some(kb.saturating_mul(1024) / 4);
            }
        }
        None
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_the_ones_in_section_6_1() {
        let p = Parameters::default();
        assert_eq!(p.k, 31);
        assert_eq!(p.densities, [16384, 2048, 256]);
        assert_eq!(p.sketch_level, 2);
        assert_eq!(p.cluster_dist, 0.005);
        assert_eq!(p.max_span_ratio, 10.0);
        assert_eq!(p.max_span_extra, 50_000);
        assert_eq!(p.dp_max_len, 5000);
        assert_eq!(p.min_identity, 0.70);
        assert_eq!(p.xdrop, 15);
        p.validate().expect("the documented defaults must be valid");
    }

    #[test]
    fn every_switch_starts_on() {
        let p = Parameters::default();
        assert!(p.block_skip && p.hierarchical && p.piece_cache && p.profile_grouping);
        assert!(!p.single_cluster);
    }

    #[test]
    fn k_must_be_odd_and_small() {
        let mut p = Parameters {
            k: 30,
            ..Default::default()
        };
        assert!(p.validate().is_err());
        p.k = 33;
        assert!(p.validate().is_err());
        p.k = 21;
        assert!(p.validate().is_ok());
    }

    #[test]
    fn densities_must_descend_and_divide() {
        let mut p = Parameters {
            densities: [2048, 16384, 256],
            ..Default::default()
        };
        assert!(p.validate().is_err());
        p.densities = [16384, 3000, 256];
        assert!(p.validate().is_err(), "3000 does not divide 16384");
        p.densities = [16384, 2048, 0];
        assert!(p.validate().is_err());
        p.densities = [16384, 2048, 256];
        assert!(p.validate().is_ok());
    }

    #[test]
    fn small_genomes_switch_densities_and_say_so() {
        let mut p = Parameters::default();
        let msg = p.apply_small_genome_rule(499_999);
        assert!(msg.is_some());
        assert_eq!(p.densities, SMALL_GENOME_DENSITIES);

        let mut p2 = Parameters::default();
        assert!(p2.apply_small_genome_rule(500_000).is_none());
        assert_eq!(p2.densities, [16384, 2048, 256]);

        // Applying it twice must not re-log.
        assert!(p.apply_small_genome_rule(100).is_none());
    }

    #[test]
    fn thread_count_zero_means_all_cores() {
        let p = Parameters {
            threads: 0,
            ..Default::default()
        };
        assert!(p.resolved_threads() >= 1);
        let p2 = Parameters {
            threads: 3,
            ..Default::default()
        };
        assert_eq!(p2.resolved_threads(), 3);
    }
}
