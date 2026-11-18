//! The "Resolution" and "Accuracy (SNP rules)" settings: how far a run refines
//! and which SNPs it trusts.
//!
//! These mirror parameters the engine is gaining (`resolution`,
//! `snp_thresholds`, `refine_backbone`, `snp_filter`). Until they exist in
//! `infinity_core::params::Parameters`, the window keeps them here, remembers
//! them in the preferences, shows them in the Settings window, and writes them
//! into the parameter set at one place: [`ResolutionSettings::apply_to`]. Wiring
//! them is then a matter of filling in that one function.
//!
//! Nothing in this module touches egui, so it is tested like the rest of the
//! model.

use serde::{Deserialize, Serialize};

use infinity_core::params::Parameters;

/// How far iterative refinement goes (section 5's focus levels as a stop).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Resolution {
    /// Focus 0: cluster the genomes from their sketches and stop.
    Sketch,
    /// Focus 1: also find *where* genomes differ, block by block.
    Blocks,
    /// Focus 2: resolve every difference to exact SNPs and indels.
    #[default]
    ExactSnps,
}

impl Resolution {
    pub const ALL: [Resolution; 3] = [
        Resolution::Sketch,
        Resolution::Blocks,
        Resolution::ExactSnps,
    ];

    /// The short name used on the radio buttons.
    pub fn label(self) -> &'static str {
        match self {
            Resolution::Sketch => "Sketch",
            Resolution::Blocks => "Blocks",
            Resolution::ExactSnps => "Exact SNPs",
        }
    }

    /// What it means, in plain words.
    pub fn explain(self) -> &'static str {
        match self {
            Resolution::Sketch => {
                "Clusters only. Fastest: groups genomes and draws a blurry tree, \
                 no per-base comparison."
            }
            Resolution::Blocks => {
                "Where genomes differ. Compares blocks to find the differing regions; \
                 counts are estimates (≈)."
            }
            Resolution::ExactSnps => {
                "Every difference resolved to exact SNPs and indels. Slowest, and the \
                 only level with exact SNP distances."
            }
        }
    }

    /// The focus level a run with this resolution stops at.
    pub fn focus(self) -> u8 {
        match self {
            Resolution::Sketch => crate::fmt::FOCUS_SKETCH,
            Resolution::Blocks => crate::fmt::FOCUS_BLOCKS,
            Resolution::ExactSnps => crate::fmt::FOCUS_EXACT,
        }
    }
}

/// Rules that decide which SNPs are trustworthy enough to count.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct SnpFilter {
    /// Apply the rules at all.
    pub enabled: bool,
    /// Count only SNPs in positions every genome of the cluster covers.
    pub core_only: bool,
    /// Ignore contigs shorter than this, in bp.
    pub min_contig_len: u64,
    /// Ignore SNPs this close to a contig end, in bp.
    pub contig_end_margin: u64,
    /// A contig whose depth-like copy ratio exceeds this is treated as a
    /// repeat and its SNPs are not counted.
    pub high_copy_ratio: f64,
    /// SNPs packed more densely than chance within this window (bp) are
    /// treated as one recombination or mis-assembly event.
    pub cluster_window: u64,
    /// Significance level of the SNP-cluster test.
    pub alpha: f64,
}

impl Default for SnpFilter {
    fn default() -> Self {
        Self {
            enabled: true,
            core_only: true,
            min_contig_len: 5_000,
            contig_end_margin: 150,
            high_copy_ratio: 1.8,
            cluster_window: 1_000,
            alpha: 0.01,
        }
    }
}

/// Everything the "Resolution" and "Accuracy (SNP rules)" sections edit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ResolutionSettings {
    pub resolution: Resolution,
    /// SNP cut-offs for outbreak-style SNP clusters, ascending.
    pub snp_thresholds: Vec<u32>,
    /// Also sharpen the distances *between* clusters, not only inside them.
    pub refine_backbone: bool,
    pub snp_filter: SnpFilter,
}

impl Default for ResolutionSettings {
    fn default() -> Self {
        Self {
            resolution: Resolution::ExactSnps,
            snp_thresholds: vec![5, 10, 20],
            refine_backbone: true,
            snp_filter: SnpFilter::default(),
        }
    }
}

/// The largest SNP threshold the Settings window accepts. Larger cut-offs
/// would merge whole lineages, which is what the cluster distance is for.
pub const MAX_SNP_THRESHOLD: u32 = 100_000;

impl ResolutionSettings {
    /// Sort the thresholds, drop duplicates and zeros. Called after every edit
    /// so that the list on screen is always the list the engine would get.
    pub fn normalise(&mut self) {
        self.snp_thresholds
            .retain(|&t| t > 0 && t <= MAX_SNP_THRESHOLD);
        self.snp_thresholds.sort_unstable();
        self.snp_thresholds.dedup();
    }

    /// Add one threshold from text typed by the user. Returns a message for
    /// the user when the text is not a usable cut-off.
    pub fn add_threshold(&mut self, text: &str) -> Result<(), String> {
        let t = text.trim();
        let v: u32 = t
            .parse()
            .map_err(|_| format!("\"{t}\" is not a whole number of SNPs"))?;
        if v == 0 {
            return Err("a threshold of 0 SNPs would split identical genomes".to_string());
        }
        if v > MAX_SNP_THRESHOLD {
            return Err(format!("at most {MAX_SNP_THRESHOLD} SNPs"));
        }
        self.snp_thresholds.push(v);
        self.normalise();
        Ok(())
    }

    /// Reject values that would make the rules meaningless. Messages are
    /// written for the person reading the Settings window.
    pub fn validate(&self) -> Result<(), String> {
        let f = &self.snp_filter;
        if !(f.high_copy_ratio >= 1.0 && f.high_copy_ratio.is_finite()) {
            return Err(format!(
                "high-copy contig ratio must be at least 1, got {}",
                f.high_copy_ratio
            ));
        }
        if !(f.alpha > 0.0 && f.alpha < 1.0) {
            return Err(format!("alpha must be between 0 and 1, got {}", f.alpha));
        }
        if f.cluster_window == 0 {
            return Err("the SNP-cluster window must be at least 1 bp".to_string());
        }
        if self.snp_thresholds.windows(2).any(|w| w[0] >= w[1]) {
            return Err("SNP thresholds must ascend without repeats".to_string());
        }
        Ok(())
    }

    /// Write these settings into a run's parameter set.
    ///
    /// The engine does not have these fields yet, so this is the single place
    /// to wire them: everything else in the window already treats them as real.
    pub fn apply_to(&self, params: &mut Parameters) {
        // WIRE: once `infinity_core::params::Parameters` has these fields:
        //   params.resolution = self.resolution.into();
        //   params.snp_thresholds = self.snp_thresholds.clone();
        //   params.refine_backbone = self.refine_backbone;
        //   params.snp_filter = self.snp_filter.clone().into();
        let _ = params;
    }
}

/// Cluster distance as the percent identity a reader thinks in: 0.005 is
/// 99.5 % identity.
pub fn identity_pct(cluster_dist: f64) -> f64 {
    (1.0 - cluster_dist) * 100.0
}

/// The inverse of [`identity_pct`], clamped to the range the slider allows.
pub fn dist_from_identity(pct: f64) -> f64 {
    (1.0 - pct / 100.0).clamp(0.000_01, 0.99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_the_documented_ones() {
        let s = ResolutionSettings::default();
        assert_eq!(s.resolution, Resolution::ExactSnps);
        assert_eq!(s.snp_thresholds, vec![5, 10, 20]);
        assert!(s.refine_backbone);
        let f = &s.snp_filter;
        assert!(f.enabled && f.core_only);
        assert_eq!(f.min_contig_len, 5_000);
        assert_eq!(f.contig_end_margin, 150);
        assert_eq!(f.high_copy_ratio, 1.8);
        assert_eq!(f.cluster_window, 1_000);
        assert_eq!(f.alpha, 0.01);
        s.validate().expect("the defaults are valid");
    }

    #[test]
    fn thresholds_stay_sorted_and_unique() {
        let mut s = ResolutionSettings::default();
        s.add_threshold("15").unwrap();
        s.add_threshold(" 5 ").unwrap();
        s.add_threshold("2").unwrap();
        assert_eq!(s.snp_thresholds, vec![2, 5, 10, 15, 20]);
        assert!(s.add_threshold("0").is_err());
        assert!(s.add_threshold("ten").is_err());
        assert!(s.add_threshold("-3").is_err());
        assert!(s.add_threshold("1000000").is_err());
        assert_eq!(s.snp_thresholds, vec![2, 5, 10, 15, 20]);
        s.snp_thresholds = vec![20, 0, 5, 5];
        s.normalise();
        assert_eq!(s.snp_thresholds, vec![5, 20]);
    }

    #[test]
    fn nonsense_rules_are_rejected_with_a_reason() {
        let mut s = ResolutionSettings::default();
        s.snp_filter.alpha = 0.0;
        assert!(s.validate().unwrap_err().contains("alpha"));
        s.snp_filter.alpha = 0.05;
        s.snp_filter.high_copy_ratio = 0.5;
        assert!(s.validate().unwrap_err().contains("ratio"));
        s.snp_filter.high_copy_ratio = 2.0;
        s.snp_filter.cluster_window = 0;
        assert!(s.validate().unwrap_err().contains("window"));
        s.snp_filter.cluster_window = 10;
        s.snp_thresholds = vec![10, 5];
        assert!(s.validate().unwrap_err().contains("ascend"));
    }

    #[test]
    fn identity_and_distance_are_inverses() {
        assert!((identity_pct(0.005) - 99.5).abs() < 1e-9);
        for d in [0.0001, 0.005, 0.02, 0.5] {
            assert!((dist_from_identity(identity_pct(d)) - d).abs() < 1e-12);
        }
        // Out-of-range identities are clamped to what the slider allows.
        assert_eq!(dist_from_identity(100.0), 0.000_01);
        assert_eq!(dist_from_identity(0.0), 0.99);
    }

    #[test]
    fn resolutions_map_to_focus_levels_and_have_words() {
        assert_eq!(Resolution::Sketch.focus(), 0);
        assert_eq!(Resolution::Blocks.focus(), 1);
        assert_eq!(Resolution::ExactSnps.focus(), 2);
        for r in Resolution::ALL {
            assert!(!r.label().is_empty() && r.explain().len() > 20);
        }
    }

    #[test]
    fn settings_round_trip_and_tolerate_missing_fields() {
        let mut s = ResolutionSettings::default();
        s.resolution = Resolution::Blocks;
        s.snp_filter.core_only = false;
        let j = serde_json::to_string(&s).unwrap();
        assert!(j.contains("\"blocks\""), "{j}");
        let back: ResolutionSettings = serde_json::from_str(&j).unwrap();
        assert_eq!(back, s);
        // A file from before the filter existed still loads, with defaults.
        let old: ResolutionSettings = serde_json::from_str(r#"{"resolution":"sketch"}"#).unwrap();
        assert_eq!(old.resolution, Resolution::Sketch);
        assert_eq!(old.snp_filter, SnpFilter::default());
    }

    #[test]
    fn applying_leaves_the_existing_parameters_valid() {
        let s = ResolutionSettings::default();
        let mut p = Parameters::default();
        s.apply_to(&mut p);
        p.validate().expect("still valid");
    }
}
