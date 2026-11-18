//! What a simulation run is asked for (section 10.1).
//!
//! Every field is serde-serialisable so the CLI can drop the whole configuration
//! into `run.json` and the benchmark can replay a dataset exactly.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

/// Where the reference sequence comes from.
///
/// Serialised with serde's external tagging, so the manifest holds
/// `{"random":{"length":2000000,"gc":0.5}}` or `{"file":"/data/ref.fasta"}`.
/// Internal tagging would be prettier to read but cannot express a newtype
/// variant like `File`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RefSpec {
    /// Draw an i.i.d. random sequence. This is what tests T1-T5 use: it has no
    /// repeats, so every difference the engine finds is one the simulator put
    /// there.
    Random {
        /// Length in bases.
        length: usize,
        /// GC fraction in `[0, 1]`.
        gc: f64,
    },
    /// Read a FASTA (plain or gzip) with [`mtato_core::fasta::read_genome`]. Its
    /// contigs are concatenated in file order into one linear reference, because
    /// every simulated genome is cut into [`SimConfig::contigs`] contigs anyway.
    File(PathBuf),
}

impl Default for RefSpec {
    fn default() -> Self {
        RefSpec::Random {
            length: 2_000_000,
            gc: 0.5,
        }
    }
}

/// Gain and loss of large segments, optional per branch.
///
/// A *gain* inserts a novel random segment; a *loss* deletes a stretch of the
/// genome the branch inherited. Both are recorded in the truth with their exact
/// length, which is what test T4 measures the engine against.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SegmentEvents {
    /// Mean number of gained segments per branch (Poisson). 0 disables gains.
    pub mean_gains_per_branch: f64,
    /// Mean number of lost segments per branch (Poisson). 0 disables losses.
    pub mean_losses_per_branch: f64,
    /// Shortest segment, in bases.
    pub min_len: usize,
    /// Longest segment, in bases.
    pub max_len: usize,
}

impl Default for SegmentEvents {
    fn default() -> Self {
        Self {
            mean_gains_per_branch: 0.0,
            mean_losses_per_branch: 0.0,
            min_len: 500,
            max_len: 5000,
        }
    }
}

/// A complete simulator configuration.
///
/// `SimConfig::default()` is the documented default of section 10.1: a 2 Mb
/// random reference at 50% GC, one cluster, ten leaves, 2% divergence between
/// cluster roots, three SNPs per branch, no indels, no segment events, and one
/// forward contig per genome.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SimConfig {
    /// Seed of the whole run. Same seed and same configuration give
    /// byte-identical FASTA and truth files.
    pub seed: u64,
    /// The sequence every cluster root is derived from.
    pub reference: RefSpec,
    /// Number of clusters, `m` in section 10.1.
    pub clusters: usize,
    /// Number of leaf genomes per cluster, `n`.
    pub leaves_per_cluster: usize,
    /// Substitution divergence between the reference and each cluster root, as a
    /// fraction in `[0, 1)`. 0.02 means 2% of reference positions differ.
    pub divergence: f64,
    /// Mean SNPs per branch of a cluster tree (Poisson). Branch lengths in the
    /// truth Newick are these counts, so a pairwise distance is a path sum.
    pub mean_snps_per_branch: f64,
    /// Mean small indels per branch (Poisson). Each one is an insertion or a
    /// deletion with equal probability. 0 disables them.
    pub mean_indels_per_branch: f64,
    /// Shortest small indel, in bases.
    pub indel_min_len: usize,
    /// Longest small indel, in bases.
    pub indel_max_len: usize,
    /// Gain and loss of 500-5000 bp segments, off by default.
    pub segment_events: SegmentEvents,
    /// Cut every genome into this many contigs, in random order. 1 keeps each
    /// genome as a single contig.
    pub contigs: usize,
    /// Reverse-complement each contig independently with probability one half.
    /// Off by default, so an unfragmented run emits plain forward sequences;
    /// T3-style runs want this together with `contigs > 1`.
    pub random_orientation: bool,
    /// GC fraction of the novel sequence the simulator invents (insertions and
    /// gained segments). `None` uses the GC measured on the reference, which is
    /// what keeps inserted segments from standing out as composition outliers.
    pub gc: Option<f64>,
    /// Prefix of every generated name: `<prefix>_c<cluster>_s<leaf>`.
    pub name_prefix: String,
    /// Bases per line in the FASTA files.
    pub fasta_width: usize,
}

impl Default for SimConfig {
    fn default() -> Self {
        Self {
            seed: 0,
            reference: RefSpec::default(),
            clusters: 1,
            leaves_per_cluster: 10,
            divergence: 0.02,
            mean_snps_per_branch: 3.0,
            mean_indels_per_branch: 0.0,
            indel_min_len: 1,
            indel_max_len: 10,
            segment_events: SegmentEvents::default(),
            contigs: 1,
            random_orientation: false,
            gc: None,
            name_prefix: "sim".to_string(),
            fasta_width: 60,
        }
    }
}

/// Above this a Poisson draw costs more than it is worth: the mean would be
/// larger than any bacterial genome, so the configuration is a mistake.
const MAX_MEAN_PER_BRANCH: f64 = 1.0e7;

impl SimConfig {
    /// Reject a configuration that cannot produce a meaningful dataset.
    ///
    /// Pure check on the numbers: nothing is read from disk, so a CLI can call
    /// this before doing any work. Problems that depend on the reference file
    /// (it does not exist, it is shorter than the contig count) are reported by
    /// [`crate::simulate`].
    pub fn validate(&self) -> Result<(), String> {
        if self.clusters == 0 {
            return Err("clusters must be at least 1".to_string());
        }
        if self.leaves_per_cluster == 0 {
            return Err("leaves_per_cluster must be at least 1".to_string());
        }
        if !(0.0..1.0).contains(&self.divergence) {
            return Err(format!(
                "divergence must be in [0, 1), got {}",
                self.divergence
            ));
        }
        check_mean("mean_snps_per_branch", self.mean_snps_per_branch)?;
        check_mean("mean_indels_per_branch", self.mean_indels_per_branch)?;
        check_mean(
            "segment_events.mean_gains_per_branch",
            self.segment_events.mean_gains_per_branch,
        )?;
        check_mean(
            "segment_events.mean_losses_per_branch",
            self.segment_events.mean_losses_per_branch,
        )?;
        if self.indel_min_len == 0 || self.indel_min_len > self.indel_max_len {
            return Err(format!(
                "indel lengths must satisfy 1 <= min <= max, got {} and {}",
                self.indel_min_len, self.indel_max_len
            ));
        }
        if self.segment_events.min_len == 0
            || self.segment_events.min_len > self.segment_events.max_len
        {
            return Err(format!(
                "segment lengths must satisfy 1 <= min <= max, got {} and {}",
                self.segment_events.min_len, self.segment_events.max_len
            ));
        }
        if self.contigs == 0 {
            return Err("contigs must be at least 1".to_string());
        }
        if let RefSpec::Random { length, gc } = self.reference {
            if length == 0 {
                return Err("reference length must be at least 1".to_string());
            }
            if length < self.contigs {
                return Err(format!(
                    "reference length {length} cannot be cut into {} contigs",
                    self.contigs
                ));
            }
            if !(0.0..=1.0).contains(&gc) {
                return Err(format!("reference gc must be in [0, 1], got {gc}"));
            }
        }
        if let Some(gc) = self.gc {
            if !(0.0..=1.0).contains(&gc) {
                return Err(format!("gc must be in [0, 1], got {gc}"));
            }
        }
        if self.fasta_width == 0 {
            return Err("fasta_width must be at least 1".to_string());
        }
        check_name_prefix(&self.name_prefix)?;
        Ok(())
    }

    /// Total number of leaf genomes the run produces.
    pub fn genome_count(&self) -> usize {
        self.clusters * self.leaves_per_cluster
    }
}

fn check_mean(what: &str, mean: f64) -> Result<(), String> {
    if !(mean.is_finite() && mean >= 0.0) {
        return Err(format!("{what} must be finite and >= 0, got {mean}"));
    }
    if mean > MAX_MEAN_PER_BRANCH {
        return Err(format!("{what} = {mean} is absurdly large"));
    }
    Ok(())
}

fn check_name_prefix(prefix: &str) -> Result<(), String> {
    if prefix.is_empty() {
        return Err("name_prefix must not be empty".to_string());
    }
    // Names become file names and Newick labels, so anything structural or
    // platform-hostile is refused up front rather than quoted into obscurity.
    let bad = prefix
        .chars()
        .find(|c| !c.is_ascii_alphanumeric() && !matches!(c, '_' | '-' | '.'));
    match bad {
        Some(c) => Err(format!(
            "name_prefix may only contain ASCII letters, digits, '_', '-' and '.', got {c:?} in {prefix:?}"
        )),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_the_outline() {
        let c = SimConfig::default();
        assert_eq!(c.seed, 0);
        assert_eq!(
            c.reference,
            RefSpec::Random {
                length: 2_000_000,
                gc: 0.5
            }
        );
        assert_eq!(c.clusters, 1);
        assert_eq!(c.leaves_per_cluster, 10);
        assert_eq!(c.divergence, 0.02);
        assert_eq!(c.mean_snps_per_branch, 3.0);
        assert_eq!(c.mean_indels_per_branch, 0.0);
        assert_eq!(c.indel_min_len, 1);
        assert_eq!(c.indel_max_len, 10);
        assert_eq!(c.segment_events.min_len, 500);
        assert_eq!(c.segment_events.max_len, 5000);
        assert_eq!(c.contigs, 1);
        assert!(c.validate().is_ok());
        assert_eq!(c.genome_count(), 10);
    }

    #[test]
    fn every_outline_test_configuration_validates() {
        // T1
        let t1 = SimConfig {
            reference: RefSpec::Random {
                length: 2_000_000,
                gc: 0.5,
            },
            clusters: 1,
            leaves_per_cluster: 50,
            mean_snps_per_branch: 3.0,
            ..SimConfig::default()
        };
        assert!(t1.validate().is_ok());
        // T2
        let t2 = SimConfig {
            mean_indels_per_branch: 0.5,
            ..t1.clone()
        };
        assert!(t2.validate().is_ok());
        // T3
        let t3 = SimConfig {
            contigs: 100,
            ..t2.clone()
        };
        assert!(t3.validate().is_ok());
        // T4
        let t4 = SimConfig {
            segment_events: SegmentEvents {
                mean_gains_per_branch: 0.2,
                mean_losses_per_branch: 0.2,
                min_len: 500,
                max_len: 5000,
            },
            ..t2.clone()
        };
        assert!(t4.validate().is_ok());
        // T5
        let t5 = SimConfig {
            clusters: 5,
            leaves_per_cluster: 20,
            divergence: 0.02,
            ..t2
        };
        assert!(t5.validate().is_ok());
        assert_eq!(t5.genome_count(), 100);
    }

    #[test]
    fn validation_rejects_nonsense() {
        let cases: Vec<(&str, SimConfig)> = vec![
            (
                "no clusters",
                SimConfig {
                    clusters: 0,
                    ..SimConfig::default()
                },
            ),
            (
                "no leaves",
                SimConfig {
                    leaves_per_cluster: 0,
                    ..SimConfig::default()
                },
            ),
            (
                "divergence 1",
                SimConfig {
                    divergence: 1.0,
                    ..SimConfig::default()
                },
            ),
            (
                "divergence -1",
                SimConfig {
                    divergence: -0.1,
                    ..SimConfig::default()
                },
            ),
            (
                "nan divergence",
                SimConfig {
                    divergence: f64::NAN,
                    ..SimConfig::default()
                },
            ),
            (
                "negative snps",
                SimConfig {
                    mean_snps_per_branch: -1.0,
                    ..SimConfig::default()
                },
            ),
            (
                "infinite snps",
                SimConfig {
                    mean_snps_per_branch: f64::INFINITY,
                    ..SimConfig::default()
                },
            ),
            (
                "huge indel mean",
                SimConfig {
                    mean_indels_per_branch: 1.0e9,
                    ..SimConfig::default()
                },
            ),
            (
                "indel min 0",
                SimConfig {
                    indel_min_len: 0,
                    ..SimConfig::default()
                },
            ),
            (
                "indel min > max",
                SimConfig {
                    indel_min_len: 20,
                    indel_max_len: 3,
                    ..SimConfig::default()
                },
            ),
            (
                "segment min > max",
                SimConfig {
                    segment_events: SegmentEvents {
                        min_len: 5000,
                        max_len: 500,
                        ..SegmentEvents::default()
                    },
                    ..SimConfig::default()
                },
            ),
            (
                "no contigs",
                SimConfig {
                    contigs: 0,
                    ..SimConfig::default()
                },
            ),
            (
                "zero length",
                SimConfig {
                    reference: RefSpec::Random { length: 0, gc: 0.5 },
                    ..SimConfig::default()
                },
            ),
            (
                "too many contigs",
                SimConfig {
                    reference: RefSpec::Random {
                        length: 10,
                        gc: 0.5,
                    },
                    contigs: 20,
                    ..SimConfig::default()
                },
            ),
            (
                "bad gc",
                SimConfig {
                    reference: RefSpec::Random {
                        length: 1000,
                        gc: 1.5,
                    },
                    ..SimConfig::default()
                },
            ),
            (
                "bad novel gc",
                SimConfig {
                    gc: Some(-0.2),
                    ..SimConfig::default()
                },
            ),
            (
                "zero width",
                SimConfig {
                    fasta_width: 0,
                    ..SimConfig::default()
                },
            ),
            (
                "empty prefix",
                SimConfig {
                    name_prefix: String::new(),
                    ..SimConfig::default()
                },
            ),
            (
                "slash in prefix",
                SimConfig {
                    name_prefix: "a/b".into(),
                    ..SimConfig::default()
                },
            ),
            (
                "space in prefix",
                SimConfig {
                    name_prefix: "a b".into(),
                    ..SimConfig::default()
                },
            ),
        ];
        for (label, cfg) in cases {
            assert!(cfg.validate().is_err(), "{label} should be rejected");
        }
    }

    #[test]
    fn config_round_trips_through_json() {
        let cfg = SimConfig {
            seed: 12345,
            reference: RefSpec::File(PathBuf::from("/data/ref.fasta")),
            gc: Some(0.38),
            ..SimConfig::default()
        };
        let json = serde_json::to_string(&cfg).unwrap();
        let back: SimConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back, cfg);
        assert_eq!(
            back.reference,
            RefSpec::File(PathBuf::from("/data/ref.fasta"))
        );
        // Floats survive exactly, which the truth manifest depends on.
        let cfg2 = SimConfig {
            divergence: 0.013,
            mean_snps_per_branch: 2.7,
            ..SimConfig::default()
        };
        let back2: SimConfig =
            serde_json::from_str(&serde_json::to_string(&cfg2).unwrap()).unwrap();
        assert_eq!(back2.divergence, 0.013);
        assert_eq!(back2.mean_snps_per_branch, 2.7);
    }

    #[test]
    fn unknown_fields_are_rejected() {
        let bad = r#"{"seed":1,"clusters":1,"leaves_per_cluster":2,"divergence":0.02,
            "mean_snps_per_branch":3.0,"mean_indels_per_branch":0.0,"indel_min_len":1,
            "indel_max_len":10,"segment_events":{"mean_gains_per_branch":0.0,
            "mean_losses_per_branch":0.0,"min_len":500,"max_len":5000},"contigs":1,
            "random_orientation":true,"gc":null,"name_prefix":"sim","fasta_width":60,
            "reference":{"random":{"length":1000,"gc":0.5}},"typo":true}"#;
        assert!(serde_json::from_str::<SimConfig>(bad).is_err());
    }
}
