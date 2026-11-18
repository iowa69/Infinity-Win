//! The run summary written to `run.json` (section 9) and carried by the
//! `Finished` event.

use serde::{Deserialize, Serialize};

use crate::counters::Counters;
use crate::params::Parameters;

/// Milliseconds per stage. `None` means the stage never completed, which is
/// different from zero and must never be printed as a number.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StageTimings {
    pub scan_ms: Option<u64>,
    pub cluster_ms: Option<u64>,
    pub backbone_ms: Option<u64>,
    pub refine_ms: Option<u64>,
    pub output_ms: Option<u64>,
    pub total_ms: Option<u64>,
    /// Elapsed time from start to the first tree the user can see (section 11:
    /// "time to first tree").
    pub first_tree_ms: Option<u64>,
    /// Elapsed time until every cluster reached Focus 2 ("time to exact tree").
    /// `None` when the run was cancelled or some cluster never got exact.
    pub exact_tree_ms: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunSummary {
    /// Software version, so a result can be traced back to the code that made it.
    pub version: String,
    pub params: Parameters,
    /// Resolved parameter values that differ from the request, with reasons.
    pub notes: Vec<String>,
    pub genomes: u32,
    pub clusters: u32,
    pub total_length: u64,
    /// How many genomes reached each focus level: index 0 is Focus 0.
    pub focus_counts: [u32; 3],
    pub timings: StageTimings,
    pub counters: Counters,
    /// Peak resident set size, when the platform reports it.
    pub peak_memory_bytes: Option<u64>,
    pub cancelled: bool,
}

impl RunSummary {
    pub fn new(version: &str, params: Parameters) -> Self {
        Self {
            version: version.to_string(),
            params,
            notes: Vec::new(),
            genomes: 0,
            clusters: 0,
            total_length: 0,
            focus_counts: [0, 0, 0],
            timings: StageTimings::default(),
            counters: Counters::default(),
            peak_memory_bytes: peak_memory_bytes(),
            cancelled: false,
        }
    }
}

/// Peak resident set size of this process, or `None` where the platform does
/// not expose it cheaply. Section 11 needs it, and section 0 rule 8 means an
/// unknown value must stay unknown rather than be reported as zero.
pub fn peak_memory_bytes() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let txt = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in txt.lines() {
            if let Some(rest) = line.strip_prefix("VmHWM:") {
                let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
                return Some(kb.saturating_mul(1024));
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
    fn summary_round_trips_through_json() {
        let mut s = RunSummary::new("test", Parameters::default());
        s.genomes = 3;
        s.clusters = 1;
        s.timings.first_tree_ms = Some(120);
        let j = serde_json::to_string(&s).unwrap();
        let back: RunSummary = serde_json::from_str(&j).unwrap();
        assert_eq!(back, s);
    }

    #[test]
    fn unknown_timing_stays_unknown() {
        let s = RunSummary::new("test", Parameters::default());
        assert_eq!(s.timings.exact_tree_ms, None);
        let j = serde_json::to_string(&s).unwrap();
        assert!(j.contains("\"exact_tree_ms\":null"), "{j}");
    }
}
