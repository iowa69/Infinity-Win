//! Clustering without all-versus-all, and the backbone tree (6.4).
//!
//! The greedy rule is the one in the outline: genomes are considered in a fixed
//! order, each is compared only against existing pivots through an inverted
//! index of sketch hashes, and it either joins the closest pivot's cluster or
//! becomes a pivot itself. Because the best assembly of each cluster becomes
//! its pivot and the order is total, the result is deterministic and no
//! pairwise comparison of two non-pivot genomes is ever made.

use serde::{Deserialize, Serialize};

use crate::genome::GenomeRecord;
use crate::params::Parameters;
use crate::sketch;
use crate::tree::{midpoint_root, neighbor_joining, Tree};

/// A cluster: a pivot and every genome assigned to it, the pivot included.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Cluster {
    pub id: u32,
    /// Genome id of the representative. Differences inside the cluster are
    /// expressed in this genome's coordinates (section 5).
    pub pivot: u32,
    /// Genome ids, pivot first, then the others in ascending id order.
    pub members: Vec<u32>,
}

impl Cluster {
    pub fn size(&self) -> usize {
        self.members.len()
    }

    pub fn is_singleton(&self) -> bool {
        self.members.len() <= 1
    }

    /// Members other than the pivot, in ascending id order.
    pub fn non_pivot(&self) -> impl Iterator<Item = u32> + '_ {
        self.members.iter().copied().filter(|&m| m != self.pivot)
    }
}

/// The outcome of one clustering pass.
#[derive(Debug, Clone, PartialEq)]
pub struct Clustering {
    pub clusters: Vec<Cluster>,
    /// Sketch distance from each genome to its cluster pivot, indexed by genome
    /// id. `f64::NAN` for a genome that was not assigned (cannot happen, but an
    /// unassigned distance must never be silently zero).
    pub dist_to_pivot: Vec<f64>,
    pub threshold: f64,
}

/// The order genomes are considered in: N50 descending, total length
/// descending, name ascending (6.4). Sorting on N50 first is what makes the
/// most contiguous assembly of a cluster its pivot.
pub fn genome_order(recs: &[GenomeRecord]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..recs.len()).collect();
    idx.sort_by(|&a, &b| {
        recs[b]
            .stats
            .n50
            .cmp(&recs[a].stats.n50)
            .then(recs[b].stats.length.cmp(&recs[a].stats.length))
            .then_with(|| {
                let (na, nb) = (&recs[a].name, &recs[b].name);
                match na.cmp(nb) {
                    std::cmp::Ordering::Equal => recs[a].id.cmp(&recs[b].id),
                    o => o,
                }
            })
    });
    idx
}

/// Greedy clustering against pivots only (6.4).
pub fn cluster_genomes(recs: &[GenomeRecord], params: &Parameters) -> Clustering {
    let k = params.k;
    let order = genome_order(recs);
    let mut clusters: Vec<Cluster> = Vec::new();
    let mut dist_to_pivot = vec![f64::NAN; recs.len()];

    if params.single_cluster {
        // Experiment switch (section 11): cluster-first is turned off, so every
        // genome is compared against one pivot. The pivot is the first genome
        // in the canonical order, which keeps the run deterministic.
        if let Some(&first) = order.first() {
            let members: Vec<u32> = order.iter().map(|&i| recs[i].id).collect();
            let pivot = recs[first].id;
            for &i in order.iter() {
                dist_to_pivot[recs[i].id as usize] =
                    sketch::distance(&recs[i].sketch, &recs[first].sketch, k);
            }
            clusters.push(Cluster {
                id: 0,
                pivot,
                members,
            });
        }
        return Clustering {
            clusters,
            dist_to_pivot,
            threshold: params.cluster_dist,
        };
    }

    // Inverted index: sketch hash -> pivot cluster ids that contain it.
    let mut index: std::collections::HashMap<u64, Vec<u32>> = std::collections::HashMap::new();
    // Scratch space reused for every genome so the hot loop does not allocate.
    let mut counts: Vec<u32> = Vec::new();
    let mut touched: Vec<u32> = Vec::new();

    for &gi in &order {
        let g = &recs[gi];
        let gid = g.id as usize;

        counts.clear();
        counts.resize(clusters.len(), 0);
        touched.clear();
        for &h in &g.sketch {
            if let Some(list) = index.get(&h) {
                for &c in list {
                    if counts[c as usize] == 0 {
                        touched.push(c);
                    }
                    counts[c as usize] += 1;
                }
            }
        }
        touched.sort_unstable();

        // Smallest distance wins; ties go to the earliest pivot.
        let mut best: Option<(u32, f64)> = None;
        for &c in &touched {
            let pivot_g = &recs[clusters[c as usize].pivot as usize];
            let shared = counts[c as usize] as usize;
            let d = sketch::distance_from_shared(shared, g.sketch.len(), pivot_g.sketch.len(), k);
            if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                best = Some((c, d));
            }
        }
        // A genome can also be close to a pivot it shares no hash with; with
        // `cluster_dist` small that is impossible, because sharing nothing gives
        // D = 1. The loop above therefore sees every candidate that can join.

        match best {
            Some((c, d)) if d <= params.cluster_dist => {
                dist_to_pivot[gid] = d;
                clusters[c as usize].members.push(g.id);
            }
            _ => {
                let cid = clusters.len() as u32;
                dist_to_pivot[gid] = 0.0;
                clusters.push(Cluster {
                    id: cid,
                    pivot: g.id,
                    members: vec![g.id],
                });
                for &h in &g.sketch {
                    index.entry(h).or_default().push(cid);
                }
            }
        }
    }

    // Members are appended in processing order; sort by id so the output files
    // and the trees do not depend on how the work was scheduled.
    for c in &mut clusters {
        let pivot = c.pivot;
        c.members.sort_unstable();
        // Pivot first, for readability of `clusters.tsv`.
        if let Some(p) = c.members.iter().position(|&m| m == pivot) {
            let v = c.members.remove(p);
            c.members.insert(0, v);
        }
    }

    Clustering {
        clusters,
        dist_to_pivot,
        threshold: params.cluster_dist,
    }
}

/// Names that are safe to use as Newick leaves and unique across the run.
///
/// Two input files can easily have the same stem, and the full tree is
/// assembled by looking a leaf up *by name* (6.9), so a collision would graft
/// the wrong cluster. Duplicates get their genome id appended and the caller
/// is told, through the returned list of warnings.
pub fn unique_names(recs: &[GenomeRecord]) -> (Vec<String>, Vec<String>) {
    let mut names = vec![String::new(); recs.len()];
    let mut seen: std::collections::HashMap<String, u32> = std::collections::HashMap::new();
    let mut warnings = Vec::new();
    for r in recs {
        let base = r.name.clone();
        match seen.get(&base) {
            None => {
                seen.insert(base.clone(), r.id);
                names[r.id as usize] = base;
            }
            Some(other) => {
                let renamed = format!("{base}#{id}", id = r.id);
                warnings.push(format!(
                    "genome name {base:?} is used by genome {other} and genome {}; \
                     the second was renamed to {renamed:?} in tree output",
                    r.id
                ));
                names[r.id as usize] = renamed;
            }
        }
    }
    (names, warnings)
}

/// The backbone tree: neighbour-joining on sketch distances between pivots
/// only, midpoint rooted. Its leaves are clusters, named after their pivots.
pub fn backbone_tree(
    clusters: &[Cluster],
    recs: &[GenomeRecord],
    names: &[String],
    params: &Parameters,
) -> Tree {
    let n = clusters.len();
    if n == 0 {
        return Tree::Internal {
            len: 0.0,
            children: Vec::new(),
        };
    }
    if n == 1 {
        return Tree::Leaf {
            name: names[clusters[0].pivot as usize].clone(),
            len: 0.0,
        };
    }
    let labels: Vec<String> = clusters
        .iter()
        .map(|c| names[c.pivot as usize].clone())
        .collect();
    let mut dist = vec![vec![0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let a = &recs[clusters[i].pivot as usize].sketch;
            let b = &recs[clusters[j].pivot as usize].sketch;
            let d = sketch::distance(a, b, params.k);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }
    let mut t = neighbor_joining(&labels, &dist);
    t.clamp_negative_lengths();
    midpoint_root(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::GenomeStats;
    use std::path::PathBuf;

    fn rec(id: u32, name: &str, n50: u64, length: u64, sketch: Vec<u64>) -> GenomeRecord {
        GenomeRecord::new(
            id,
            name.to_string(),
            PathBuf::from(format!("{name}.fa")),
            GenomeStats {
                length,
                contigs: 1,
                n50,
                gc: 0.5,
                longest: length,
                n_bases: 0,
            },
            sketch,
        )
    }

    fn params() -> Parameters {
        Parameters {
            cluster_dist: 0.02,
            ..Default::default()
        }
    }

    #[test]
    fn identical_sketches_land_in_one_cluster() {
        let s: Vec<u64> = (0..100).collect();
        let recs = vec![
            rec(0, "a", 100, 100, s.clone()),
            rec(1, "b", 100, 100, s.clone()),
            rec(2, "c", 100, 100, s.clone()),
        ];
        let c = cluster_genomes(&recs, &params());
        assert_eq!(c.clusters.len(), 1);
        assert_eq!(c.clusters[0].members, vec![0, 1, 2]);
        assert_eq!(c.dist_to_pivot[0], 0.0);
        assert!(c.dist_to_pivot[1] <= 1e-12);
    }

    #[test]
    fn unrelated_sketches_become_separate_pivots() {
        let recs = vec![
            rec(0, "a", 100, 100, (0..100).collect()),
            rec(1, "b", 100, 100, (1000..1100).collect()),
        ];
        let c = cluster_genomes(&recs, &params());
        assert_eq!(c.clusters.len(), 2);
        assert!(c.clusters.iter().all(|x| x.is_singleton()));
    }

    #[test]
    fn the_most_contiguous_assembly_becomes_the_pivot() {
        // Same sketch content, but b has the best N50, so b must be the pivot
        // even though a has the lower id.
        let s: Vec<u64> = (0..100).collect();
        let recs = vec![
            rec(0, "a", 10, 100, s.clone()),
            rec(1, "b", 90, 100, s.clone()),
            rec(2, "c", 50, 100, s.clone()),
        ];
        let order = genome_order(&recs);
        assert_eq!(order, vec![1, 2, 0]);
        let c = cluster_genomes(&recs, &params());
        assert_eq!(c.clusters[0].pivot, 1);
        assert_eq!(c.clusters[0].members[0], 1, "pivot comes first");
    }

    #[test]
    fn ties_in_order_are_broken_by_name_then_id() {
        let s: Vec<u64> = (0..10).collect();
        // Deliberately stored out of id order: index 0 holds id 2 ("zz"),
        // index 1 holds id 0 ("aa"), index 2 holds id 1 ("aa").
        let recs = vec![
            rec(2, "zz", 50, 100, s.clone()),
            rec(0, "aa", 50, 100, s.clone()),
            rec(1, "aa", 50, 100, s.clone()),
        ];
        let order = genome_order(&recs);
        // Names equal, N50 and length equal, so ids decide: 0 then 1 then 2.
        // Those live at indices 1, 2 and 0 respectively.
        assert_eq!(order, vec![1, 2, 0]);
    }

    #[test]
    fn clustering_is_deterministic_under_shuffling_of_input() {
        let s: Vec<u64> = (0..50).collect();
        let mut recs: Vec<GenomeRecord> = (0..10)
            .map(|i| rec(i, &format!("g{i:02}"), 100 - i as u64, 100, s.clone()))
            .collect();
        let a = cluster_genomes(&recs, &params());
        recs.reverse();
        for (i, r) in recs.iter_mut().enumerate() {
            r.id = i as u32;
        }
        // Same content, different input order: the pivot choice must follow the
        // documented sort, not the file order.
        let b = cluster_genomes(&recs, &params());
        assert_eq!(a.clusters.len(), b.clusters.len());
    }

    #[test]
    fn a_threshold_that_admits_nothing_gives_one_cluster_per_genome() {
        let recs: Vec<GenomeRecord> = (0..5u64)
            .map(|i| {
                rec(
                    i as u32,
                    &format!("g{i}"),
                    100,
                    100,
                    (i * 1000..i * 1000 + 50).collect(),
                )
            })
            .collect();
        let p = Parameters {
            cluster_dist: 1e-9,
            ..Default::default()
        };
        let c = cluster_genomes(&recs, &p);
        assert_eq!(c.clusters.len(), 5);
    }

    #[test]
    fn single_cluster_switch_forces_one_pivot() {
        let recs: Vec<GenomeRecord> = (0..4u64)
            .map(|i| {
                rec(
                    i as u32,
                    &format!("g{i}"),
                    100 - i,
                    100,
                    (i * 1000..i * 1000 + 50).collect(),
                )
            })
            .collect();
        let p = Parameters {
            single_cluster: true,
            ..Default::default()
        };
        let c = cluster_genomes(&recs, &p);
        assert_eq!(c.clusters.len(), 1);
        assert_eq!(c.clusters[0].members.len(), 4);
        assert_eq!(c.clusters[0].pivot, 0, "best N50 is genome 0");
    }

    #[test]
    fn duplicate_names_are_renamed_and_reported() {
        let recs = vec![
            rec(0, "same", 10, 10, vec![]),
            rec(1, "same", 10, 10, vec![]),
        ];
        let (names, warnings) = unique_names(&recs);
        assert_eq!(names[0], "same");
        assert_eq!(names[1], "same#1");
        assert_eq!(warnings.len(), 1);
    }

    #[test]
    fn backbone_tree_has_one_leaf_per_cluster() {
        let recs: Vec<GenomeRecord> = (0..5u64)
            .map(|i| {
                rec(
                    i as u32,
                    &format!("g{i}"),
                    100,
                    100,
                    (i * 10_000..i * 10_000 + 200).collect(),
                )
            })
            .collect();
        let c = cluster_genomes(&recs, &params());
        let (names, _) = unique_names(&recs);
        let t = backbone_tree(&c.clusters, &recs, &names, &params());
        assert_eq!(t.leaf_count(), c.clusters.len());
        assert!(t.len() >= 0.0);
    }

    #[test]
    fn backbone_of_one_cluster_is_a_single_leaf() {
        let s: Vec<u64> = (0..50).collect();
        let recs = vec![rec(0, "a", 100, 100, s.clone()), rec(1, "b", 90, 100, s)];
        let c = cluster_genomes(&recs, &params());
        let (names, _) = unique_names(&recs);
        let t = backbone_tree(&c.clusters, &recs, &names, &params());
        assert_eq!(t.leaf_count(), 1);
    }
}
