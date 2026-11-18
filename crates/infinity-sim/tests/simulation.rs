//! Truth tests for the simulator: the nine properties that make the datasets it
//! writes usable as ground truth, plus the scanner property of section 10.2 seen
//! from the simulator side.
//!
//! Everything here goes through the public API only, the same way a benchmark
//! script would.

use std::collections::{btree_map, BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use infinity_core::fasta::read_genome;
use infinity_core::params::Parameters;
use infinity_core::seq::revcomp;
use infinity_core::tree::parse_newick;
use infinity_sim::{
    pairwise_from_newick, simulate, MutationKind, RefSource, RefSpec, SegmentEvents, SimConfig,
    Truth,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// A distinct event of a cluster's pedigree, keyed by where it starts and what
/// kind it is: the end of the reference interval it touches, and the bases it
/// removed and added.
type EventKey = (usize, MutationKind);
type EventSeen = (usize, Vec<u8>, Vec<u8>);

fn tmpdir(tag: &str) -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("infinity-sim-it-{}-{}", std::process::id(), tag));
    let _ = fs::remove_dir_all(&p);
    fs::create_dir_all(&p).unwrap();
    p
}

/// Every file under `dir`, keyed by its path relative to `dir`.
fn all_files(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut out = BTreeMap::new();
    collect(dir, dir, &mut out);
    out
}

fn collect(root: &Path, dir: &Path, out: &mut BTreeMap<String, Vec<u8>>) {
    let mut entries: Vec<PathBuf> = fs::read_dir(dir)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", dir.display()))
        .map(|e| e.unwrap().path())
        .collect();
    entries.sort();
    for path in entries {
        if path.is_dir() {
            collect(root, &path, out);
        } else {
            let rel = path
                .strip_prefix(root)
                .unwrap()
                .to_string_lossy()
                .replace('\\', "/");
            out.insert(rel, fs::read(&path).unwrap());
        }
    }
}

/// Concatenate a genome's contigs, so that two genomes can be compared base by
/// base without caring how they were cut up.
fn concat(path: &Path) -> Vec<u8> {
    let data = read_genome(path).unwrap();
    let mut out = Vec::with_capacity(data.total_length() as usize);
    for c in &data.contigs {
        out.extend_from_slice(&c.seq);
    }
    out
}

/// First offset of `needle` in `haystack`, anchored on its first bytes so that a
/// genome-sized haystack does not cost quadratic time.
fn find_substring(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    let anchor = &needle[..8.min(needle.len())];
    let mut from = 0usize;
    while from + anchor.len() <= haystack.len() {
        let rel = haystack[from..]
            .windows(anchor.len())
            .position(|w| w == anchor)?;
        let start = from + rel;
        if start + needle.len() <= haystack.len()
            && haystack[start..start + needle.len()] == *needle
        {
            return Some(start);
        }
        from = start + 1;
    }
    None
}

fn occurs(haystack: &[u8], needle: &[u8]) -> bool {
    find_substring(haystack, needle).is_some()
}

fn occurs_either_strand(haystack: &[u8], needle: &[u8]) -> bool {
    occurs(haystack, needle) || occurs(haystack, &revcomp(needle))
}

fn mismatches(a: &[u8], b: &[u8]) -> usize {
    assert_eq!(a.len(), b.len(), "comparing sequences of different length");
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

/// The reference interval an event touches: its own base for a substitution or an
/// insertion anchor, the whole removed stretch for a deletion.
fn touched(m: &infinity_sim::Mutation) -> (usize, usize) {
    match m.kind {
        MutationKind::Snp => (m.pos, m.pos + 1),
        MutationKind::Insertion | MutationKind::Gain => (m.pos, m.pos + 1),
        MutationKind::Deletion | MutationKind::Loss => (m.pos, m.pos + m.ref_bases.len()),
    }
}

fn small(seed: u64) -> SimConfig {
    SimConfig {
        seed,
        reference: RefSpec::Random {
            length: 100_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 8,
        mean_snps_per_branch: 3.0,
        random_orientation: false,
        ..SimConfig::default()
    }
}

// ---------------------------------------------------------------------------
// 1. Determinism
// ---------------------------------------------------------------------------

#[test]
fn same_seed_is_byte_identical_and_a_different_seed_is_not() {
    let cfg = SimConfig {
        seed: 4242,
        reference: RefSpec::Random {
            length: 60_000,
            gc: 0.42,
        },
        clusters: 2,
        leaves_per_cluster: 5,
        divergence: 0.02,
        mean_snps_per_branch: 4.0,
        mean_indels_per_branch: 0.5,
        segment_events: SegmentEvents {
            mean_gains_per_branch: 0.3,
            mean_losses_per_branch: 0.3,
            min_len: 500,
            max_len: 2000,
        },
        contigs: 4,
        random_orientation: true,
        ..SimConfig::default()
    };

    let a_dir = tmpdir("det-a");
    let b_dir = tmpdir("det-b");
    let c_dir = tmpdir("det-c");
    let a = simulate(&cfg, &a_dir).unwrap();
    let b = simulate(&cfg, &b_dir).unwrap();

    // Every file, byte for byte: FASTA, truth, manifest.
    let fa = all_files(&a_dir);
    let fb = all_files(&b_dir);
    assert_eq!(fa.keys().collect::<Vec<_>>(), fb.keys().collect::<Vec<_>>());
    assert!(!fa.is_empty());
    for (name, bytes) in &fa {
        assert_eq!(
            bytes, &fb[name],
            "{name} differs between two runs of the same seed"
        );
    }
    assert_eq!(a.truth, b.truth);
    assert_eq!(a.genomes.len(), 10);
    // The output is not empty either: 10 genomes, 2 roots, a reference, truth.
    assert!(fa.keys().any(|k| k.starts_with("genomes/")));
    assert!(fa.keys().any(|k| k.starts_with("truth/")));
    assert!(fa.contains_key("reference.fasta"));

    // A different seed changes the data but not the shape of the output.
    let other = SimConfig {
        seed: 4243,
        ..cfg.clone()
    };
    let c = simulate(&other, &c_dir).unwrap();
    let fc = all_files(&c_dir);
    assert_eq!(fc.keys().collect::<Vec<_>>(), fa.keys().collect::<Vec<_>>());
    let differing = fa
        .iter()
        .filter(|(name, bytes)| fc.get(*name).map(|b| b != *bytes).unwrap_or(true))
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    assert!(
        differing.iter().any(|n| n.starts_with("genomes/")),
        "a new seed must change the genomes, changed: {differing:?}"
    );
    assert_ne!(a.truth, c.truth);
    assert_ne!(a.truth.pairwise, c.truth.pairwise);

    for d in [&a_dir, &b_dir, &c_dir] {
        fs::remove_dir_all(d).ok();
    }
}

#[test]
fn rerunning_into_the_same_directory_is_stable() {
    let cfg = small(99);
    let dir = tmpdir("det-same-dir");
    let first = simulate(&cfg, &dir).unwrap();
    let before = all_files(&dir);
    let second = simulate(&cfg, &dir).unwrap();
    let after = all_files(&dir);
    assert_eq!(before, after);
    assert_eq!(first.truth, second.truth);
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// 2. Mutation exclusivity
// ---------------------------------------------------------------------------

#[test]
fn no_two_events_in_a_cluster_touch_the_same_position() {
    let cfg = SimConfig {
        seed: 31337,
        reference: RefSpec::Random {
            length: 120_000,
            gc: 0.5,
        },
        clusters: 2,
        leaves_per_cluster: 9,
        divergence: 0.02,
        mean_snps_per_branch: 4.0,
        mean_indels_per_branch: 0.8,
        segment_events: SegmentEvents {
            mean_gains_per_branch: 0.4,
            mean_losses_per_branch: 0.4,
            min_len: 500,
            max_len: 3000,
        },
        ..SimConfig::default()
    };
    let dir = tmpdir("exclusive");
    let out = simulate(&cfg, &dir).unwrap();
    let truth = &out.truth;

    for &cluster in &truth.cluster_ids() {
        let members = truth.members(cluster);
        let root = truth.root(cluster).unwrap();

        // The root's own substitutions are all at distinct positions.
        let root_positions: BTreeSet<usize> = root.mutations.iter().map(|m| m.pos).collect();
        assert_eq!(root_positions.len(), root.mutations.len());
        assert!(root.mutations.iter().all(|m| m.kind == MutationKind::Snp));

        // Every distinct event of this cluster, whichever leaves inherited it.
        // An event on a shared ancestral branch is in every descendant's log, so
        // the same event legitimately appears several times; what must never
        // happen is two *different* events touching one position.
        let mut distinct: BTreeMap<EventKey, EventSeen> = BTreeMap::new();
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut snps_and_anchors: BTreeSet<usize> = BTreeSet::new();
        for name in &members {
            let log = truth.mutations_of(name).unwrap();
            for m in log {
                let (s, e) = touched(m);
                assert!(s < e, "{name}: empty event range");
                match distinct.entry((m.pos, m.kind)) {
                    btree_map::Entry::Occupied(seen) => {
                        let (se, sref, salt) = seen.get();
                        assert_eq!(*se, e, "two different events start at {}", m.pos);
                        assert_eq!(
                            sref, &m.ref_bases,
                            "the same event is recorded twice, differently"
                        );
                        assert_eq!(salt, &m.alt_bases);
                    }
                    btree_map::Entry::Vacant(slot) => {
                        slot.insert((e, m.ref_bases.clone(), m.alt_bases.clone()));
                        ranges.push((s, e));
                    }
                }
                if !m.kind.is_deletion() {
                    snps_and_anchors.insert(m.pos);
                }
            }
            // Within one lineage no position is touched twice, which is the
            // property section 10.1 states outright.
            let positions: BTreeSet<usize> = log.iter().map(|m| m.pos).collect();
            assert_eq!(
                positions.len(),
                log.len(),
                "{name}: two events share a position"
            );
            assert!(
                log.windows(2).all(|w| w[0].pos <= w[1].pos),
                "{name}: the log is not sorted by position"
            );
        }
        // Across the whole cluster: sort by start and check nothing overlaps.
        ranges.sort_unstable();
        for w in ranges.windows(2) {
            assert!(
                w[0].1 <= w[1].0,
                "cluster {cluster}: events [{},{}) and [{},{}) overlap",
                w[0].0,
                w[0].1,
                w[1].0,
                w[1].1
            );
        }
        // A branch event never lands on a position the root already substituted.
        // (A segment loss may *cover* one: those substitutions are shared by the
        // whole cluster and cancel out of every within-cluster distance.)
        for &p in &snps_and_anchors {
            assert!(
                !root_positions.contains(&p),
                "cluster {cluster}: branch event at {p} is also a root substitution"
            );
        }
    }

    // The payoff: two genomes differ at exactly the symmetric difference of
    // their private SNP positions, so that count is the path length in the tree.
    for &cluster in &truth.cluster_ids() {
        let members = truth.members(cluster);
        for (i, a) in members.iter().enumerate() {
            let sa: BTreeSet<usize> = truth
                .mutations_of(a)
                .unwrap()
                .iter()
                .filter(|m| m.kind == MutationKind::Snp)
                .map(|m| m.pos)
                .collect();
            for b in &members[i + 1..] {
                let sb: BTreeSet<usize> = truth
                    .mutations_of(b)
                    .unwrap()
                    .iter()
                    .filter(|m| m.kind == MutationKind::Snp)
                    .map(|m| m.pos)
                    .collect();
                let symdiff = sa.union(&sb).count() - sa.intersection(&sb).count();
                assert_eq!(
                    symdiff as u64,
                    truth.pairwise_snps(a, b).unwrap(),
                    "symmetric difference is not the truth distance for {a} / {b}"
                );
            }
        }
    }

    fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_single_snp_per_branch_lands_where_the_log_says() {
    // One SNP per branch, no noise: the log must describe the sequence exactly.
    let cfg = SimConfig {
        seed: 8,
        reference: RefSpec::Random {
            length: 40_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 5,
        divergence: 0.0,
        mean_snps_per_branch: 1.0,
        random_orientation: false,
        ..SimConfig::default()
    };
    let dir = tmpdir("one-snp");
    let out = simulate(&cfg, &dir).unwrap();
    let reference = concat(&out.reference);
    for path in &out.genomes {
        let name = infinity_core::fasta::genome_name(path);
        let seq = concat(path);
        assert_eq!(seq.len(), reference.len());
        let log = out.truth.mutations_of(&name).unwrap();
        let mut expected = reference.clone();
        for m in log {
            assert_eq!(m.kind, MutationKind::Snp);
            assert_eq!(m.ref_bases, vec![reference[m.pos]]);
            assert_ne!(m.alt_bases, m.ref_bases);
            expected[m.pos] = m.alt_bases[0];
        }
        assert_eq!(
            seq, expected,
            "{name}: the log does not reproduce the genome"
        );
    }
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// 3. Truth self-consistency
// ---------------------------------------------------------------------------

#[test]
fn pairwise_truth_equals_the_path_length_in_the_truth_tree() {
    let cfg = SimConfig {
        seed: 2718,
        reference: RefSpec::Random {
            length: 80_000,
            gc: 0.5,
        },
        clusters: 3,
        leaves_per_cluster: 11,
        mean_snps_per_branch: 5.0,
        ..SimConfig::default()
    };
    let dir = tmpdir("selfconsistent");
    let out = simulate(&cfg, &dir).unwrap();
    let truth = &out.truth;

    let mut rows_from_trees = Vec::new();
    for &cluster in &truth.cluster_ids() {
        let newick = truth.tree(cluster).unwrap();

        // The engine's own Newick reader must accept it, and see a binary tree
        // over exactly this cluster's leaves with integer, non-negative lengths.
        let tree = parse_newick(newick).unwrap();
        let mut leaves = tree.leaf_names();
        leaves.sort();
        assert_eq!(leaves, truth.members(cluster));
        assert_eq!(tree.leaf_count(), cfg.leaves_per_cluster);
        check_lengths_are_snp_counts(&tree);

        // And the path sums recomputed from the file are the published matrix.
        rows_from_trees.extend(pairwise_from_newick(newick).unwrap());
    }
    rows_from_trees.sort();
    let mut published = truth.pairwise.clone();
    published.sort();
    assert_eq!(rows_from_trees, published);

    // The expected number of rows: every within-cluster pair, nothing else.
    let per_cluster = cfg.leaves_per_cluster * (cfg.leaves_per_cluster - 1) / 2;
    assert_eq!(published.len(), cfg.clusters * per_cluster);
    for (a, b, _) in &published {
        assert!(a < b);
        assert_eq!(truth.cluster_of(a), truth.cluster_of(b));
    }
    // Across clusters there is deliberately no truth value.
    assert!(truth
        .pairwise_snps(&truth.members(0)[0], &truth.members(1)[0])
        .is_none());
    // A genome against itself is zero by definition.
    assert_eq!(
        truth.pairwise_snps(&published[0].0, &published[0].0),
        Some(0)
    );

    fs::remove_dir_all(&dir).ok();
}

/// Branch lengths are SNP counts, so every one of them is a non-negative integer
/// and the root carries none.
fn check_lengths_are_snp_counts(tree: &infinity_core::tree::Tree) {
    let mut stack = vec![(tree, true)];
    while let Some((node, is_root)) = stack.pop() {
        let len = node.len();
        assert!(len >= 0.0 && len.is_finite(), "branch length {len}");
        assert!(
            (len - len.round()).abs() < 1e-9,
            "branch length {len} is not a SNP count"
        );
        if is_root {
            assert_eq!(len, 0.0, "the cluster root must not carry a branch length");
        }
        for c in node.children() {
            stack.push((c, false));
        }
    }
}

#[test]
fn observed_differences_match_the_truth_distance_exactly_without_indels() {
    // T1's premise: with substitutions only, the number of positions at which
    // two genomes actually differ is the truth number, not an approximation.
    let cfg = SimConfig {
        seed: 555,
        reference: RefSpec::Random {
            length: 150_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 12,
        mean_snps_per_branch: 3.0,
        random_orientation: false,
        ..SimConfig::default()
    };
    let dir = tmpdir("t1-premise");
    let out = simulate(&cfg, &dir).unwrap();
    let mut seqs: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    for path in &out.genomes {
        let name = infinity_core::fasta::genome_name(path);
        seqs.insert(name, concat(path));
    }
    let names = out.truth.genome_names();
    let mut checked = 0usize;
    for (i, a) in names.iter().enumerate() {
        for b in &names[i + 1..] {
            let observed = mismatches(&seqs[a], &seqs[b]) as u64;
            let want = out.truth.pairwise_snps(a, b).unwrap();
            assert_eq!(observed, want, "{a} vs {b}");
            checked += 1;
        }
    }
    assert_eq!(checked, 66);
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// 4. Reference round-trip
// ---------------------------------------------------------------------------

#[test]
fn a_supplied_reference_gives_roots_at_the_requested_divergence() {
    let dir = tmpdir("refroundtrip");
    let ref_path = dir.join("supplied.fasta");
    // A reference with an N run and lowercase input, to check normalisation too.
    let mut ref_seq: Vec<u8> = (0..50_000)
        .map(|i| b"acgtACGT"[(i * 7 + i / 13) % 8])
        .collect();
    ref_seq[10_000..10_200].fill(b'n');
    write_single_fasta(&ref_path, "supplied", &ref_seq);

    let out_dir = dir.join("out");
    let divergence = 0.02;
    let cfg = SimConfig {
        seed: 1234,
        reference: RefSpec::File(ref_path.clone()),
        clusters: 3,
        leaves_per_cluster: 4,
        divergence,
        mean_snps_per_branch: 2.0,
        random_orientation: false,
        ..SimConfig::default()
    };
    let out = simulate(&cfg, &out_dir).unwrap();
    let truth = &out.truth;

    assert_eq!(truth.reference.source, RefSource::File);
    assert_eq!(truth.reference.length, 50_000);
    let reference = concat(&out.reference);
    assert_eq!(reference.len(), ref_seq.len());
    // The reference is normalised on the way in, exactly as the engine reads it.
    assert!(reference
        .iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')));
    assert!(reference[10_000..10_200].iter().all(|&b| b == b'N'));

    // The N run cannot carry a substitution, so the achievable count is the
    // requested fraction of the whole length, capped by what is substitutable.
    let requested = (divergence * reference.len() as f64).round() as usize;
    for root_path in &out.roots {
        let root = concat(root_path);
        assert_eq!(
            root.len(),
            reference.len(),
            "a root only substitutes, never indels"
        );
        let differences = mismatches(&reference, &root);
        assert_eq!(differences, requested);
        // No substitution sits on an N.
        for (i, (r, s)) in reference.iter().zip(root.iter()).enumerate() {
            if r != s {
                assert_ne!(*r, b'N', "position {i} was an N");
            }
        }
    }
    for r in &truth.roots {
        assert_eq!(r.substitutions as usize, requested);
        assert_eq!(r.length, reference.len() as u64);
        // Tolerance: the count is an integer, so the realised divergence can be
        // at most half a substitution away from the requested fraction.
        let realised = r.substitutions as f64 / r.length as f64;
        assert!(
            (realised - divergence).abs() <= 0.5 / r.length as f64,
            "realised {realised} vs requested {divergence}"
        );
    }

    // With no indels and no segment events every leaf keeps the reference length.
    for g in &truth.genomes {
        assert_eq!(g.length, reference.len() as u64, "{}", g.name);
        assert_eq!(g.contigs, 1);
    }
    for path in &out.genomes {
        assert_eq!(concat(path).len(), reference.len());
    }
    fs::remove_dir_all(&dir).ok();
}

fn write_single_fasta(path: &Path, name: &str, seq: &[u8]) {
    let mut text = format!(">{name}\n");
    for line in seq.chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    fs::write(path, text).unwrap();
}

// ---------------------------------------------------------------------------
// 5. Indels
// ---------------------------------------------------------------------------

#[test]
fn indel_lengths_account_for_every_base_of_length_change() {
    let cfg = SimConfig {
        seed: 606,
        reference: RefSpec::Random {
            length: 100_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 12,
        mean_snps_per_branch: 2.0,
        mean_indels_per_branch: 1.5,
        indel_min_len: 1,
        indel_max_len: 10,
        random_orientation: false,
        ..SimConfig::default()
    };
    let dir = tmpdir("indels");
    let out = simulate(&cfg, &dir).unwrap();
    let truth = &out.truth;
    let root_len = truth.root(0).unwrap().length as i64;
    let root_seq = concat(&out.roots[0]);
    assert_eq!(root_seq.len() as i64, root_len);

    let mut insertions = 0usize;
    let mut deletions = 0usize;
    for g in &truth.genomes {
        let mut delta = 0i64;
        for m in &g.mutations {
            match m.kind {
                MutationKind::Insertion => {
                    insertions += 1;
                    assert!(m.ref_bases.is_empty());
                    assert!(
                        (1..=10).contains(&m.alt_bases.len()),
                        "insertion of {}",
                        m.alt_bases.len()
                    );
                }
                MutationKind::Deletion => {
                    deletions += 1;
                    assert!(m.alt_bases.is_empty());
                    assert!(
                        (1..=10).contains(&m.ref_bases.len()),
                        "deletion of {}",
                        m.ref_bases.len()
                    );
                }
                MutationKind::Snp => {
                    assert_eq!(m.ref_bases.len(), 1);
                    assert_eq!(m.alt_bases.len(), 1);
                }
                other => panic!("unexpected event {other:?} with indels only"),
            }
            delta += m.length_delta();
        }
        assert_eq!(
            g.length as i64,
            root_len + delta,
            "{}: length is not the root length plus the recorded indels",
            g.name
        );
        let seq = concat(
            out.genomes
                .iter()
                .find(|p| p.ends_with(format!("{}.fasta", g.name)))
                .unwrap(),
        );
        assert_eq!(
            seq.len() as i64,
            g.length as i64,
            "{}: FASTA and truth disagree",
            g.name
        );
    }
    assert!(
        insertions > 0 && deletions > 0,
        "expected both kinds of indel"
    );

    // The mutation log is a replay script: applying it to the root, left to
    // right in reference coordinates, reproduces the genome.
    for g in &truth.genomes {
        let rebuilt = replay(&root_seq, &g.mutations);
        let path = out
            .genomes
            .iter()
            .find(|p| p.ends_with(format!("{}.fasta", g.name)))
            .unwrap();
        assert_eq!(
            rebuilt,
            concat(path),
            "{}: replay does not reproduce the genome",
            g.name
        );
    }
    fs::remove_dir_all(&dir).ok();
}

/// Apply a mutation log to a sequence, left to right in reference coordinates.
///
/// Independent of the simulator's own machinery on purpose: if the log and the
/// FASTA agree through two different routes, the log really does describe the
/// genome. It only works because events sit at disjoint reference positions, so
/// at most one of them can start at any given base.
fn replay(root: &[u8], log: &[infinity_sim::Mutation]) -> Vec<u8> {
    let mut by_pos: BTreeMap<usize, &infinity_sim::Mutation> = BTreeMap::new();
    for m in log {
        assert!(
            by_pos.insert(m.pos, m).is_none(),
            "two events start at {}",
            m.pos
        );
    }
    let mut out: Vec<u8> = Vec::with_capacity(root.len());
    let mut pos = 0usize;
    while pos < root.len() {
        match by_pos.get(&pos) {
            None => {
                out.push(root[pos]);
                pos += 1;
            }
            Some(m) => match m.kind {
                MutationKind::Snp => {
                    assert_eq!(
                        m.ref_bases,
                        vec![root[pos]],
                        "the log's ref base is wrong at {pos}"
                    );
                    out.push(m.alt_bases[0]);
                    pos += 1;
                }
                // An insertion sits in front of the base it is anchored to, which
                // stays where it is.
                MutationKind::Insertion | MutationKind::Gain => {
                    out.extend_from_slice(&m.alt_bases);
                    out.push(root[pos]);
                    pos += 1;
                }
                MutationKind::Deletion | MutationKind::Loss => {
                    pos += m.ref_bases.len();
                }
            },
        }
    }
    out
}

// ---------------------------------------------------------------------------
// 6. Segment gain and loss
// ---------------------------------------------------------------------------

#[test]
fn gained_segments_are_new_and_lost_ones_are_gone() {
    let cfg = SimConfig {
        seed: 909,
        reference: RefSpec::Random {
            length: 300_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 10,
        mean_snps_per_branch: 2.0,
        segment_events: SegmentEvents {
            mean_gains_per_branch: 0.9,
            mean_losses_per_branch: 0.9,
            min_len: 500,
            max_len: 5000,
        },
        random_orientation: false,
        ..SimConfig::default()
    };
    let dir = tmpdir("segments");
    let out = simulate(&cfg, &dir).unwrap();
    let truth = &out.truth;
    let root_seq = concat(&out.roots[0]);
    let root_len = truth.root(0).unwrap().length;

    let mut gains = 0usize;
    let mut losses = 0usize;
    for g in &truth.genomes {
        let path = out
            .genomes
            .iter()
            .find(|p| p.ends_with(format!("{}.fasta", g.name)))
            .unwrap();
        let seq = concat(path);
        let mut delta = 0i64;
        for m in &g.mutations {
            delta += m.length_delta();
            match m.kind {
                MutationKind::Gain => {
                    gains += 1;
                    assert!(
                        (500..=5000).contains(&m.alt_bases.len()),
                        "gain of {}",
                        m.alt_bases.len()
                    );
                    assert!(m.ref_bases.is_empty());
                    assert!(
                        occurs(&seq, &m.alt_bases),
                        "{}: gained segment is missing",
                        g.name
                    );
                    assert!(
                        !occurs(&root_seq, &m.alt_bases),
                        "{}: gained segment was already in the root",
                        g.name
                    );
                }
                MutationKind::Loss => {
                    losses += 1;
                    assert!(
                        (500..=5000).contains(&m.ref_bases.len()),
                        "loss of {}",
                        m.ref_bases.len()
                    );
                    assert!(m.alt_bases.is_empty());
                    assert!(
                        occurs(&root_seq, &m.ref_bases),
                        "{}: lost segment is not in the root",
                        g.name
                    );
                    assert!(
                        !occurs(&seq, &m.ref_bases),
                        "{}: lost segment is still there",
                        g.name
                    );
                    // The loss really is at the recorded reference coordinate.
                    assert_eq!(
                        find_substring(&root_seq, &m.ref_bases),
                        Some(m.pos),
                        "{}: loss recorded at the wrong position",
                        g.name
                    );
                }
                MutationKind::Snp | MutationKind::Insertion | MutationKind::Deletion => {}
            }
        }
        assert_eq!(
            g.length as i64,
            root_len as i64 + delta,
            "{}: length",
            g.name
        );
        assert_eq!(seq.len() as u64, g.length);
    }
    assert!(gains > 0, "no gained segments were drawn");
    assert!(losses > 0, "no lost segments were drawn");
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// 7. Fragmentation
// ---------------------------------------------------------------------------

#[test]
fn fragmentation_preserves_the_sequence_up_to_order_and_orientation() {
    let contigs = 25;
    let base = SimConfig {
        seed: 1717,
        reference: RefSpec::Random {
            length: 250_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 6,
        mean_snps_per_branch: 3.0,
        mean_indels_per_branch: 0.4,
        ..SimConfig::default()
    };
    // Fragmentation has its own derived stream, so the two runs produce the same
    // underlying sequences and only differ in how they are cut.
    let plain_cfg = SimConfig {
        contigs: 1,
        random_orientation: false,
        ..base.clone()
    };
    let cut_cfg = SimConfig {
        contigs,
        random_orientation: true,
        ..base.clone()
    };

    let plain_dir = tmpdir("frag-plain");
    let cut_dir = tmpdir("frag-cut");
    let plain = simulate(&plain_cfg, &plain_dir).unwrap();
    let cut = simulate(&cut_cfg, &cut_dir).unwrap();

    assert_eq!(plain.genomes.len(), cut.genomes.len());
    let mut flipped_total = 0usize;
    for (pp, cp) in plain.genomes.iter().zip(cut.genomes.iter()) {
        assert_eq!(
            infinity_core::fasta::genome_name(pp),
            infinity_core::fasta::genome_name(cp)
        );
        let whole = concat(pp);
        let data = read_genome(cp).unwrap();
        assert_eq!(data.contigs.len(), contigs, "{}: contig count", data.name);
        assert_eq!(
            data.contigs.iter().map(|c| c.seq.len()).sum::<usize>(),
            whole.len(),
            "{}: fragmentation changed the total length",
            data.name
        );
        assert!(data.contigs.iter().all(|c| !c.seq.is_empty()));
        // Contig names are numbered in output order.
        for (i, c) in data.contigs.iter().enumerate() {
            assert_eq!(c.name, format!("{}_ctg{}", data.name, i + 1));
        }
        // Each contig is a piece of the genome, possibly reverse-complemented.
        for c in &data.contigs {
            assert!(
                occurs_either_strand(&whole, &c.seq),
                "{}: contig {} is not a piece of the genome",
                data.name,
                c.name
            );
            if !occurs(&whole, &c.seq) {
                flipped_total += 1;
            }
        }
        // The pieces are disjoint and tile the genome exactly once: sorting them
        // by where they came from must walk from 0 to the end with no gap.
        let mut placed: Vec<(usize, usize)> = data
            .contigs
            .iter()
            .map(|c| {
                let forward = find_substring(&whole, &c.seq);
                let at = match forward {
                    Some(at) => at,
                    None => find_substring(&whole, &revcomp(&c.seq))
                        .expect("a contig that is neither a piece nor its reverse complement"),
                };
                (at, c.seq.len())
            })
            .collect();
        placed.sort_unstable();
        let mut cursor = 0usize;
        for (at, len) in &placed {
            assert_eq!(
                *at, cursor,
                "{}: a gap or an overlap before offset {cursor}",
                data.name
            );
            cursor += len;
        }
        assert_eq!(
            cursor,
            whole.len(),
            "{}: the contigs do not cover the genome",
            data.name
        );
        assert_eq!(placed.iter().map(|&(_, l)| l).sum::<usize>(), whole.len());
    }
    // With 150 contigs and a fair coin, seeing none flipped would be a bug.
    assert!(flipped_total > 0, "no contig was reverse-complemented");

    // Same underlying sequences, so the truth distances are untouched.
    assert_eq!(plain.truth.pairwise, cut.truth.pairwise);
    assert_eq!(plain.truth.trees, cut.truth.trees);
    // Lengths and contig counts in the truth describe what was written.
    for g in &cut.truth.genomes {
        assert_eq!(g.contigs as usize, contigs);
        assert_eq!(g.length, plain.truth.genome(&g.name).unwrap().length);
    }

    fs::remove_dir_all(&plain_dir).ok();
    fs::remove_dir_all(&cut_dir).ok();
}

#[test]
fn one_contig_and_no_orientation_leaves_the_genome_alone() {
    let cfg = SimConfig {
        contigs: 1,
        random_orientation: false,
        ..small(23)
    };
    let dir = tmpdir("frag-none");
    let out = simulate(&cfg, &dir).unwrap();
    for path in &out.genomes {
        let data = read_genome(path).unwrap();
        assert_eq!(data.contigs.len(), 1);
        assert_eq!(data.contigs[0].name, data.name);
    }
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// 8. Small-scale sanity, pinned numbers
// ---------------------------------------------------------------------------

#[test]
fn three_leaves_with_two_snips_per_branch_give_the_pinned_distances() {
    let cfg = SimConfig {
        seed: 11,
        reference: RefSpec::Random {
            length: 50_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 3,
        divergence: 0.0,
        mean_snps_per_branch: 2.0,
        contigs: 1,
        random_orientation: false,
        ..SimConfig::default()
    };
    let dir = tmpdir("pinned");
    let out = simulate(&cfg, &dir).unwrap();
    let truth = &out.truth;

    // Pinned from a run of this exact configuration, then checked by hand, from
    // the tree, from the logs and from the sequences. The topology is
    //
    //            root
    //          /      \
    //      s002:2    (internal):2
    //                /         \
    //            s000:2       s001:3
    //
    // so the path sums are
    //   d(s000,s001) = 2 + 3          = 5
    //   d(s000,s002) = 2 + 2 + 2      = 6
    //   d(s001,s002) = 3 + 2 + 2      = 7
    let newick = truth.tree(0).unwrap();
    assert_eq!(
        newick,
        "(sim_c00_s002:2,(sim_c00_s000:2,sim_c00_s001:3):2);"
    );
    let expected = vec![
        ("sim_c00_s000".to_string(), "sim_c00_s001".to_string(), 5u64),
        ("sim_c00_s000".to_string(), "sim_c00_s002".to_string(), 6),
        ("sim_c00_s001".to_string(), "sim_c00_s002".to_string(), 7),
    ];
    assert_eq!(truth.pairwise, expected);
    assert_eq!(pairwise_from_newick(newick).unwrap(), expected);
    assert_eq!(truth.pairwise_snps("sim_c00_s000", "sim_c00_s001"), Some(5));
    assert_eq!(truth.pairwise_snps("sim_c00_s002", "sim_c00_s000"), Some(6));
    assert_eq!(truth.pairwise_snps("sim_c00_s001", "sim_c00_s002"), Some(7));

    // Each leaf's log is the sum of the branches from the root down to it:
    // s000 and s001 inherit the internal branch's 2 SNPs, s002 hangs off the
    // root directly. Divergence is zero, so the root is the reference itself.
    let reference = concat(&out.reference);
    assert_eq!(concat(&out.roots[0]), reference);
    for (name, own, on_path) in [
        ("sim_c00_s000", 2usize, 4usize),
        ("sim_c00_s001", 3, 5),
        ("sim_c00_s002", 2, 2),
    ] {
        let log = truth.mutations_of(name).unwrap();
        assert_eq!(log.len(), on_path, "{name}: log length");
        assert!(log.iter().all(|m| m.kind == MutationKind::Snp));
        let path = out
            .genomes
            .iter()
            .find(|p| p.ends_with(format!("{name}.fasta")))
            .unwrap();
        let seq = concat(path);
        assert_eq!(seq.len(), reference.len());
        assert_eq!(
            mismatches(&reference, &seq),
            on_path,
            "{name}: observed SNPs vs the root"
        );
        // The log names exactly the positions that differ, with the right bases.
        let mut rebuilt = reference.clone();
        for m in log {
            assert_eq!(
                m.ref_bases,
                vec![reference[m.pos]],
                "{name}: wrong ref base"
            );
            assert_ne!(
                m.alt_bases, m.ref_bases,
                "{name}: a substitution that changes nothing"
            );
            rebuilt[m.pos] = m.alt_bases[0];
        }
        assert_eq!(
            rebuilt, seq,
            "{name}: the log does not reproduce the genome"
        );
        assert!(own <= on_path);
    }

    // The three genomes are all different, and the observed differences reproduce
    // the pinned matrix base for base.
    let seqs: BTreeMap<String, Vec<u8>> = out
        .genomes
        .iter()
        .map(|p| (infinity_core::fasta::genome_name(p), concat(p)))
        .collect();
    assert_eq!(seqs.len(), 3);
    for (a, b, want) in &expected {
        assert_eq!(mismatches(&seqs[a], &seqs[b]) as u64, *want, "{a} vs {b}");
    }
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// 9. Truth round-trip through the files
// ---------------------------------------------------------------------------

#[test]
fn the_truth_read_back_from_disk_is_the_truth_returned() {
    let cfg = SimConfig {
        seed: 4711,
        reference: RefSpec::Random {
            length: 70_000,
            gc: 0.55,
        },
        clusters: 2,
        leaves_per_cluster: 6,
        divergence: 0.03,
        mean_snps_per_branch: 3.0,
        mean_indels_per_branch: 0.4,
        segment_events: SegmentEvents {
            mean_gains_per_branch: 0.2,
            mean_losses_per_branch: 0.2,
            min_len: 500,
            max_len: 1500,
        },
        contigs: 3,
        gc: Some(0.4),
        ..SimConfig::default()
    };
    let dir = tmpdir("roundtrip");
    let out = simulate(&cfg, &dir).unwrap();

    let back = Truth::read_from(&dir).unwrap();
    assert_eq!(back, out.truth);
    assert_eq!(back.config, cfg);
    assert_eq!(back.config.seed, 4711);
    assert_eq!(back.reference.gc, out.truth.reference.gc);
    assert_eq!(back.genomes.len(), 12);
    assert_eq!(back.roots.len(), 2);
    assert_eq!(back.trees.len(), 2);

    // Writing what was read gives the same bytes again.
    let again_dir = tmpdir("roundtrip-again");
    back.write_to(&again_dir).unwrap();
    let first = all_files(&dir.join("truth"));
    let second = all_files(&again_dir.join("truth"));
    assert_eq!(first, second);

    // The files a script would read are all there and well formed.
    let manifest = fs::read_to_string(dir.join("truth/manifest.json")).unwrap();
    assert!(manifest.contains("\"format_version\": 1"));
    assert!(manifest.contains("\"seed\": 4711"));
    let clusters_tsv = fs::read_to_string(dir.join("truth/clusters.tsv")).unwrap();
    assert_eq!(clusters_tsv.lines().count(), 13);
    assert!(clusters_tsv.starts_with("genome\tcluster_id\tfile\tlength\tcontigs\n"));
    let pairs = fs::read_to_string(dir.join("truth/pairwise_snps.tsv")).unwrap();
    assert!(pairs.starts_with("genome_a\tgenome_b\tsnps\n"));
    assert_eq!(pairs.lines().count(), 1 + 2 * 15);
    assert!(pairs.ends_with('\n') && !pairs.ends_with("\n\n"));
    for cluster in [0, 1] {
        assert!(dir.join(format!("truth/tree_c{cluster:02}.nwk")).is_file());
    }
    for name in back.genome_names() {
        assert!(dir.join(format!("truth/mutations_{name}.tsv")).is_file());
    }
    assert!(dir.join("truth/mutations_root_c00.tsv").is_file());

    for d in [&dir, &again_dir] {
        fs::remove_dir_all(d).ok();
    }
}

// ---------------------------------------------------------------------------
// 10. Property 4 of section 10.2, from the simulator side
// ---------------------------------------------------------------------------

#[test]
fn a_simulated_genome_gives_the_engine_the_expected_anchor_density() {
    let length = 2_000_000;
    let cfg = SimConfig {
        seed: 10_000,
        reference: RefSpec::Random { length, gc: 0.5 },
        clusters: 1,
        leaves_per_cluster: 2,
        mean_snps_per_branch: 3.0,
        contigs: 1,
        random_orientation: false,
        ..SimConfig::default()
    };
    let dir = tmpdir("anchors");
    let out = simulate(&cfg, &dir).unwrap();

    let params = Parameters::default();
    assert!(params.validate().is_ok());
    // 2 Mb is above the small-genome threshold, so the default densities stand.
    assert!((length as u64) > infinity_core::params::SMALL_GENOME_LENGTH);

    for path in &out.genomes {
        let data = read_genome(path).unwrap();
        let scan = infinity_core::genome::scan_data(&data, &params);
        assert_eq!(scan.stats.length, length as u64);
        assert_eq!(scan.stats.n_bases, 0, "a simulated genome has no N");
        assert_eq!(scan.stats.contigs, 1);

        // Property 4: the number of level-3 anchors is within 10% of
        // length / densities[2]. Every level-3 anchor is in `scan.anchors`.
        let expected = length as f64 / params.densities[2] as f64;
        let got = scan.anchors.len() as f64;
        assert!(
            (got - expected).abs() <= 0.10 * expected,
            "{}: {} anchors, expected about {expected:.0}",
            data.name,
            scan.anchors.len()
        );
        assert!(scan.anchors.iter().all(|a| a.level <= 3 && a.level >= 1));
        // Anchors are ordered by position within their contig, as 6.2 requires.
        assert!(scan.anchors.windows(2).all(|w| w[0].pos <= w[1].pos));

        // The sparsest level scales the same way; it is far fewer anchors, so
        // the same 10% band would be too tight to be a fair test there.
        let level1 = scan.anchors.iter().filter(|a| a.level == 1).count();
        assert!(level1 > 0 && level1 < scan.anchors.len());

        // And a reverse complement of the same genome gives the same anchor
        // hashes, which is property 1 and the reason fragmentation is safe.
        let flipped = infinity_core::fasta::GenomeData {
            name: data.name.clone(),
            contigs: data
                .contigs
                .iter()
                .map(|c| infinity_core::fasta::Contig {
                    name: c.name.clone(),
                    seq: revcomp(&c.seq),
                })
                .collect(),
        };
        let scan2 = infinity_core::genome::scan_data(&flipped, &params);
        assert_eq!(scan2.anchors.len(), scan.anchors.len());
        let mut h1: Vec<u64> = scan.anchors.iter().map(|a| a.h).collect();
        let mut h2: Vec<u64> = scan2.anchors.iter().map(|a| a.h).collect();
        h1.sort_unstable();
        h2.sort_unstable();
        assert_eq!(
            h1, h2,
            "a genome and its reverse complement must sketch alike"
        );
    }
    fs::remove_dir_all(&dir).ok();
}

// ---------------------------------------------------------------------------
// Shapes the CLI will be asked for
// ---------------------------------------------------------------------------

#[test]
fn five_clusters_of_twenty_are_separable_and_labelled() {
    // T5's shape, scaled down so the test stays quick: clusters must be far
    // enough apart that no two of them could be confused, and every genome must
    // carry the right label.
    let cfg = SimConfig {
        seed: 55,
        reference: RefSpec::Random {
            length: 100_000,
            gc: 0.5,
        },
        clusters: 5,
        leaves_per_cluster: 4,
        divergence: 0.02,
        mean_snps_per_branch: 3.0,
        ..SimConfig::default()
    };
    let dir = tmpdir("t5-shape");
    let out = simulate(&cfg, &dir).unwrap();
    let truth = &out.truth;
    assert_eq!(truth.genomes.len(), 20);
    assert_eq!(truth.roots.len(), 5);
    assert_eq!(truth.cluster_ids(), vec![0, 1, 2, 3, 4]);
    for &c in &truth.cluster_ids() {
        assert_eq!(truth.members(c).len(), 4);
    }
    // Within a cluster genomes are a handful of SNPs apart; between roots they
    // are 2% apart, four orders of magnitude more.
    let roots: Vec<Vec<u8>> = out.roots.iter().map(|p| concat(p)).collect();
    for i in 0..roots.len() {
        for j in (i + 1)..roots.len() {
            let d = mismatches(&roots[i], &roots[j]);
            assert!(
                (3500..4500).contains(&d),
                "roots {i} and {j} are {d} apart, expected about 4000"
            );
        }
    }
    let mut max_within = 0usize;
    for &c in &truth.cluster_ids() {
        let members = truth.members(c);
        for (i, a) in members.iter().enumerate() {
            for b in &members[i + 1..] {
                let pa = out
                    .genomes
                    .iter()
                    .find(|p| p.ends_with(format!("{a}.fasta")))
                    .unwrap();
                let pb = out
                    .genomes
                    .iter()
                    .find(|p| p.ends_with(format!("{b}.fasta")))
                    .unwrap();
                max_within = max_within.max(mismatches(&concat(pa), &concat(pb)));
            }
        }
    }
    assert!(
        max_within < 200,
        "the widest within-cluster distance is {max_within} SNPs; it must stay an order of \
         magnitude below the ~4000 that separate two roots, or T5 could not ask for an exact \
         cluster recovery"
    );
    // Names encode the cluster, and the truth agrees with the names.
    for g in &truth.genomes {
        assert!(
            g.name.starts_with(&format!("sim_c{:02}_s", g.cluster)),
            "{}",
            g.name
        );
        assert_eq!(truth.cluster_of(&g.name), Some(g.cluster));
    }
    fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_single_leaf_cluster_and_a_single_snp_mean_are_not_special_cases() {
    let cfg = SimConfig {
        seed: 1,
        reference: RefSpec::Random {
            length: 20_000,
            gc: 0.5,
        },
        clusters: 1,
        leaves_per_cluster: 1,
        mean_snps_per_branch: 0.0,
        divergence: 0.01,
        ..SimConfig::default()
    };
    let dir = tmpdir("degenerate");
    let out = simulate(&cfg, &dir).unwrap();
    assert_eq!(out.genomes.len(), 1);
    assert!(out.truth.pairwise.is_empty());
    assert_eq!(out.truth.tree(0).unwrap(), "sim_c00_s000;");
    assert!(out.truth.mutations_of("sim_c00_s000").unwrap().is_empty());
    // The genome is the root, and the root is 1% away from the reference.
    let reference = concat(&out.reference);
    let root = concat(&out.roots[0]);
    assert_eq!(mismatches(&reference, &root), 200);
    let leaf = concat(&out.genomes[0]);
    assert_eq!(leaf, root);
    fs::remove_dir_all(&dir).ok();
}
