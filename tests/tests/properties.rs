//! The properties of section 10.2 that must always hold.
//!
//! Each test below is numbered after the property it covers. They are
//! integration tests: they drive the real engine over real files, so a
//! regression anywhere in the pipeline shows up here.

use std::collections::HashSet;
use std::path::PathBuf;

use infinity_core::blocks::{block_ranges, hash_block};
use infinity_core::engine::RunResults;
use infinity_core::params::Parameters;
use infinity_core::seq::{contig_anchors, Thresholds, LEVELS};
use infinity_tests::{
    diff_output_dirs, random_dna, reported_pairs, revcomp as revcomp_str, run_engine, with_snps,
    write_genomes, TempDir,
};

/// Parameters for a small, fast test genome. Densities are denser than the
/// defaults so that a few hundred kilobases still produce plenty of anchors;
/// the genome is kept above 500 kb so the small-genome rule of 6.1 does not
/// override them behind the test's back.
fn test_params() -> Parameters {
    Parameters {
        k: 31,
        densities: [1024, 128, 16],
        threads: 4,
        ..Default::default()
    }
}

const GENOME_LEN: usize = 600_000;

fn one_contig_genome(seed: u64) -> String {
    random_dna(GENOME_LEN, seed)
}

fn split_into(seq: &str, n: usize) -> Vec<String> {
    let step = seq.len() / n;
    let mut out: Vec<String> = (0..n)
        .map(|i| seq[i * step..(i + 1) * step].to_string())
        .collect();
    // The tail belongs to the last contig.
    if let Some(last) = out.last_mut() {
        last.push_str(&seq[n * step..]);
    }
    out
}

// ---------------------------------------------------------------------------
// Property 1
// ---------------------------------------------------------------------------

/// A genome and its reverse complement give the same set of anchor hashes and
/// the same set of block hashes.
#[test]
fn property_1_reverse_complement_is_invisible() {
    let seq = one_contig_genome(101);
    let rc = revcomp_str(&seq);
    let p = test_params();
    let th = Thresholds::new(&p.densities);

    let fwd: Vec<u8> = seq.bytes().collect();
    let back: Vec<u8> = rc.bytes().collect();

    let a1 = contig_anchors(&fwd, 0, p.k, &th);
    let a2 = contig_anchors(&back, 0, p.k, &th);
    let mut h1: Vec<u64> = a1.iter().map(|a| a.h).collect();
    let mut h2: Vec<u64> = a2.iter().map(|a| a.h).collect();
    h1.sort_unstable();
    h2.sort_unstable();
    assert_eq!(h1, h2, "the anchor hash multiset must not depend on strand");
    assert!(
        !h1.is_empty(),
        "the test genome must actually produce anchors"
    );

    // The same for block hashes, at every level.
    for level in 1..=LEVELS {
        let pick = |anchors: &[infinity_core::seq::Anchor], s: &[u8]| -> HashSet<u64> {
            let pos: Vec<u64> = anchors
                .iter()
                .filter(|a| a.level <= level as u8)
                .map(|a| a.pos)
                .collect();
            block_ranges(&pos, s.len() as u64, p.k as u64)
                .iter()
                .map(|&(s0, e0)| hash_block(&s[s0 as usize..e0 as usize]).0)
                .collect()
        };
        let b1 = pick(&a1, &fwd);
        let b2 = pick(&a2, &back);
        assert_eq!(
            b1, b2,
            "level {level} block hashes must not depend on strand"
        );
    }
}

// ---------------------------------------------------------------------------
// Property 2
// ---------------------------------------------------------------------------

fn assert_no_difference(results: &RunResults, label: &str) {
    let mut checked = 0;
    for cmp in results.comparisons.values() {
        checked += 1;
        assert!(
            cmp.variants.is_empty(),
            "{label}: expected zero variants, got {}: {:?}",
            cmp.variants.len(),
            cmp.variants
                .iter()
                .take(3)
                .map(|v| format!(
                    "{}@{}{}>{}",
                    v.kind.as_str(),
                    v.pivot_contig,
                    v.pivot_start,
                    v.alt_allele
                ))
                .collect::<Vec<_>>()
        );
    }
    assert!(
        checked > 0,
        "{label}: nothing was compared, so nothing was proven"
    );

    // Every block that was looked up must have matched.
    let c = &results.counters;
    assert!(c.blocks_looked_up > 0, "{label}: no blocks were looked up");
    assert_eq!(
        c.total_matched(),
        c.blocks_looked_up,
        "{label}: {} of {} blocks matched; identical sequence must match everywhere",
        c.total_matched(),
        c.blocks_looked_up
    );
    assert_eq!(
        c.breakpoints, 0,
        "{label}: a breakpoint means a lost mapping"
    );
}

#[test]
fn property_2a_a_genome_compared_with_itself_has_no_differences() {
    let dir = TempDir::new("p2a");
    let seq = one_contig_genome(202);
    // "aaa" sorts first and both have the same N50, so it becomes the pivot.
    let paths = write_genomes(
        &dir.join("in"),
        &[
            ("aaa_self", vec![seq.clone()]),
            ("bbb_copy", vec![seq.clone()]),
        ],
    );
    let (results, _) = run_engine(&paths, &test_params(), None);
    assert_eq!(results.records.len(), 2);
    assert_no_difference(&results, "self vs self");
}

#[test]
fn property_2b_reverse_complement_has_no_differences() {
    let dir = TempDir::new("p2b");
    let seq = one_contig_genome(303);
    let paths = write_genomes(
        &dir.join("in"),
        &[
            ("aaa_fwd", vec![seq.clone()]),
            ("bbb_rev", vec![revcomp_str(&seq)]),
        ],
    );
    let (results, _) = run_engine(&paths, &test_params(), None);
    assert_no_difference(&results, "genome vs its reverse complement");
}

#[test]
fn property_2c_shuffled_and_flipped_contigs_have_no_differences() {
    let dir = TempDir::new("p2c");
    let seq = one_contig_genome(404);
    let contigs = split_into(&seq, 8);
    // Permute and flip. The permutation is fixed, not random, so the test is
    // reproducible; what matters is that the order and orientation differ.
    let order = [5usize, 0, 7, 2, 4, 1, 6, 3];
    let flip = [true, false, true, true, false, true, false, true];
    let shuffled: Vec<String> = order
        .iter()
        .zip(flip.iter())
        .map(|(&i, &f)| {
            if f {
                revcomp_str(&contigs[i])
            } else {
                contigs[i].clone()
            }
        })
        .collect();
    let paths = write_genomes(
        &dir.join("in"),
        &[("aaa_ordered", contigs.clone()), ("bbb_shuffled", shuffled)],
    );
    let (results, _) = run_engine(&paths, &test_params(), None);
    assert_no_difference(&results, "shuffled and flipped contigs");
    // Every contig must still have been placed, so nothing is "novel".
    for cmp in results.comparisons.values() {
        assert_eq!(
            cmp.novel_bp, 0,
            "a shuffled genome must not invent novel bp"
        );
    }
}

// ---------------------------------------------------------------------------
// Property 3
// ---------------------------------------------------------------------------

#[test]
fn property_3_one_snp_gives_exactly_one_snp_at_the_right_place() {
    let dir = TempDir::new("p3");
    let seq = one_contig_genome(505);
    let (mutated, changes) = with_snps(&seq, 1, 0);
    assert_eq!(changes.len(), 1);
    let (pos, ref_base, alt_base) = changes[0];

    let paths = write_genomes(
        &dir.join("in"),
        &[("aaa_ref", vec![seq]), ("bbb_mut", vec![mutated])],
    );
    let (results, _) = run_engine(&paths, &test_params(), None);

    let cmp = results
        .comparisons
        .values()
        .next()
        .expect("one member must have been compared");
    assert_eq!(
        cmp.variants.len(),
        1,
        "exactly one variant expected, got {:?}",
        cmp.variants
    );
    let v = &cmp.variants[0];
    assert_eq!(v.kind, infinity_core::variant::VariantKind::Snp);
    // The reference is the pivot ("aaa_ref"), so the pivot coordinate is the
    // position the SNP was introduced at.
    assert_eq!(v.pivot_start as usize, pos, "wrong position");
    assert_eq!(v.ref_allele.chars().next().unwrap(), ref_base);
    assert_eq!(v.alt_allele.chars().next().unwrap(), alt_base);
    assert_eq!(v.pivot_contig, 0);
    assert_eq!(cmp.snps(), 1);
    assert_eq!(cmp.indels(), 0);
}

#[test]
fn property_3b_several_snps_are_all_found_at_their_positions() {
    let dir = TempDir::new("p3b");
    let seq = one_contig_genome(606);
    let (mutated, changes) = with_snps(&seq, 12, 7);
    let paths = write_genomes(
        &dir.join("in"),
        &[("aaa_ref", vec![seq]), ("bbb_mut", vec![mutated])],
    );
    let (results, _) = run_engine(&paths, &test_params(), None);
    let cmp = results.comparisons.values().next().unwrap();
    let snps: Vec<(u64, char, char)> = cmp
        .variants
        .iter()
        .filter(|v| v.kind == infinity_core::variant::VariantKind::Snp)
        .map(|v| {
            (
                v.pivot_start,
                v.ref_allele.chars().next().unwrap(),
                v.alt_allele.chars().next().unwrap(),
            )
        })
        .collect();
    let expected: Vec<(u64, char, char)> = changes
        .iter()
        .map(|(p, r, a)| (*p as u64, *r, *a))
        .collect();
    assert_eq!(
        snps, expected,
        "every SNP must be found, at its position, with its bases"
    );
}

// ---------------------------------------------------------------------------
// Property 4
// ---------------------------------------------------------------------------

#[test]
fn property_4_anchor_counts_match_the_densities() {
    // A long sequence, so that even the sparsest level has enough anchors for a
    // 10% tolerance to be meaningful rather than noise.
    let len = 8_000_000usize;
    let seq: Vec<u8> = random_dna(len, 707).into_bytes();
    let p = test_params();
    let th = Thresholds::new(&p.densities);

    for level in 1..=LEVELS {
        let anchors = contig_anchors(&seq, 0, p.k, &th)
            .into_iter()
            .filter(|a| a.level <= level as u8)
            .count();
        // Windows containing an N are skipped; this sequence has none, so the
        // number of k-mers is len - k + 1.
        let windows = (len - p.k + 1) as f64;
        let expected = windows / p.densities[level - 1] as f64;
        let dev = (anchors as f64 - expected).abs() / expected;
        assert!(
            dev <= 0.10,
            "level {level}: {anchors} anchors, expected about {expected:.0} (within 10%), off by {:.1}%",
            dev * 100.0
        );
    }
}

// ---------------------------------------------------------------------------
// Property 5
// ---------------------------------------------------------------------------

#[test]
fn property_5_output_is_byte_identical_whatever_the_thread_count() {
    let dir = TempDir::new("p5");
    // Several genomes with real differences, so the parallel paths are actually
    // exercised: a run where nothing happens would trivially be identical.
    let base = one_contig_genome(808);
    let mut genomes = vec![("g00_base".to_string(), vec![base.clone()])];
    for i in 1..6u32 {
        let (m, _) = with_snps(&base, (i * 3) as usize, i as usize * 101);
        genomes.push((format!("g{i:02}_mut"), vec![m]));
    }
    let refs: Vec<(&str, Vec<String>)> = genomes
        .iter()
        .map(|(n, c)| (n.as_str(), c.clone()))
        .collect();
    let paths = write_genomes(&dir.join("in"), &refs);

    let out1 = dir.join("out1");
    let out4 = dir.join("out4");
    let p1 = Parameters {
        threads: 1,
        ..test_params()
    };
    let p4 = Parameters {
        threads: 4,
        ..test_params()
    };
    run_engine(&paths, &p1, Some(&out1));
    run_engine(&paths, &p4, Some(&out4));

    let diffs = diff_output_dirs(&out1, &out4);
    assert!(
        diffs.is_empty(),
        "1 thread and 4 threads must produce identical files, but these differ: {diffs:?}"
    );
    // And there must actually have been files to compare.
    let files: Vec<String> = std::fs::read_dir(&out1)
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    assert!(files.contains(&"genomes.tsv".to_string()), "{files:?}");
    assert!(files.contains(&"tree.nwk".to_string()), "{files:?}");
}

// ---------------------------------------------------------------------------
// Property 6
// ---------------------------------------------------------------------------

fn variants_and_distances(paths: &[PathBuf], params: &Parameters) -> (String, String, String) {
    let (results, _) = run_engine(paths, params, None);
    // Variants, in a canonical textual form.
    let mut v = Vec::new();
    let mut ids: Vec<u32> = results.comparisons.keys().copied().collect();
    ids.sort_unstable();
    for id in ids {
        let cmp = &results.comparisons[&id];
        for var in &cmp.variants {
            v.push(format!(
                "{id}:{}:{}:{}:{}:{}",
                var.kind.as_str(),
                var.pivot_contig,
                var.pivot_start,
                var.ref_allele,
                var.alt_allele
            ));
        }
    }
    let mut pairs: Vec<String> = reported_pairs(&results)
        .into_iter()
        .map(|((a, b), s)| format!("{a}|{b}={s}"))
        .collect();
    pairs.sort();
    // Trees in SNP units, collapsed at half a SNP, compared as *splits* rather
    // than as Newick text. Section 11 uses exactly this resolution when it
    // compares trees, and it is the only honest way to compare them here for two
    // reasons: neighbour joining arranges taxa that are zero distance apart
    // arbitrarily, and the order in which children are written inside a polytomy
    // carries no information at all. Two trees with the same splits are the same
    // tree. See D15.
    let mut trees: Vec<String> = results
        .cluster_snp_trees
        .iter()
        .map(|(id, t)| format!("{id}:{}", split_signature(t)))
        .collect();
    trees.sort();
    (v.join("\n"), pairs.join("\n"), trees.join("\n"))
}

/// A canonical, order-independent signature of a tree's topology: its leaf set
/// plus its bipartitions, after collapsing branches shorter than half a SNP.
fn split_signature(t: &infinity_core::tree::Tree) -> String {
    let collapsed = infinity_core::tree::collapse_short(t.clone(), 0.5);
    let mut names = collapsed.leaf_names();
    names.sort();
    let index: std::collections::HashMap<String, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.clone(), i))
        .collect();
    let mut sp: Vec<String> = infinity_core::tree::splits(&collapsed, &index)
        .iter()
        .map(|s| {
            s.iter()
                .map(|i| names[*i].as_str())
                .collect::<Vec<_>>()
                .join("+")
        })
        .collect();
    sp.sort();
    format!("[{}]<{}>", names.join(","), sp.join("|"))
}

/// One `--switch` case of property 6: the flag's name plus the parameter
/// mutation it stands for. The closure is boxed so that the cases can live in
/// one list despite having distinct closure types.
type SwitchCase = (&'static str, Box<dyn Fn(&mut Parameters)>);

#[test]
fn property_6_the_optimisation_switches_do_not_change_the_answer() {
    let dir = TempDir::new("p6");
    let base = one_contig_genome(909);
    let mut genomes = vec![("g00_base".to_string(), vec![base.clone()])];
    for i in 1..6u32 {
        // Different SNP counts and a shared mutation, so that the piece cache
        // and the profile grouping both have something to do.
        let (m, _) = with_snps(&base, (i * 2) as usize, 500);
        genomes.push((format!("g{i:02}_mut"), vec![m]));
    }
    // Two identical mutants: profile grouping must handle them as one.
    genomes.push(("g99_twin".to_string(), genomes[1].1.clone()));
    let refs: Vec<(&str, Vec<String>)> = genomes
        .iter()
        .map(|(n, c)| (n.as_str(), c.clone()))
        .collect();
    let paths = write_genomes(&dir.join("in"), &refs);

    let baseline = variants_and_distances(&paths, &test_params());
    assert!(
        !baseline.0.is_empty(),
        "the test must produce some variants, otherwise it proves nothing"
    );

    let cases: Vec<SwitchCase> = vec![
        (
            "--flat",
            Box::new(|p: &mut Parameters| p.hierarchical = false),
        ),
        (
            "--no-piece-cache",
            Box::new(|p: &mut Parameters| p.piece_cache = false),
        ),
        (
            "--no-profile-grouping",
            Box::new(|p: &mut Parameters| p.profile_grouping = false),
        ),
    ];
    for (label, modify) in cases {
        let mut p = test_params();
        modify(&mut p);
        let got = variants_and_distances(&paths, &p);
        assert_eq!(got.0, baseline.0, "{label} changed the variants");
        assert_eq!(got.1, baseline.1, "{label} changed the distances");
        assert_eq!(got.2, baseline.2, "{label} changed the trees");
        // Collapsing must not be hiding a genome: the leaf set has to survive.
        assert!(
            got.2.contains("g99_twin") && got.2.contains("g00_base"),
            "{label} lost a leaf from the tree: {}",
            got.2
        );
    }
}

#[test]
fn property_6b_block_skip_keeps_variants_identical() {
    // `--no-block-skip` removes block hashing altogether (section 11). The
    // variants it produces must be the same, because a block match and a piece
    // comparison of identical sequence agree by construction.
    //
    // The coverage bookkeeping is NOT identical, and that difference is recorded
    // as D13 in DECISIONS.md: with block skipping on, the covered set is built
    // from matched blocks as well as from compared pieces, and the two are not
    // the same set when a member contig spans a junction between two pivot
    // contigs. This test pins the part of the property that does hold, so that a
    // regression in variant calling cannot hide behind the known difference.
    let dir = TempDir::new("p6b");
    let base = one_contig_genome(111);
    let (mutated, _) = with_snps(&base, 9, 300);
    let paths = write_genomes(
        &dir.join("in"),
        &[("aaa_ref", vec![base]), ("bbb_mut", vec![mutated])],
    );

    let with_skip = variants_and_distances(&paths, &test_params());
    let mut p = test_params();
    p.block_skip = false;
    let without = variants_and_distances(&paths, &p);

    assert_eq!(
        with_skip.0, without.0,
        "block skipping must not change which variants are called"
    );
    assert!(!with_skip.0.is_empty(), "the test found no variants at all");
}
