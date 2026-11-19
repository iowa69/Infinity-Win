//! Divergent assemblies: a genome of "another species" among close relatives is
//! flagged after scanning, written to `outliers.tsv` and `run.json`, and can be
//! left out with `Parameters::exclude`.

use infinity_core::events::Event;
use infinity_core::params::{Parameters, Resolution};
use infinity_tests::{random_dna, run_engine, with_snps, write_genomes, TempDir};

fn params() -> Parameters {
    Parameters {
        densities: [1024, 128, 16],
        threads: 2,
        resolution: Resolution::Sketch,
        ..Default::default()
    }
}

/// Four close relatives and one unrelated genome of the same length and GC.
fn five_genomes(dir: &std::path::Path) -> Vec<std::path::PathBuf> {
    let base = random_dna(600_000, 7);
    let mut genomes: Vec<(String, Vec<String>)> = (0..4)
        .map(|i| {
            let (s, _) = with_snps(&base, 3 + i, 1000 * (i + 1));
            (format!("rel{i}"), vec![s])
        })
        .collect();
    genomes.push(("stranger".to_string(), vec![random_dna(600_000, 99)]));
    let refs: Vec<(&str, Vec<String>)> = genomes
        .iter()
        .map(|(n, c)| (n.as_str(), c.clone()))
        .collect();
    write_genomes(dir, &refs)
}

#[test]
fn an_unrelated_genome_is_flagged_and_can_be_excluded() {
    let tmp = TempDir::new("outliers");
    let paths = five_genomes(&tmp.join("genomes"));
    let out = tmp.join("out");
    let (results, events) = run_engine(&paths, &params(), Some(&out));

    let flagged: Vec<_> = events
        .iter()
        .filter_map(|e| match &e.event {
            Event::Outliers(list) => Some(list.clone()),
            _ => None,
        })
        .collect();
    assert_eq!(flagged.len(), 1, "exactly one Outliers event");
    let list = &flagged[0];
    assert_eq!(list.len(), 1, "{list:?}");
    assert_eq!(list[0].name, "stranger");
    assert!(list[0].nearest_dist > 0.05);
    assert!(list[0].reasons[0].contains("no close relative"));
    // It arrives before any cluster is refined: right after the scan.
    let pos = |name: &str| events.iter().position(|e| e.event_name() == name);
    assert!(pos("outliers") < pos("clustered"));
    assert_eq!(results.outliers.len(), 1);

    let tsv = std::fs::read_to_string(out.join("outliers.tsv")).unwrap();
    let lines: Vec<&str> = tsv.lines().collect();
    assert_eq!(
        lines[0],
        "id\tname\tlength\tgc\tnearest_dist\tnearest_name\treasons"
    );
    assert_eq!(lines.len(), 2);
    assert!(lines[1].contains("\tstranger\t"));
    let run_json = std::fs::read_to_string(out.join("run.json")).unwrap();
    assert!(
        run_json.contains("divergent assembly stranger"),
        "{run_json}"
    );
    assert!(run_json.contains("--exclude stranger"));

    // Run again without it: nothing flagged, four genomes, and the exclusion
    // is in the notes and the parameters.
    let p = Parameters {
        exclude: vec!["stranger".into(), "typo".into()],
        ..params()
    };
    let out2 = tmp.join("out2");
    let (results, events) = run_engine(&paths, &p, Some(&out2));
    assert_eq!(results.records.len(), 4);
    assert!(results.outliers.is_empty());
    assert!(!events.iter().any(|e| e.event_name() == "outliers"));
    assert_eq!(results.excluded, vec!["stranger".to_string()]);
    let run_json = std::fs::read_to_string(out2.join("run.json")).unwrap();
    assert!(run_json.contains("excluded by request (1): stranger"));
    assert!(run_json.contains("match no input genome: typo"));
    let tsv = std::fs::read_to_string(out2.join("outliers.tsv")).unwrap();
    assert_eq!(tsv.lines().count(), 1, "header only");
}

#[test]
fn excluding_everything_is_an_error() {
    let tmp = TempDir::new("outliers-all");
    let paths = five_genomes(&tmp.join("genomes"));
    let p = Parameters {
        exclude: ["rel0", "rel1", "rel2", "rel3", "stranger"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
        ..params()
    };
    let cfg = infinity_core::engine::RunConfig {
        inputs: paths,
        params: p,
        outdir: None,
    };
    let err = infinity_core::engine::run_collecting(cfg).unwrap_err();
    assert!(err.contains("every input genome was excluded"), "{err}");
}
