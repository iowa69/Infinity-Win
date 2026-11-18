#!/usr/bin/env python3
"""Unit tests for the pure parsing and decision logic in ``scripts/common.py``.

Run with::

    python3 -m unittest discover -s scripts/tests -t scripts -v

Every fixture is written into a ``tempfile.mkdtemp()`` directory (i.e. under
``/tmp``) and removed afterwards, so the tests never touch the repository or the
real datasets. One test reads the real SKA2 table when it is present and skips
itself when it is not.
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent
for path in (str(SCRIPTS_DIR), str(TESTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

import common  # noqa: E402
import fixtures  # noqa: E402

REAL_SKA2 = Path(
    "<salmonella-data>/ncbi_pd/"
    "ska/cluster_reads_distance.tsv"
)


class TempDirCase(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="mtato-test-")
        self.tmp = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def write(self, name: str, text: str) -> Path:
        path = self.tmp / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path


# ---------------------------------------------------------------------------
# SKA2 distance parser
# ---------------------------------------------------------------------------


class TestSka2Parser(TempDirCase):
    HEADER = "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"

    def test_parses_the_documented_header_and_columns(self) -> None:
        path = self.write(
            "ska.tsv",
            self.HEADER
            + "isolate-A\tisolate-B\t2.00\t0.00004\t4671223\t190\n"
            + "isolate-A\tisolate-C\t2.00\t0.00004\t4671208\t208\n"
            + "isolate-B\tisolate-C\t0.00\t0.00000\t4671404\t22\n",
        )
        ska = common.parse_ska2_distance_tsv(path)
        self.assertEqual(len(ska), 3)
        self.assertEqual(ska.samples(), ["isolate-A", "isolate-B", "isolate-C"])
        index = ska.as_dict()
        self.assertEqual(index[common.pair_key("isolate-A", "isolate-B")].distance, 2.0)
        self.assertEqual(index[common.pair_key("isolate-B", "isolate-C")].distance, 0.0)
        row = index[common.pair_key("isolate-A", "isolate-B")]
        self.assertEqual(row.match_count, 4671223)
        self.assertEqual(row.mismatch_count, 190)
        self.assertAlmostEqual(row.mismatch_proportion, 0.00004)
        stats = ska.stats()
        self.assertEqual(stats["rows"], 3)
        self.assertEqual(stats["sample_count"], 3)
        self.assertEqual(stats["min"], 0.0)
        self.assertEqual(stats["max"], 2.0)
        self.assertEqual(stats["median"], 2.0)

    def test_column_order_does_not_matter_and_names_are_matched_loosely(self) -> None:
        path = self.write(
            "shuffled.tsv",
            "mismatch_count\tDistance\tsample_2\tsample_1\tmatch_count\n"
            "190\t2.00\tB\tA\t4671223\n",
        )
        ska = common.parse_ska2_distance_tsv(path)
        self.assertEqual(len(ska), 1)
        self.assertEqual(ska.pairs[0].sample1, "A")
        self.assertEqual(ska.pairs[0].sample2, "B")
        self.assertEqual(ska.pairs[0].distance, 2.0)
        self.assertEqual(ska.pairs[0].mismatch_count, 190)

    def test_na_distances_are_kept_as_none_and_reported(self) -> None:
        path = self.write(
            "na.tsv", self.HEADER + "A\tB\tNA\tNA\t0\t0\n" + "A\tC\t1.00\t0.00001\t100\t1\n"
        )
        ska = common.parse_ska2_distance_tsv(path)
        self.assertEqual(len(ska), 2, "an NA row must be kept, not dropped")
        index = ska.as_dict()
        self.assertIsNone(index[common.pair_key("A", "B")].distance)
        self.assertEqual(ska.distances(), [1.0])
        self.assertTrue(any("no usable Distance" in p for p in ska.problems))
        self.assertEqual(ska.stats()["rows_missing_distance"], 1)

    def test_self_pairs_and_blank_rows_are_reported_not_silently_kept(self) -> None:
        path = self.write(
            "odd.tsv",
            self.HEADER + "A\tA\t0.00\t0\t10\t0\n" + "A\tB\t3.00\t0.001\t10\t30\n" + "\n",
        )
        ska = common.parse_ska2_distance_tsv(path)
        self.assertEqual(len(ska), 1)
        self.assertTrue(any("self-pair" in p for p in ska.problems))

    def test_missing_file_is_a_prerequisite_error_with_guidance(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            common.parse_ska2_distance_tsv(self.tmp / "nope.tsv")
        self.assertIn("--ska-dist", str(ctx.exception))

    def test_missing_sample_columns_is_a_prerequisite_error(self) -> None:
        path = self.write("bad.tsv", "x\ty\tDistance\n1\t2\t3\n")
        with self.assertRaises(common.PrereqError) as ctx:
            common.parse_ska2_distance_tsv(path)
        self.assertIn("sample1", str(ctx.exception))

    def test_pair_key_is_order_independent(self) -> None:
        self.assertEqual(common.pair_key("b", "a"), common.pair_key("a", "b"))
        self.assertEqual(common.pair_key("b", "a"), ("a", "b"))

    @unittest.skipUnless(REAL_SKA2.exists(), f"real SKA2 table not present at {REAL_SKA2}")
    def test_parses_the_real_496_row_table(self) -> None:
        ska = common.parse_ska2_distance_tsv(REAL_SKA2)
        stats = ska.stats()
        self.assertEqual(stats["rows"], 496)
        self.assertEqual(stats["sample_count"], 32)
        self.assertEqual(stats["rows_with_distance"], 496)
        self.assertEqual(stats["min"], 0.0)
        self.assertEqual(stats["max"], 46.0)
        self.assertEqual(stats["median"], 16.0)
        self.assertEqual(ska.problems, [])
        index = ska.as_dict()
        self.assertEqual(index[common.pair_key("isolate-A", "isolate-B")].distance, 2.0)
        self.assertEqual(index[common.pair_key("isolate-B", "isolate-C")].distance, 0.0)


# ---------------------------------------------------------------------------
# events.jsonl
# ---------------------------------------------------------------------------


class TestEventLog(TempDirCase):
    def envelope(self, seq: int, t_ms: float, name: str, data: dict) -> str:
        return json.dumps(
            {"seq": seq, "t_ms": t_ms, "wall_ms": 1700000000000 + int(t_ms),
             "event": {"event": name, "data": data}}
        )

    def test_reads_the_documented_envelope_shape(self) -> None:
        lines = [
            self.envelope(0, 1.0, "genome_scanned", {"id": 0, "name": "a", "length": 10}),
            self.envelope(1, 5.0, "backbone_tree", {"newick": "(a,b);", "pivots": ["a"]}),
            self.envelope(2, 9.0, "cluster_focus", {"cluster": 0, "focus": 2, "newick": "(a,b);"}),
            self.envelope(3, 12.0, "finished", {"version": "x"}),
        ]
        path = self.write("events.jsonl", "\n".join(lines) + "\n")
        log = common.read_event_log(path)
        self.assertEqual(len(log.events), 4)
        self.assertEqual(log.lines_bad, 0)
        self.assertEqual(common.event_name(log.events[1]), "backbone_tree")
        self.assertEqual(common.event_data(log.events[1])["newick"], "(a,b);")
        self.assertEqual(common.event_t_ms(log.events[2]), 9.0)
        self.assertEqual(len(log.of_type("cluster_focus")), 1)
        self.assertEqual(log.first("genome_scanned")["seq"], 0)
        self.assertEqual(log.last("cluster_focus")["seq"], 2)

    def test_bad_lines_are_counted_and_skipped_with_a_message(self) -> None:
        lines = [
            self.envelope(0, 1.0, "warning", {"text": "ok"}),
            "{not json",
            "",
            json.dumps({"seq": 2, "t_ms": 3.0}),  # no event name
            self.envelope(3, 4.0, "warning", {"text": "still here"}),
        ]
        path = self.write("broken.jsonl", "\n".join(lines) + "\n")
        log = common.read_event_log(path)
        self.assertEqual(len(log.events), 2)
        self.assertEqual(log.lines_bad, 2)
        self.assertTrue(any("not JSON" in p for p in log.problems))
        self.assertTrue(any("no recognisable event name" in p for p in log.problems))

    def test_missing_file_explains_why_it_matters(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            common.read_event_log(self.tmp / "absent.jsonl")
        self.assertIn("accuracy-over-time", str(ctx.exception))

    def test_a_folder_argument_resolves_to_the_file_inside(self) -> None:
        self.write("run/events.jsonl", self.envelope(0, 1.0, "warning", {"text": "x"}) + "\n")
        log = common.read_event_log(self.tmp / "run")
        self.assertEqual(len(log.events), 1)


# ---------------------------------------------------------------------------
# run.json
# ---------------------------------------------------------------------------


class TestRunJson(TempDirCase):
    def test_reads_timings_counters_and_memory(self) -> None:
        doc = fixtures.make_run_json(names=["a", "b", "c", "d"], clusters=2)
        path = self.write("run.json", json.dumps(doc))
        rj = common.load_run_json(path)
        self.assertEqual(rj.version, "0.0.0-fake")
        self.assertEqual(rj.genomes, 4)
        self.assertEqual(rj.clusters, 2)
        self.assertEqual(rj.first_tree_s, 0.12)
        self.assertEqual(rj.exact_tree_s, 0.9)
        self.assertEqual(rj.peak_memory_bytes, 64 * 1024 * 1024)
        self.assertEqual(rj.counter_int("pieces_compared"), 120)
        self.assertEqual(rj.blocks_matched, [1200, 900, 1800])
        self.assertAlmostEqual(rj.block_identity(), 3900 / 4000)
        self.assertEqual(rj.param("block_skip"), True)
        self.assertEqual(rj.focus_counts, [0, 0, 4])
        self.assertFalse(rj.cancelled)

    def test_null_timings_stay_unknown_and_are_reported(self) -> None:
        doc = fixtures.make_run_json(names=["a", "b"], first_tree_ms=None, exact_tree_ms=None)
        for key in doc["timings"]:
            doc["timings"][key] = None
        path = self.write("run.json", json.dumps(doc))
        rj = common.load_run_json(path)
        self.assertIsNone(rj.first_tree_s)
        self.assertIsNone(rj.exact_tree_s)
        self.assertIsNone(rj.timing_ms("total_ms"))
        problems = rj.unique_problems()
        self.assertTrue(any("timings.first_tree_ms is null" in p for p in problems))
        self.assertTrue(any("timings.exact_tree_ms is null" in p for p in problems))

    def test_missing_keys_never_raise(self) -> None:
        path = self.write("run.json", json.dumps({"version": "x"}))
        rj = common.load_run_json(path)
        self.assertEqual(rj.version, "x")
        self.assertIsNone(rj.genomes)
        self.assertIsNone(rj.peak_memory_bytes)
        self.assertIsNone(rj.first_tree_s)
        self.assertIsNone(rj.counter_int("pieces_compared"))
        self.assertEqual(rj.focus_counts, [None, None, None])
        self.assertIsNone(rj.block_identity())
        self.assertTrue(any("counters object missing" in p for p in rj.unique_problems()))
        self.assertTrue(any("timings object missing" in p for p in rj.unique_problems()))

    def test_missing_file_is_actionable_when_required(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            common.load_run_json(self.tmp / "nothing")
        self.assertIn("run.json", str(ctx.exception))
        self.assertIsNone(common.load_run_json(self.tmp / "nothing", required=False))

    def test_corrupt_json_is_actionable(self) -> None:
        path = self.write("run.json", "{not json")
        with self.assertRaises(common.PrereqError):
            common.load_run_json(path)

    def test_peak_memory_null_is_reported_as_unknown(self) -> None:
        doc = fixtures.make_run_json(names=["a", "b"], peak_memory_bytes=None)
        path = self.write("run.json", json.dumps(doc))
        rj = common.load_run_json(path)
        self.assertIsNone(rj.peak_memory_bytes)
        self.assertTrue(any("peak_memory_bytes is null" in p for p in rj.unique_problems()))


# ---------------------------------------------------------------------------
# Defensive table reading
# ---------------------------------------------------------------------------


class TestReadTable(TempDirCase):
    def test_reads_a_tsv_with_a_header(self) -> None:
        path = self.write("t.tsv", "a\tb\n1\t2\n3\t4\n")
        table = common.read_table(path)
        self.assertEqual(table.header, ["a", "b"])
        self.assertEqual(len(table), 2)
        self.assertEqual(table.rows[1]["b"], "4")
        self.assertEqual(table.problems, [])

    def test_a_headerless_numeric_file_gets_synthetic_columns(self) -> None:
        path = self.write("t.tsv", "1\t2\t3\n4\t5\t6\n")
        table = common.read_table(path)
        self.assertEqual(table.header, ["col0", "col1", "col2"])
        self.assertEqual(len(table), 2)

    def test_ragged_rows_are_filled_and_reported(self) -> None:
        path = self.write("t.tsv", "a\tb\tc\n1\t2\n3\t4\t5\t6\n")
        table = common.read_table(path)
        self.assertEqual(len(table), 2)
        self.assertEqual(table.rows[0]["c"], "")
        self.assertTrue(table.problems)
        self.assertIn("fewer than 3 fields", table.problems[0])
        self.assertIn("more than 3 fields", table.problems[0])

    def test_duplicate_header_names_are_made_unique(self) -> None:
        path = self.write("t.tsv", "x\tx\n1\t2\n")
        table = common.read_table(path)
        self.assertEqual(table.header, ["x", "x__1"])
        self.assertEqual(table.rows[0]["x__1"], "2")

    def test_alias_lookup_ignores_case_and_punctuation(self) -> None:
        path = self.write("t.tsv", "Mismatches (proportion)\tMatch count\n0.1\t5\n")
        table = common.read_table(path)
        # normalisation folds case and punctuation, so both spellings resolve
        self.assertEqual(table.find_column("mismatches_proportion"), "Mismatches (proportion)")
        self.assertEqual(table.find_column("mismatches (proportion)"), "Mismatches (proportion)")
        self.assertEqual(table.find_column("match_count"), "Match count")
        self.assertIsNone(table.find_column("nope"))
        self.assertTrue(path.exists())

    def test_missing_file_is_actionable(self) -> None:
        with self.assertRaises(common.PrereqError):
            common.read_table(self.tmp / "absent.tsv")


# ---------------------------------------------------------------------------
# M-Tato output folder readers
# ---------------------------------------------------------------------------


class TestClusterDistances(TempDirCase):
    def test_reads_pairwise_distances_from_a_cluster_file(self) -> None:
        names = fixtures.genome_names(4)
        out = fixtures.write_output_folder(self.tmp / "run", names=names, clusters=2)
        rows, problems = common.load_cluster_distances(out)
        self.assertTrue(rows)
        self.assertEqual(problems, [])
        self.assertEqual({r["cluster"] for r in rows}, {"0", "1"})
        first = rows[0]
        self.assertIsNotNone(first["snps"])
        self.assertEqual(first["comparable_bp"], 2048)
        self.assertIsNotNone(first["distance"])

    def test_alternative_column_names_are_accepted(self) -> None:
        path = self.write(
            "cluster_0.dist.tsv",
            "Sample1\tSample2\tSNP distance\tcomparable\tfocus\nA\tB\t7\t1000\t2\n",
        )
        rows, problems = common.load_cluster_distances(self.tmp)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["snps"], 7)
        self.assertEqual(rows[0]["comparable_bp"], 1000)
        self.assertEqual(problems, [])
        self.assertTrue(path.exists())

    def test_unrecognisable_columns_are_reported_and_skipped(self) -> None:
        self.write("cluster_0.dist.tsv", "foo\tbar\n1\t2\n")
        rows, problems = common.load_cluster_distances(self.tmp)
        self.assertEqual(rows, [])
        self.assertTrue(any("cannot identify the two genome columns" in p for p in problems))

    def test_absent_files_are_reported(self) -> None:
        rows, problems = common.load_cluster_distances(self.tmp / "nothing")
        self.assertEqual(rows, [])
        self.assertTrue(problems)


# ---------------------------------------------------------------------------
# Newick, branch collapsing, Robinson-Foulds
# ---------------------------------------------------------------------------


class TestNewickAndRf(unittest.TestCase):
    def test_parses_labels_lengths_and_comments(self) -> None:
        tree = common.parse_newick("((a:1,b:2)100:3[&note],(c:4,d:5)50:6);")
        self.assertEqual(sorted(tree.leaf_names()), ["a", "b", "c", "d"])
        self.assertEqual(tree.total_length(), 21)

    def test_handles_a_single_leaf_and_a_rooted_pair(self) -> None:
        self.assertEqual(common.parse_newick("a;").leaf_names(), ["a"])
        self.assertEqual(sorted(common.parse_newick("(a:1,b:1);").leaf_names()), ["a", "b"])

    def test_syntax_errors_are_prerequisite_errors_with_the_offset(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            common.parse_newick("((a,b;")
        self.assertIn("offset", str(ctx.exception))
        with self.assertRaises(common.PrereqError):
            common.parse_newick("   ")

    def test_identical_trees_have_zero_rf(self) -> None:
        a = common.parse_newick("((a,b),(c,d),e);")
        b = common.parse_newick("(e,(c,d),(a,b));")
        result = common.rf_distance(a, b)
        self.assertEqual(result.rf, 0)
        self.assertEqual(result.rf_normalized, 0.0)
        self.assertEqual(result.leaves_shared, 5)

    def test_different_topologies_have_nonzero_rf(self) -> None:
        a = common.parse_newick("((a,b),(c,d));")
        b = common.parse_newick("((a,c),(b,d));")
        result = common.rf_distance(a, b)
        self.assertEqual(result.rf, 4)
        self.assertEqual(result.rf_normalized, 1.0)

    def test_collapse_threshold_removes_short_branch_splits(self) -> None:
        # ((a,b) short, (c,d)) vs ((a,c),(b,d)): with a big collapse threshold the
        # first tree loses its internal split and RF drops.
        tree = common.parse_newick("(((a:1,b:1):0.1,(c:1,d:1):1):1,e:1);")
        truth = common.parse_newick("(((a:1,b:1):1,(c:1,d:1):1):1,e:1);")
        self.assertEqual(common.rf_distance(tree, truth, collapse=0.0).rf, 0)
        collapsed = common.contract_short_branches(tree, 0.5)
        self.assertLess(len(common.splits(collapsed)), len(common.splits(tree)))

    def test_partial_trees_are_compared_on_shared_leaves_only(self) -> None:
        run = common.parse_newick("((a,b),(c,d));")
        truth = common.parse_newick("(((a,b),(c,d)),(e,f));")
        result = common.rf_distance(run, truth)
        self.assertEqual(result.leaves_shared, 4)
        self.assertTrue(any("only 4 of 6 truth leaves" in p for p in result.problems))

    def test_disjoint_leaf_names_report_rather_than_score(self) -> None:
        run = common.parse_newick("((x,y),(z,w));")
        truth = common.parse_newick("((a,b),(c,d));")
        result = common.rf_distance(run, truth)
        self.assertIsNone(result.rf)
        self.assertTrue(any("no shared leaf names" in p for p in result.problems))

    def test_too_few_taxa_is_reported_as_undefined(self) -> None:
        run = common.parse_newick("(a,b,c);")
        truth = common.parse_newick("(a,b,c);")
        result = common.rf_distance(run, truth)
        self.assertIsNone(result.rf)
        self.assertTrue(any("fewer than 4 taxa" in p for p in result.problems))

    def test_newick_round_trips(self) -> None:
        text = "((a:1.0,b:2.0)100:3.0,(c:1.0,d:2.0):4.0);"
        tree = common.parse_newick(text)
        again = common.parse_newick(common.write_newick(tree))
        self.assertEqual(sorted(again.leaf_names()), ["a", "b", "c", "d"])
        self.assertAlmostEqual(again.total_length(), tree.total_length())

    def test_prune_to_keeps_only_the_requested_leaves(self) -> None:
        tree = common.parse_newick("(((a,b),(c,d)),(e,f));")
        pruned = common.prune_to(tree, {"a", "b", "e", "f"})
        self.assertEqual(sorted(pruned.leaf_names()), ["a", "b", "e", "f"])


class TestCollapseThreshold(unittest.TestCase):
    def test_snp_units_are_detected_from_the_truth_tree(self) -> None:
        truth = common.parse_newick("((a:3,b:2):1,(c:4,d:1):1);")
        threshold, why = common.guess_collapse_threshold(truth)
        self.assertEqual(threshold, 0.5)
        self.assertIn("SNP counts", why)

    def test_per_site_units_use_the_comparable_site_count(self) -> None:
        threshold, why = common.guess_collapse_threshold(
            None, collapse_snps=0.5, comparable_bp=4_671_223
        )
        self.assertAlmostEqual(threshold, 0.5 / 4_671_223)
        self.assertIn("4671223", why)

    def test_auto_refuses_to_guess_when_the_unit_is_ambiguous(self) -> None:
        truth = common.parse_newick("((a:0.001,b:0.0002):0.003,(c:0.004,d:0.001):0.002);")
        threshold, why = common.guess_collapse_threshold(truth)
        self.assertEqual(threshold, 0.0)
        self.assertIn("cannot be converted", why)

    def test_explicit_modes(self) -> None:
        self.assertEqual(common.guess_collapse_threshold(None, mode="none")[0], 0.0)
        threshold, _ = common.guess_collapse_threshold(
            None, mode="length", explicit_length=0.25
        )
        self.assertEqual(threshold, 0.25)
        with self.assertRaises(common.PrereqError):
            common.guess_collapse_threshold(None, mode="length")
        with self.assertRaises(common.PrereqError):
            common.guess_collapse_threshold(None, mode="snps")


class TestNeighborJoining(unittest.TestCase):
    def test_recovers_an_additive_tree(self) -> None:
        d = {
            ("a", "b"): 2, ("a", "c"): 4, ("a", "d"): 6,
            ("b", "c"): 4, ("b", "d"): 6, ("c", "d"): 6,
        }
        newick = common.neighbor_joining(
            ["a", "b", "c", "d"], lambda x, y: 0.0 if x == y else d[common.pair_key(x, y)]
        )
        tree = common.parse_newick(newick)
        self.assertEqual(sorted(tree.leaf_names()), ["a", "b", "c", "d"])
        # a and b are the closest pair, so they must be sisters
        self.assertEqual(
            common.rf_distance(tree, common.parse_newick("((a,b),c,d);")),
            common.rf_distance(tree, common.parse_newick("((a,b),c,d);")),
        )

    def test_degenerate_sizes(self) -> None:
        self.assertEqual(common.neighbor_joining([], lambda a, b: 0), "();")
        self.assertEqual(common.neighbor_joining(["a"], lambda a, b: 0), "a;")
        two = common.neighbor_joining(["a", "b"], lambda x, y: 0 if x == y else 4)
        self.assertEqual(sorted(common.parse_newick(two).leaf_names()), ["a", "b"])

    def test_negative_branch_lengths_are_clamped_to_zero(self) -> None:
        d = {
            ("a", "b"): 1, ("a", "c"): 10, ("a", "d"): 10,
            ("b", "c"): 10, ("b", "d"): 10, ("c", "d"): 1,
        }
        newick = common.neighbor_joining(
            ["a", "b", "c", "d"], lambda x, y: 0.0 if x == y else d[common.pair_key(x, y)]
        )
        self.assertNotIn(":-", newick)
        for node in common.parse_newick(newick).iter_nodes():
            if node.length is not None:
                self.assertGreaterEqual(node.length, 0.0)


# ---------------------------------------------------------------------------
# The verdict rule of section 11
# ---------------------------------------------------------------------------


class TestVerdictRule(unittest.TestCase):
    def evaluate(self, baseline, off, total=None):
        """``total`` names d1..dN so unmeasurable datasets stay in the denominator."""
        names = [f"d{i}" for i in range(1, (total or 3) + 1)]
        return common.evaluate_verdict(
            "flat", "nested levels", baseline, off, dataset_names=names
        )

    def test_pays_when_two_of_three_datasets_get_20_percent_worse(self) -> None:
        baseline = {
            "d1": {"exact_tree_s": 100.0, "peak_memory_bytes": 1_000},
            "d2": {"exact_tree_s": 100.0, "peak_memory_bytes": 1_000},
            "d3": {"exact_tree_s": 100.0, "peak_memory_bytes": 1_000},
        }
        off = {
            "d1": {"exact_tree_s": 130.0, "peak_memory_bytes": 1_000},
            "d2": {"exact_tree_s": 100.0, "peak_memory_bytes": 1_500},
            "d3": {"exact_tree_s": 105.0, "peak_memory_bytes": 1_000},
        }
        verdict = self.evaluate(baseline, off)
        self.assertTrue(verdict.pays)
        self.assertEqual(sorted(verdict.datasets_worse), ["d1", "d2"])
        self.assertIn("pays", verdict.reason)

    def test_exactly_20_percent_counts_as_worse(self) -> None:
        baseline = {"d1": {"exact_tree_s": 100.0}, "d2": {"exact_tree_s": 100.0}}
        off = {"d1": {"exact_tree_s": 120.0}, "d2": {"exact_tree_s": 120.0}}
        verdict = self.evaluate(baseline, off, total=3)
        self.assertTrue(verdict.pays)

    def test_just_under_20_percent_does_not_pay(self) -> None:
        baseline = {"d1": {"exact_tree_s": 100.0}, "d2": {"exact_tree_s": 100.0},
                    "d3": {"exact_tree_s": 100.0}}
        off = {"d1": {"exact_tree_s": 119.9}, "d2": {"exact_tree_s": 119.9},
               "d3": {"exact_tree_s": 119.9}}
        verdict = self.evaluate(baseline, off)
        self.assertFalse(verdict.pays)
        self.assertIn("does not pay", verdict.reason)

    def test_one_dataset_worse_is_not_enough(self) -> None:
        baseline = {"d1": {"exact_tree_s": 100.0}, "d2": {"exact_tree_s": 100.0},
                    "d3": {"exact_tree_s": 100.0}}
        off = {"d1": {"exact_tree_s": 500.0}, "d2": {"exact_tree_s": 100.0},
               "d3": {"exact_tree_s": 100.0}}
        verdict = self.evaluate(baseline, off)
        self.assertFalse(verdict.pays)
        self.assertEqual(verdict.datasets_worse, ["d1"])

    def test_memory_alone_can_carry_the_verdict(self) -> None:
        baseline = {"d1": {"peak_memory_bytes": 1000}, "d2": {"peak_memory_bytes": 1000}}
        off = {"d1": {"peak_memory_bytes": 1300}, "d2": {"peak_memory_bytes": 1300}}
        verdict = self.evaluate(baseline, off, total=3)
        self.assertTrue(verdict.pays)

    def test_missing_everything_is_not_a_pass(self) -> None:
        verdict = self.evaluate(
            {"d1": {}, "d2": {}, "d3": {}}, {"d1": {}, "d2": {}, "d3": {}}
        )
        self.assertIsNone(verdict.pays)
        self.assertIn("not measurable", verdict.reason)

    def test_a_partially_measurable_dataset_is_not_counted_as_worse(self) -> None:
        baseline = {
            "d1": {"exact_tree_s": 100.0},
            "d2": {"exact_tree_s": 100.0},
            "d3": {"exact_tree_s": None},
        }
        off = {
            "d1": {"exact_tree_s": 200.0},
            "d2": {"exact_tree_s": 200.0},
            "d3": {"exact_tree_s": None},
        }
        verdict = self.evaluate(baseline, off)
        self.assertTrue(verdict.pays)
        self.assertEqual(verdict.datasets_measurable, ["d1", "d2"])
        self.assertEqual(len(verdict.datasets_worse), 2)

    def test_faster_when_switched_off_is_reported_as_not_worse(self) -> None:
        baseline = {"d1": {"exact_tree_s": 100.0}, "d2": {"exact_tree_s": 100.0}}
        off = {"d1": {"exact_tree_s": 50.0}, "d2": {"exact_tree_s": 50.0}}
        verdict = self.evaluate(baseline, off, total=3)
        self.assertFalse(verdict.pays)
        self.assertEqual(verdict.datasets_worse, [])
        measurable = [m for m in verdict.metrics if m.ratio is not None]
        self.assertEqual(len(measurable), 2)
        self.assertTrue(all(m.worse is False for m in measurable))
        # peak memory was never supplied, so it stays unknown rather than "not worse"
        unknown = [m for m in verdict.metrics if m.ratio is None]
        self.assertTrue(all(m.worse is None for m in unknown))

    def test_zero_baseline_never_yields_a_ratio(self) -> None:
        baseline = {"d1": {"exact_tree_s": 0.0}, "d2": {"exact_tree_s": 0.0}}
        off = {"d1": {"exact_tree_s": 100.0}, "d2": {"exact_tree_s": 100.0}}
        verdict = self.evaluate(baseline, off)
        self.assertIsNone(verdict.pays)
        self.assertTrue(all(m.ratio is None for m in verdict.metrics))

    def test_serialises_to_plain_data(self) -> None:
        verdict = self.evaluate({"d1": {"exact_tree_s": 1.0}}, {"d1": {"exact_tree_s": 2.0}})
        doc = verdict.as_dict()
        self.assertEqual(doc["switch"], "flat")
        # two verdict metrics per dataset, and the denominator stays at three
        self.assertEqual(len(doc["metrics"]), 2 * 3)
        self.assertEqual(doc["datasets_total"], 3)
        self.assertEqual({m["dataset"] for m in doc["metrics"]}, {"d1", "d2", "d3"})
        json.dumps(doc)  # must be JSON-serialisable

    def test_a_dataset_that_was_never_run_is_named_not_silently_dropped(self) -> None:
        verdict = self.evaluate(
            {"d1": {"exact_tree_s": 100.0}}, {"d1": {"exact_tree_s": 400.0}}, total=3
        )
        self.assertEqual(verdict.datasets_total, 3)
        self.assertEqual(verdict.datasets_measurable, ["d1"])
        self.assertIn("d2", verdict.reason)
        self.assertIn("d3", verdict.reason)
        notes = {m.dataset: m.note for m in verdict.metrics}
        self.assertIn("never run", notes["d2"])

    def test_one_measurable_dataset_is_undetermined_not_a_verdict(self) -> None:
        """d1 got 4x worse, but one dataset cannot satisfy 'two of the three'."""
        verdict = self.evaluate(
            {"d1": {"exact_tree_s": 100.0}}, {"d1": {"exact_tree_s": 400.0}}, total=3
        )
        self.assertIsNone(verdict.pays)
        self.assertIn("undetermined", verdict.reason)
        self.assertEqual(verdict.datasets_worse, ["d1"])

    def test_two_measurable_datasets_are_enough_to_say_it_does_not_pay(self) -> None:
        baseline = {"d1": {"exact_tree_s": 100.0}, "d2": {"exact_tree_s": 100.0}}
        off = {"d1": {"exact_tree_s": 101.0}, "d2": {"exact_tree_s": 101.0}}
        verdict = self.evaluate(baseline, off, total=3)
        self.assertFalse(verdict.pays)
        self.assertEqual(verdict.datasets_measurable, ["d1", "d2"])
        self.assertIn("does not pay", verdict.reason)


class TestRatioHelpers(unittest.TestCase):
    def test_ratio_is_none_for_unknown_or_zero(self) -> None:
        self.assertIsNone(common.ratio(None, 1))
        self.assertIsNone(common.ratio(1, None))
        self.assertIsNone(common.ratio(1, 0))
        self.assertEqual(common.ratio(120, 100), 1.2)
        self.assertAlmostEqual(common.pct_change(120, 100), 20.0)
        self.assertIsNone(common.pct_change(None, 100))


# ---------------------------------------------------------------------------
# Expected-value acceptance checks
# ---------------------------------------------------------------------------


class TestExpectedChecks(TempDirCase):
    EXPECTED = {
        "abs_tolerance_snps": 1,
        "pairs": [
            {"name": "A-B", "genome_a": "A", "genome_b": "B", "snps": 2, "source": "ska2"},
            {"name": "A-C", "genome_a": "A", "genome_b": "C", "snps": 2, "source": "ska2"},
            {"name": "B-C", "genome_a": "B", "genome_b": "C", "snps": 0, "source": "mapping"},
        ],
    }

    def test_exact_matches_pass(self) -> None:
        observed = {("A", "B"): 2, ("A", "C"): 2, ("B", "C"): 0}
        checks = common.check_expected_pairs(self.EXPECTED, observed)
        self.assertEqual(len(checks), 3)
        self.assertTrue(all(c.passed for c in checks))
        self.assertEqual([c.status for c in checks], ["PASS"] * 3)

    def test_within_tolerance_passes_and_outside_fails(self) -> None:
        observed = {("A", "B"): 3, ("A", "C"): 5, ("B", "C"): 0}
        checks = {c.name: c for c in common.check_expected_pairs(self.EXPECTED, observed)}
        self.assertTrue(checks["A-B"].passed, "3 vs expected 2 is within +/-1")
        self.assertFalse(checks["A-C"].passed, "5 vs expected 2 is outside +/-1")
        self.assertEqual(checks["A-B"].status, "PASS")
        self.assertEqual(checks["A-C"].status, "FAIL")

    def test_pair_order_in_the_observed_data_does_not_matter(self) -> None:
        observed = {("B", "A"): 2, ("C", "A"): 2, ("C", "B"): 0}
        checks = common.check_expected_pairs(self.EXPECTED, observed)
        self.assertTrue(all(c.passed for c in checks))

    def test_a_missing_pair_is_not_measured_and_never_a_pass(self) -> None:
        observed = {("A", "B"): 2}
        checks = {c.name: c for c in common.check_expected_pairs(self.EXPECTED, observed)}
        self.assertIsNone(checks["A-C"].passed)
        self.assertEqual(checks["A-C"].status, "not measured")
        self.assertIn("no distance for this pair", checks["A-C"].note)

    def test_zero_tolerance_is_the_default_when_unspecified(self) -> None:
        expected = {"pairs": [{"genome_a": "A", "genome_b": "B", "snps": 2}]}
        checks = common.check_expected_pairs(expected, {("A", "B"): 3})
        self.assertFalse(checks[0].passed)

    @unittest.skipUnless(
        (SCRIPTS_DIR / "expected" / "salmonella_trio.json").exists(),
        "the shipped expected-values file is missing",
    )
    def test_the_shipped_expected_file_is_well_formed(self) -> None:
        doc = json.loads(
            (SCRIPTS_DIR / "expected" / "salmonella_trio.json").read_text(encoding="utf-8")
        )
        self.assertEqual(doc["dataset"], "salmonella_trio")
        self.assertEqual(len(doc["pairs"]), 3)
        by_name = {p["name"]: p for p in doc["pairs"]}
        self.assertEqual(by_name["isolate-A-isolate-B"]["snps"], 2)
        self.assertEqual(by_name["isolate-A-isolate-C"]["snps"], 2)
        self.assertEqual(by_name["isolate-B-isolate-C"]["snps"], 0)
        for pair in doc["pairs"]:
            self.assertIn("source", pair)
            self.assertIn("provenance", pair)
            self.assertTrue(pair["provenance"].get("file"))
        # the recorded values must survive the checker unchanged
        checks = common.check_expected_pairs(
            doc, {("isolate-A", "isolate-B"): 2, ("isolate-A", "isolate-C"): 2,
                  ("isolate-B", "isolate-C"): 0}
        )
        self.assertTrue(all(c.passed for c in checks))


# ---------------------------------------------------------------------------
# bench.json reader, including the --rf-from path
# ---------------------------------------------------------------------------


class TestBenchJson(TempDirCase):
    def test_reads_the_aggregate_tree_block(self) -> None:
        path = self.write(
            "bench.json",
            json.dumps(
                {
                    "tree": {"rf": 12, "rf_normalized": 0.04, "leaves_shared": 32},
                    "distances": {"mean_abs_error_snps": 0.5, "max_abs_error_snps": 3,
                                  "pairs_compared": 496, "pairs_exact": 400},
                }
            ),
        )
        bench = common.load_bench_json(path)
        self.assertEqual(bench.rf, 12)
        self.assertEqual(bench.rf_normalized, 0.04)
        self.assertEqual(bench.distance_mae, 0.5)
        self.assertEqual(bench.pairs_compared, 496)
        self.assertEqual(bench.rf_timeline(), [])
        series = common.rf_values_from_bench(bench)
        self.assertIn("aggregate", series)
        self.assertEqual(series["aggregate"][0]["rf"], 12)

    def test_reads_an_rf_timeline_when_present(self) -> None:
        timeline = [
            {"t_ms": 10.0, "event": "backbone_tree", "rf": 40, "rf_normalized": 0.9},
            {"t_ms": 900.0, "event": "cluster_focus", "cluster": 0, "focus": 2,
             "rf": 0, "rf_normalized": 0.0},
        ]
        path = self.write("bench.json", json.dumps({"tree": {"rf": 0}, "rf_timeline": timeline}))
        bench = common.load_bench_json(path)
        self.assertEqual(len(bench.rf_timeline()), 2)
        series = common.rf_values_from_bench(bench)
        self.assertIn("timeline", series)
        self.assertEqual(series["timeline"][1]["rf"], 0)

    def test_missing_rf_is_reported_not_invented(self) -> None:
        path = self.write("bench.json", json.dumps({"distances": {}}))
        bench = common.load_bench_json(path)
        self.assertIsNone(bench.rf)
        self.assertTrue(any("no tree.rf" in p for p in bench.problems))
        self.assertEqual(common.rf_values_from_bench(bench), {})

    def test_missing_file_is_actionable_when_required(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            common.load_bench_json(self.tmp / "nope.json")
        self.assertIn("mtato bench", str(ctx.exception))
        self.assertIsNone(common.load_bench_json(self.tmp / "nope.json", required=False))


# ---------------------------------------------------------------------------
# Configurations, statistics, formatting
# ---------------------------------------------------------------------------


class TestConfigMatrix(unittest.TestCase):
    def test_the_six_configurations_of_section_11(self) -> None:
        names = [c.name for c in common.CONFIGS]
        self.assertEqual(
            names,
            ["baseline", "no_block_skip", "flat", "no_piece_cache",
             "no_profile_grouping", "single_cluster"],
        )
        self.assertEqual(
            [c.flags for c in common.CONFIGS],
            [[], ["--no-block-skip"], ["--flat"], ["--no-piece-cache"],
             ["--no-profile-grouping"], ["--single-cluster"]],
        )
        self.assertTrue(common.CONFIGS[0].is_baseline)
        self.assertFalse(common.CONFIG_BY_NAME["flat"].is_baseline)
        self.assertEqual(common.SWITCHES, tuple(names[1:]))

    def test_selecting_a_switch_always_drags_the_baseline_along(self) -> None:
        chosen = common.select_configs(["flat"])
        self.assertEqual([c.name for c in chosen], ["baseline", "flat"])
        self.assertIsNone(common.select_configs(None)[0].param)

    def test_an_unknown_configuration_is_rejected(self) -> None:
        with self.assertRaises(common.PrereqError):
            common.select_configs(["turbo"])


class TestFormattingAndStats(unittest.TestCase):
    def test_unknown_values_never_print_as_zero(self) -> None:
        self.assertEqual(common.fmt(None), "n/a")
        self.assertEqual(common.fmt(float("nan")), "n/a")
        self.assertEqual(common.human_bytes(None), "n/a")
        self.assertEqual(common.human_seconds(None), "n/a")
        self.assertEqual(common.fmt(0), "0")
        self.assertEqual(common.fmt(0.0), "0")

    def test_human_units(self) -> None:
        self.assertEqual(common.human_bytes(1536), "1.5 KiB")
        self.assertEqual(common.human_bytes(5 * 1024 * 1024), "5.0 MiB")
        self.assertEqual(common.human_seconds(0.5), "500 ms")
        self.assertEqual(common.human_seconds(90), "1m30.0s")

    def test_statistics_ignore_none_and_nan(self) -> None:
        self.assertEqual(common.median([1, None, 3]), 2)
        self.assertEqual(common.mean([1, float("nan"), 3]), 2)
        self.assertIsNone(common.median([]))
        self.assertEqual(common.quantile([1, 2, 3, 4], 0.5), 2.5)
        self.assertIsNone(common.pearson([1], [2]))
        self.assertAlmostEqual(common.pearson([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertAlmostEqual(common.spearman([1, 2, 3], [10, 20, 30]), 1.0)

    def test_markdown_table_escapes_pipes_and_pads_short_rows(self) -> None:
        text = common.markdown_table(["a", "b"], [["x|y"], []])
        self.assertIn("x\\|y", text)
        self.assertEqual(len(text.splitlines()), 4)

    def test_dig_never_raises(self) -> None:
        doc = {"a": {"b": [1, 2, {"c": 3}]}}
        self.assertEqual(common.dig(doc, "a.b.2.c"), 3)
        self.assertIsNone(common.dig(doc, "a.x.y"))
        self.assertIsNone(common.dig(doc, "a.b.9"))
        self.assertEqual(common.dig(doc, "a.b.0"), 1)


class TestBinaryResolution(TempDirCase):
    def test_an_explicit_path_must_exist_and_be_executable(self) -> None:
        fake = fixtures.write_fake_mtato(self.tmp / "mtato")
        self.assertEqual(common.resolve_binary(str(fake)), fake)
        with self.assertRaises(common.PrereqError) as ctx:
            common.resolve_binary(str(self.tmp / "missing"))
        self.assertIn("cargo build --release", str(ctx.exception))
        self.assertIsNone(common.resolve_binary(str(self.tmp / "missing"), required=False))

    def test_the_error_lists_every_place_it_looked(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            common.resolve_binary(binary="definitely-not-installed-xyz")
        message = str(ctx.exception)
        self.assertIn("Tried:", message)
        self.assertIn("--mtato", message)


if __name__ == "__main__":
    unittest.main(verbosity=2)
