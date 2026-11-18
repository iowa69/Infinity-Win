#!/usr/bin/env python3
"""End-to-end tests for the benchmark scripts, driven by a fake ``mtato``.

The real binary does not exist yet, so ``scripts/tests/fixtures.py`` writes an
executable double that speaks the documented CLI contract (PROJECT.md section 9)
and emits every output file. That makes it possible to test the parts of this
deliverable that matter most:

* a configuration that fails is **reported as failed, never omitted**
  (PROJECT.md section 2 anti-goal);
* the matrix is **resumable**;
* ``--dry-run`` changes nothing;
* a missing prerequisite produces **one actionable line, not a traceback**;
* the accuracy curves come out of ``events.jsonl``, and the same numbers can be
  read from a ``bench.json`` instead (``--rf-from``);
* the verdict rule and the acceptance check reach ``REPORT.md``.

Run with::

    python3 -m unittest discover -s scripts/tests -t scripts -v

All fixtures live in ``tempfile.mkdtemp()`` directories under ``/tmp``.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional
from unittest import mock

TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent
for path in (str(SCRIPTS_DIR), str(TESTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

import accuracy_over_time  # noqa: E402
import baselines  # noqa: E402
import common  # noqa: E402
import fixtures  # noqa: E402
import make_report  # noqa: E402
import real_datasets  # noqa: E402
import run_matrix  # noqa: E402

CONFIG_NAMES = [c.name for c in common.CONFIGS]


class ScriptCase(unittest.TestCase):
    """Base class: a temp root, a fake mtato and a fixture bench-work tree."""

    genomes = 8
    dataset_names = ("fix_small",)

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="mtato-itest-")
        self.root = Path(self._tmp.name)
        self.work = self.root / "work"
        self.matrix = self.work / "matrix"
        self.mtato = fixtures.write_fake_mtato(self.root / "mtato")
        fixtures.make_work_dir(
            self.work, dataset_names=self.dataset_names, genomes_per_dataset=self.genomes
        )
        self._env_backup = dict(os.environ)
        os.environ["FAKE_MTATO_SLEEP"] = "0.01"

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self._env_backup)
        self._tmp.cleanup()

    # -- helpers ---------------------------------------------------------
    def truth_tree(self, dataset: str = "fix_small") -> Path:
        return self.work / dataset / "truth" / "truth.nwk"

    def run_dir(self, config: str, dataset: str = "fix_small") -> Path:
        return self.matrix / dataset / config

    def summary(self) -> Dict[str, Any]:
        return json.loads((self.matrix / "summary.json").read_text(encoding="utf-8"))

    def cells(self, summary: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
        return (summary or self.summary()).get("cells", [])

    def cell(self, config: str, summary: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        for c in self.cells(summary):
            if c.get("config") == config:
                return c
        raise AssertionError(f"no cell for {config!r} in the summary - it was omitted")

    def matrix_argv(self, *extra: str) -> List[str]:
        return [
            "--work-dir", str(self.work),
            "--matrix-dir", str(self.matrix),
            "--mtato", str(self.mtato),
            *extra,
        ]


# ---------------------------------------------------------------------------
# run_matrix.py
# ---------------------------------------------------------------------------


class TestRunMatrix(ScriptCase):
    def test_every_configuration_runs_and_is_recorded(self) -> None:
        code = run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        self.assertEqual(code, 0)
        summary = self.summary()
        self.assertEqual(len(summary["cells"]), 6)
        self.assertEqual([c["config"] for c in summary["cells"]], CONFIG_NAMES)
        for c in summary["cells"]:
            self.assertEqual(c["status"], "ok", f"{c['config']}: {c.get('error')}")
            self.assertIsNotNone(c["wall_s"], "wall time must be recorded per configuration")
            self.assertGreater(c["wall_s"], 0)
            self.assertTrue((Path(c["out_dir"]) / "run.json").exists())
            self.assertEqual(c["run"]["genomes"], self.genomes)
            self.assertEqual(c["bench"]["status"], "ok")
            self.assertEqual(c["run"]["switch_verified"], True)

    def test_the_baseline_command_line_has_no_switch_and_each_switch_has_its_own(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        summary = self.summary()
        by_config = {c["config"]: c for c in summary["cells"]}
        self.assertNotIn("--flat", by_config["baseline"]["run_command"])
        self.assertEqual(by_config["flat"]["run_command"][-1], "--flat")
        self.assertEqual(by_config["no_block_skip"]["run_command"][-1], "--no-block-skip")
        self.assertEqual(by_config["single_cluster"]["run_command"][-1], "--single-cluster")
        for c in summary["cells"]:
            self.assertIn("-o", c["run_command"])
            self.assertEqual(c["run_command"][1], "run")

    def test_a_failing_configuration_is_reported_not_omitted(self) -> None:
        os.environ["FAKE_MTATO_FAIL"] = "1"
        code = run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        self.assertEqual(code, 1, "a failed matrix must exit non-zero")
        summary = self.summary()
        self.assertEqual(len(summary["cells"]), 6, "failed cells must stay in the summary")
        for c in summary["cells"]:
            self.assertEqual(c["status"], "failed")
            self.assertIn("17", str(c.get("error")), "the exit code must be reported")
        self.assertEqual(len(summary["failures"]), 6)
        self.assertEqual(summary["totals"]["failed"], 6)

    def test_exit_zero_without_run_json_counts_as_a_failure(self) -> None:
        os.environ["FAKE_MTATO_NO_RUN_JSON"] = "1"
        code = run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        self.assertEqual(code, 1)
        c = self.cell("baseline")
        self.assertEqual(c["status"], "failed")
        self.assertIn("no run.json", str(c.get("error")))

    def test_a_failing_bench_does_not_hide_a_successful_run(self) -> None:
        """The run still produced timings; only the accuracy measurement is lost.

        So the cell keeps its run status, the bench failure is recorded on the
        cell *and* in the summary failure list, and the exit code is non-zero.
        """
        os.environ["FAKE_MTATO_BENCH_FAIL"] = "1"
        code = run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "flat"))
        self.assertEqual(code, 1)
        c = self.cell("flat")
        self.assertEqual(c["status"], "ok", "the run itself succeeded")
        self.assertTrue(c["run"]["run_json_present"], "the run itself did succeed")
        self.assertIsNotNone(c["exact_tree_s"], "timings survive a bench failure")
        self.assertEqual(c["bench"]["status"], "failed")
        self.assertIsNone(c["bench"].get("rf"))
        self.assertIn("simulated failure", c["bench"]["reason"])
        self.assertIn("fix_small/flat: bench failed", self.summary()["failures"])

    def test_a_switch_the_cli_ignored_is_flagged(self) -> None:
        os.environ["FAKE_MTATO_IGNORE_SWITCH"] = "1"
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "flat"))
        c = self.cell("flat")
        self.assertIs(c["run"]["switch_verified"], False)
        self.assertTrue(any("SWITCH NOT APPLIED" in p for p in c["problems"]))

    def test_missing_timings_fall_back_to_labelled_proxies(self) -> None:
        """summary.json records what the tool said; the report derives the proxies.

        Keeping the two apart matters: a proxy must never be written into the raw
        record as though the tool had reported it (PROJECT.md section 0 rule 8).
        """
        os.environ["FAKE_MTATO_NO_TIMINGS"] = "1"
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        c = self.cell("baseline")
        self.assertEqual(c["status"], "ok")
        self.assertIsNone(c["run"]["timings"]["exact_tree_ms"])
        self.assertIsNone(c["first_tree_s"], "the raw record must not invent a number")
        self.assertIsNone(c["exact_tree_s"])
        # the fallbacks live in make_report, and they are labelled as such
        m = make_report.measurements_for_cell(c)
        self.assertIsNotNone(m["first_tree_s"])
        self.assertEqual(m["first_tree_source"], "first tree event in events.jsonl")
        self.assertIsNotNone(m["exact_tree_s"])
        self.assertIn("proxy", m["exact_tree_source"])
        self.assertTrue(any("wall time" in p for p in m["problems"]))

    def test_missing_peak_memory_is_recorded_as_a_floor_not_as_zero(self) -> None:
        os.environ["FAKE_MTATO_NO_PEAKMEM"] = "1"
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        c = self.cell("baseline")
        self.assertIsNone(c["run"]["peak_memory_bytes"])
        self.assertIsNotNone(c["memory"]["bytes"], "the /proc fallback should have sampled something")
        self.assertTrue(c["memory"]["is_floor"])
        self.assertIn("VmHWM", c["memory"]["source"])

    def test_a_second_run_is_resumable_and_does_not_re_execute(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        first = self.cell("baseline")
        marker = self.run_dir("baseline") / "mtato.run.log"
        self.assertTrue(marker.exists())
        before = marker.stat().st_mtime_ns
        code = run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        self.assertEqual(code, 0)
        second = self.cell("baseline")
        self.assertEqual(second["status"], "cached")
        self.assertTrue(second["reused"])
        self.assertEqual(marker.stat().st_mtime_ns, before, "the run must not have been repeated")
        self.assertEqual(second["run"]["genomes"], first["run"]["genomes"])

    def test_force_re_runs_a_cached_configuration(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline", "--force"))
        c = self.cell("baseline")
        self.assertEqual(c["status"], "ok")
        self.assertFalse(c["reused"])

    def test_dry_run_prints_commands_and_writes_nothing(self) -> None:
        code = run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--dry-run"))
        self.assertEqual(code, 0)
        self.assertFalse((self.matrix / "summary.json").exists())
        self.assertFalse(self.run_dir("baseline").exists())

    def test_a_missing_binary_is_an_actionable_error(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            run_matrix.main(
                ["--work-dir", str(self.work), "--mtato", str(self.root / "no-such-mtato")]
            )
        self.assertIn("cargo build --release", str(ctx.exception))

    def test_a_dataset_with_no_truth_skips_bench_with_a_reason(self) -> None:
        fixtures.make_work_dir(
            self.work, dataset_names=("no_truth",), genomes_per_dataset=6, with_truth=False
        )
        code = run_matrix.main(self.matrix_argv("--datasets", "no_truth", "--configs", "baseline"))
        self.assertEqual(code, 0, "no truth is a skipped measurement, not a failure")
        summary = self.summary()
        cell = [c for c in summary["cells"] if c["dataset"] == "no_truth"][0]
        self.assertEqual(cell["status"], "ok")
        self.assertEqual(cell["bench"]["status"], "skipped")
        self.assertIn("no truth reference", cell["bench"]["reason"])

    def test_selecting_a_switch_also_runs_the_baseline(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "flat"))
        self.assertEqual(
            sorted(c["config"] for c in self.cells()), ["baseline", "flat"],
            "the verdict rule compares against the baseline, so it must always be run",
        )

    def test_the_summary_is_written_even_when_later_cells_fail(self) -> None:
        # First configuration succeeds, then the fake starts failing.
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))
        os.environ["FAKE_MTATO_FAIL"] = "1"
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "flat"))
        summary = self.summary()
        statuses = {c["config"]: c["status"] for c in summary["cells"]}
        self.assertEqual(statuses["baseline"], "cached")
        self.assertEqual(statuses["flat"], "failed")


# ---------------------------------------------------------------------------
# accuracy_over_time.py
# ---------------------------------------------------------------------------


class TestAccuracyOverTime(ScriptCase):
    def setUp(self) -> None:
        super().setUp()
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline", "flat"))

    def test_internal_rf_replays_the_event_stream_into_a_curve(self) -> None:
        out = self.root / "curves"
        code = accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.work),
                "--truth-tree", str(self.truth_tree()),
                "--out-dir", str(out),
                "--rf-mode", "internal",
                "--no-plot",
            ]
        )
        self.assertEqual(code, 0)
        csv_path = out / "accuracy_over_time__fix_small__baseline.csv"
        self.assertTrue(csv_path.exists())
        rows = self.read_csv(csv_path)
        self.assertGreaterEqual(len(rows), 4)
        # the fake run tree equals truth, so the final RF must be exactly zero
        self.assertEqual(float(rows[-1]["rf"]), 0.0)
        self.assertEqual(float(rows[-1]["rf_normalized"]), 0.0)
        # accuracy grows: the last point covers every taxon
        self.assertEqual(int(rows[-1]["leaves_shared"]), self.genomes)
        self.assertTrue(all(float(r["t_ms"]) >= 0 for r in rows))
        times = [float(r["t_ms"]) for r in rows]
        self.assertEqual(times, sorted(times), "the curve must be in time order")

    def test_a_wrong_tree_gives_a_nonzero_rf(self) -> None:
        """The internal RF path must actually detect a topology difference."""
        truth = common.read_newick_file(self.truth_tree())
        wrong_truth = self.root / "wrong_truth.nwk"
        wrong_truth.write_text(
            fixtures.perturb_tree(common.write_newick(truth), swaps=2) + "\n", encoding="utf-8"
        )
        self.assertNotEqual(
            wrong_truth.read_text(encoding="utf-8").strip(),
            self.truth_tree().read_text(encoding="utf-8").strip(),
            "the perturbed tree must differ from the real one",
        )
        out = self.root / "curves2"
        accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("flat")),
                "--work-dir", str(self.work),
                "--truth-tree", str(wrong_truth),
                "--out-dir", str(out),
                "--rf-mode", "internal",
                "--no-plot",
            ]
        )
        rows = self.read_csv(out / "accuracy_over_time__fix_small__flat.csv")
        self.assertGreater(float(rows[-1]["rf"]), 0.0)
        self.assertGreater(float(rows[-1]["rf_normalized"]), 0.0)

    def test_rf_can_come_from_a_bench_json_instead_of_being_recomputed(self) -> None:
        bench_json = self.run_dir("baseline") / "bench.json"
        self.assertTrue(bench_json.exists(), "run_matrix should have produced a bench.json")
        out = self.root / "frombench"
        code = accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.work),
                "--rf-from", str(bench_json),
                "--rf-mode", "bench",
                "--out-dir", str(out),
                "--no-plot",
            ]
        )
        self.assertEqual(code, 0)
        rows = self.read_csv(out / "accuracy_over_time__fix_small__baseline.csv")
        self.assertTrue(rows)
        self.assertTrue(all(r["rf_source"].startswith("bench.json") for r in rows))
        self.assertEqual(float(rows[-1]["rf"]), 0.0)

    def test_rf_mode_bench_without_a_bench_json_is_an_actionable_error(self) -> None:
        (self.run_dir("baseline") / "bench.json").unlink()
        out = self.root / "nobench"
        res = accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.work),
                "--truth-tree", str(self.truth_tree()),
                "--rf-mode", "bench",
                "--out-dir", str(out),
                "--no-plot",
            ]
        )
        # a single run with no points exits non-zero rather than pretending success
        self.assertEqual(res, 1)

    def test_matrix_mode_writes_one_curve_per_configuration(self) -> None:
        out = self.root / "matrix-curves"
        code = accuracy_over_time.main(
            [
                "--matrix-dir", str(self.matrix),
                "--work-dir", str(self.work),
                "--dataset", "fix_small",
                "--out-dir", str(out),
                "--no-plot",
            ]
        )
        self.assertEqual(code, 0)
        combined = out / "accuracy_over_time.csv"
        self.assertTrue(combined.exists())
        rows = self.read_csv(combined)
        self.assertEqual({r["config"] for r in rows}, {"baseline", "flat"})
        self.assertTrue((out / "accuracy_over_time__fix_small__baseline.csv").exists())
        self.assertTrue((out / "accuracy_over_time__fix_small__flat.csv").exists())

    def test_missing_events_file_is_reported_clearly(self) -> None:
        (self.run_dir("baseline") / "events.jsonl").unlink()
        out = self.root / "noevents"
        code = accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.work),
                "--truth-tree", str(self.truth_tree()),
                "--out-dir", str(out),
                "--rf-mode", "internal",
                "--no-plot",
            ]
        )
        self.assertEqual(code, 1)
        csv_path = out / "accuracy_over_time__fix_small__baseline.csv"
        self.assertFalse(csv_path.exists() and len(self.read_csv(csv_path)) > 0)

    def test_no_truth_and_no_bench_json_yields_no_points_and_says_why(self) -> None:
        (self.run_dir("baseline") / "bench.json").unlink()
        out = self.root / "notruth"
        code = accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.root / "empty-work"),
                "--truth-tree", str(self.root / "absent.nwk"),
                "--out-dir", str(out),
                "--rf-mode", "auto",
                "--no-plot",
            ]
        )
        self.assertEqual(code, 1, "with no RF source at all the script must not claim success")
        # the header-only CSV still exists, so the gap is visible as an empty file
        csv_path = out / "accuracy_over_time__fix_small__baseline.csv"
        self.assertTrue(csv_path.exists())
        self.assertEqual(self.read_csv(csv_path), [])

    def test_the_bench_timeline_is_used_when_no_truth_tree_is_available(self) -> None:
        out = self.root / "benchonly"
        code = accuracy_over_time.main(
            [
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.root / "empty-work"),
                "--truth-tree", str(self.root / "absent.nwk"),
                "--out-dir", str(out),
                "--rf-mode", "auto",
                "--no-plot",
            ]
        )
        self.assertEqual(code, 0)
        rows = self.read_csv(out / "accuracy_over_time__fix_small__baseline.csv")
        self.assertTrue(rows)
        self.assertTrue(all(r["rf_source"].startswith("bench.json") for r in rows))

    def test_the_figure_is_written_when_matplotlib_is_available(self) -> None:
        plt, reason = common.try_import_matplotlib()
        out = self.root / "fig"
        accuracy_over_time.main(
            [
                "--matrix-dir", str(self.matrix),
                "--work-dir", str(self.work),
                "--dataset", "fix_small",
                "--out-dir", str(out),
            ]
        )
        if plt is None:
            self.assertFalse((out / "accuracy_over_time.png").exists())
        else:
            self.assertTrue((out / "accuracy_over_time.png").exists(), reason)
        # the numbers survive either way
        self.assertTrue((out / "accuracy_over_time.csv").exists())

    def test_a_plotting_failure_still_leaves_the_csv(self) -> None:
        out = self.root / "figfail"
        with mock.patch.object(common, "try_import_matplotlib", return_value=(None, "simulated")):
            accuracy_over_time.main(
                [
                    "--matrix-dir", str(self.matrix),
                    "--work-dir", str(self.work),
                    "--dataset", "fix_small",
                    "--out-dir", str(out),
                ]
            )
        self.assertTrue((out / "accuracy_over_time.csv").exists())
        self.assertGreater(len(self.read_csv(out / "accuracy_over_time.csv")), 0)

    @staticmethod
    def read_csv(path: Path) -> List[Dict[str, str]]:
        import csv

        with open(path, newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# baselines.py
# ---------------------------------------------------------------------------


class TestBaselines(ScriptCase):
    def setUp(self) -> None:
        super().setUp()
        run_matrix.main(self.matrix_argv("--datasets", "fix_small", "--configs", "baseline"))

    def test_ska2_consumes_a_precomputed_table_when_the_tool_is_absent(self) -> None:
        with mock.patch.object(baselines.shutil, "which", return_value=None):
            code = baselines.main(
                ["ska2", "--dataset", "fix_small", "--work-dir", str(self.work), "--tree"]
            )
        self.assertEqual(code, 0)
        out = self.work / "baselines" / "fix_small"
        doc = json.loads((out / "ska2.json").read_text(encoding="utf-8"))
        self.assertEqual(doc["status"], "ok")
        self.assertFalse(doc["ran_here"], "the tool is absent, so it must not claim to have run")
        self.assertEqual(doc["stats"]["rows"], len(fixtures.pair_distances(fixtures.genome_names(self.genomes))))
        self.assertTrue((out / "ska2.dist.tsv").exists())
        self.assertTrue((out / "ska2.nj.tree.nwk").exists())

    def test_mash_degrades_to_a_reported_skip_instead_of_failing(self) -> None:
        with mock.patch.object(baselines.shutil, "which", return_value=None):
            code = baselines.main(
                ["mash", "--dataset", "fix_small", "--work-dir", str(self.work)]
            )
        self.assertEqual(code, 0, "a missing outside tool is a reported gap, not an error")
        doc = json.loads(
            (self.work / "baselines" / "fix_small" / "mash.json").read_text(encoding="utf-8")
        )
        self.assertEqual(doc["status"], "unavailable")
        self.assertTrue(any("not installed" in n for n in doc["notes"]))

    def test_mash_refuses_to_guess_when_told_not_to_run_and_nothing_is_cached(self) -> None:
        with mock.patch.object(baselines.shutil, "which", return_value=None):
            with self.assertRaises(common.PrereqError) as ctx:
                baselines.main(
                    ["mash", "--dataset", "fix_small", "--work-dir", str(self.work), "--no-run"]
                )
        self.assertIn("--mash-dist", str(ctx.exception))

    def test_compare_joins_ska2_against_the_mtato_cluster_distances(self) -> None:
        code = baselines.main(
            [
                "compare", "--dataset", "fix_small",
                "--run-dir", str(self.run_dir("baseline")),
                "--work-dir", str(self.work),
            ]
        )
        self.assertEqual(code, 0)
        out = self.work / "baselines" / "fix_small"
        doc = json.loads((out / "compare__baseline.json").read_text(encoding="utf-8"))
        summary = doc["summary"]
        self.assertEqual(summary["pairs_in_both"], 12, "8 genomes in 2 clusters = 2 x 6 pairs")
        self.assertEqual(summary["pairs_only_in_mtato"], 0)
        self.assertEqual(summary["exact_matches"], 12, "the fixture truth and run agree exactly")
        self.assertEqual(summary["mean_abs_diff_snps"], 0.0)
        self.assertAlmostEqual(summary["pearson"], 1.0)
        tsv = out / "ska2_vs_mtato__baseline.tsv"
        self.assertTrue(tsv.exists())
        self.assertEqual(len(common.read_table(tsv).rows), 12)
        tree_cmp = doc["tree_vs_reference"]
        self.assertEqual(tree_cmp["status"], "ok")
        self.assertEqual(tree_cmp["rf"], 0)

    def test_compare_reports_a_run_folder_with_no_distances(self) -> None:
        empty = self.root / "empty-run"
        empty.mkdir()
        with self.assertRaises(common.PrereqError) as ctx:
            baselines.main(
                ["compare", "--dataset", "fix_small", "--run-dir", str(empty),
                 "--work-dir", str(self.work)]
            )
        self.assertIn("cluster_<id>.dist.tsv", str(ctx.exception))

    def test_compare_without_a_reference_table_is_actionable(self) -> None:
        fixtures.make_work_dir(
            self.work, dataset_names=("no_truth",), genomes_per_dataset=6, with_truth=False
        )
        with self.assertRaises(common.PrereqError) as ctx:
            baselines.main(
                ["compare", "--dataset", "no_truth",
                 "--run-dir", str(self.run_dir("baseline")), "--work-dir", str(self.work)]
            )
        self.assertIn("--ska-dist", str(ctx.exception))

    def test_mash_table_parsing_handles_the_headerless_layout(self) -> None:
        """`mash dist -t` writes no header, and sample IDs are not numeric."""
        path = self.root / "mash.dist.tab"
        path.write_text("a\tb\t0.001\t0.0\t999\nc\td\t0.02\t0.0\t980\n", encoding="utf-8")
        rows, problems = baselines.parse_mash_table(path)
        self.assertEqual(len(rows), 2, "the first data row must not be eaten as a header")
        self.assertEqual([r["a"] for r in rows], ["a", "c"])
        self.assertAlmostEqual(rows[0]["distance"], 0.001)
        self.assertEqual(rows[0]["matching"], 999)
        self.assertTrue(any("positionally" in p for p in problems))

    def test_mash_table_parsing_also_handles_a_headered_table(self) -> None:
        path = self.root / "mash.dist.tsv"
        path.write_text(
            "ID\tID2\tdistance\tp-value\tmatching-hashes\na\tb\t0.001\t0.0\t999\n",
            encoding="utf-8",
        )
        rows, problems = baselines.parse_mash_table(path)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]["distance"], 0.001)
        self.assertFalse(any("positionally" in p for p in problems))

    def test_an_unparsable_mash_table_is_actionable(self) -> None:
        path = self.root / "mash.broken"
        path.write_text("only_one_column\nstill_one\n", encoding="utf-8")
        with self.assertRaises(common.PrereqError) as ctx:
            baselines.parse_mash_table(path)
        self.assertIn("mash dist -t", str(ctx.exception))

    def test_nj_tree_from_a_distance_table(self) -> None:
        rows = [
            {"a": "x", "b": "y", "distance": 1.0},
            {"a": "x", "b": "z", "distance": 5.0},
            {"a": "y", "b": "z", "distance": 5.0},
        ]
        newick, notes = baselines.nj_tree_from_rows(rows)
        self.assertEqual(sorted(common.parse_newick(newick).leaf_names()), ["x", "y", "z"])
        self.assertEqual(notes, [])

    def test_status_runs_without_any_tool_installed(self) -> None:
        with mock.patch.object(baselines.shutil, "which", return_value=None):
            self.assertEqual(baselines.main(["status", "--work-dir", str(self.work)]), 0)


# ---------------------------------------------------------------------------
# make_report.py
# ---------------------------------------------------------------------------


class TestMakeReport(ScriptCase):
    dataset_names = ("fix_small", "fix_pairs")
    genomes = 8

    def report(self, *extra: str) -> str:
        out = self.root / "REPORT.md"
        make_report.main(
            ["--work-dir", str(self.work), "--matrix-dir", str(self.matrix),
             "--out", str(out), *extra]
        )
        return out.read_text(encoding="utf-8")

    def test_the_report_covers_every_configuration_and_dataset(self) -> None:
        run_matrix.main(self.matrix_argv())
        text = self.report()
        self.assertIn("# M-Tato benchmark report", text)
        for config in CONFIG_NAMES:
            self.assertIn(config, text)
        for dataset in self.dataset_names:
            self.assertIn(dataset, text)
        self.assertIn("Verdict", text)
        self.assertIn("does not pay", text)
        self.assertIn("time to exact tree", text)
        self.assertIn("peak memory", text)

    def test_failed_cells_appear_in_the_report(self) -> None:
        os.environ["FAKE_MTATO_FAIL"] = "1"
        run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        del os.environ["FAKE_MTATO_FAIL"]
        run_matrix.main(self.matrix_argv("--datasets", "fix_pairs"))
        text = self.report()
        self.assertIn("FAILED", text)
        self.assertIn("Failures, gaps and data-quality notes", text)
        self.assertIn("simulated crash", text)
        # the failing dataset is still tabulated, not dropped
        self.assertIn("fix_small", text)

    def test_the_report_is_written_even_when_everything_failed(self) -> None:
        os.environ["FAKE_MTATO_FAIL"] = "1"
        run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        text = self.report()
        self.assertIn("failed", text)
        self.assertIn("not measurable", text.lower())

    def test_the_verdict_states_the_rule_and_never_claims_a_win_it_cannot_show(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        text = self.report()
        self.assertIn("20% worse", text)
        self.assertIn("at least 2 of the three datasets", text.replace("**", ""))
        # the fake makes every configuration identical, so nothing can "pay"
        self.assertNotIn("**pays**", text)

    def test_the_figure_and_its_numbers_are_both_referenced(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        text = self.report()
        self.assertIn("tree accuracy versus elapsed time", text.lower())
        self.assertIn("accuracy_over_time.csv", text)
        plt, _ = common.try_import_matplotlib()
        if plt is not None:
            self.assertIn(".png", text)

    def test_no_figures_flag_still_produces_a_complete_report(self) -> None:
        run_matrix.main(self.matrix_argv("--datasets", "fix_small"))
        text = self.report("--no-figures")
        self.assertIn("Robinson-Foulds", text)
        self.assertIn("figure generation disabled", text)

    def test_a_missing_summary_is_an_actionable_error(self) -> None:
        with self.assertRaises(common.PrereqError) as ctx:
            make_report.main(
                ["--work-dir", str(self.work), "--matrix-dir", str(self.root / "nope"),
                 "--out", str(self.root / "REPORT.md")]
            )
        self.assertIn("run_matrix.py", str(ctx.exception))


class TestAcceptanceCheck(ScriptCase):
    """The published trio numbers must reach the report as a PASS/FAIL check."""

    dataset_names = ("salmonella_trio",)
    genomes = 3
    expected_file = SCRIPTS_DIR / "expected" / "salmonella_trio.json"

    def setUp(self) -> None:
        super().setUp()
        names = ["isolate-A", "isolate-B", "isolate-C"]
        # rebuild the fixture dataset with the real genome names and truth
        import shutil

        shutil.rmtree(self.work / "salmonella_trio", ignore_errors=True)
        distances = [("isolate-A", "isolate-B", 2.0), ("isolate-A", "isolate-C", 2.0),
                     ("isolate-B", "isolate-C", 0.0)]
        fixtures.write_dataset(
            self.work,
            "salmonella_trio",
            names,
            truth_tree=fixtures.balanced_tree(names),
            truth_distances=distances,
            kind="real",
        )
        doc = json.loads((self.work / "salmonella_trio" / "dataset.json").read_text(encoding="utf-8"))
        doc["truth"]["expected_values"] = str(self.expected_file)
        (self.work / "salmonella_trio" / "dataset.json").write_text(json.dumps(doc, indent=2))
        self.dist_json = self.root / "trio_distances.json"
        self.dist_json.write_text(json.dumps([list(d) for d in distances]), encoding="utf-8")
        os.environ["FAKE_MTATO_DIST_JSON"] = str(self.dist_json)

    def report(self) -> str:
        out = self.root / "REPORT.md"
        make_report.main(
            ["--work-dir", str(self.work), "--matrix-dir", str(self.matrix),
             "--out", str(out), "--no-figures"]
        )
        return out.read_text(encoding="utf-8")

    def test_matching_distances_pass(self) -> None:
        self.assertTrue(self.expected_file.exists(), "the shipped expected-values file is required")
        run_matrix.main(self.matrix_argv("--datasets", "salmonella_trio", "--configs", "baseline"))
        text = self.report()
        self.assertIn("Acceptance check against published values", text)
        self.assertIn("isolate-A-isolate-B", text)
        self.assertIn("**PASS**", text)
        self.assertNotIn("**FAIL**", text)

    def test_wrong_distances_fail_loudly(self) -> None:
        self.dist_json.write_text(
            json.dumps([["isolate-A", "isolate-B", 40.0],
                        ["isolate-A", "isolate-C", 41.0],
                        ["isolate-B", "isolate-C", 39.0]]),
            encoding="utf-8",
        )
        run_matrix.main(self.matrix_argv("--datasets", "salmonella_trio", "--configs", "baseline"))
        text = self.report()
        self.assertIn("**FAIL**", text)
        self.assertNotIn("**PASS**", text)

    def test_a_pair_that_was_never_compared_is_not_a_pass(self) -> None:
        # Only one pair reported: the other two must be "not measured".
        self.dist_json.write_text(
            json.dumps([["isolate-A", "isolate-B", 2.0]]), encoding="utf-8"
        )
        run_matrix.main(self.matrix_argv("--datasets", "salmonella_trio", "--configs", "baseline"))
        text = self.report()
        self.assertIn("not measured", text)
        self.assertIn("PARTIAL", text)


# ---------------------------------------------------------------------------
# real_datasets.py
# ---------------------------------------------------------------------------


class TestRealDatasets(ScriptCase):
    def test_it_copies_verifies_and_records_provenance(self) -> None:
        src = self.root / "source" / "asm"
        src.mkdir(parents=True)
        names = ["alpha", "beta", "gamma", "delta"]
        for i, name in enumerate(names):
            fixtures.write_fasta(src / f"{name}.fasta", name, length=1024, seed=i)
        dest = self.root / "work2"
        with mock.patch.object(real_datasets, "NCBI_PD_ASM", src), \
             mock.patch.object(real_datasets, "SKA2_DISTANCE", self.root / "absent.tsv"), \
             mock.patch.object(real_datasets, "IQTREE_TREE", self.root / "absent.nwk"), \
             mock.patch.object(real_datasets, "IQTREE_BIONJ", self.root / "absent.bionj"):
            code = real_datasets.main(
                ["--datasets", "ncbi_pd_cluster", "--work-dir", str(dest)]
            )
        self.assertEqual(code, 0)
        doc = json.loads((dest / "ncbi_pd_cluster" / "dataset.json").read_text(encoding="utf-8"))
        self.assertEqual(doc["genome_count"], 4)
        self.assertTrue(doc["copy_verified"])
        copied = sorted(p.name for p in (dest / "ncbi_pd_cluster" / "asm").iterdir())
        self.assertEqual(copied, sorted(f"{n}.fasta" for n in names))
        # sources are untouched
        self.assertEqual(
            sorted(p.name for p in src.iterdir()), sorted(f"{n}.fasta" for n in names)
        )
        for entry in doc["genomes"]:
            self.assertEqual(entry["sha256"], entry["source_sha256"])
            self.assertEqual(entry["bytes"], entry["source_bytes"])
            self.assertTrue(entry["copy_verified"])
            self.assertEqual(entry["source"], str(src / entry["file"]))
        # missing reference files are reported, not ignored
        self.assertTrue(any("SKA2 distance table not present" in p for p in doc["problems"]))

    def test_a_missing_source_folder_is_reported_per_dataset(self) -> None:
        """One unusable source must not abort the others, and must not pass silently."""
        absent = self.root / "not-here" / "asm"
        dest = self.root / "w-missing"
        with mock.patch.object(real_datasets, "TRIO_ASM", absent):
            code = real_datasets.main(
                ["--datasets", "salmonella_trio", "--work-dir", str(dest)]
            )
        self.assertEqual(code, 1)
        self.assertFalse((dest / "salmonella_trio" / "dataset.json").exists())

    def test_sampling_is_deterministic_and_capped(self) -> None:
        closed = self.root / "closed"
        for cc, count in (("CC1", 4), ("CC2", 3)):
            for i in range(count):
                folder = closed / f"{cc}_GCF_{i:09d}.1"
                folder.mkdir(parents=True)
                fixtures.write_fasta(folder / f"{cc}_GCF_{i:09d}.1.fna", f"{cc}_{i}", length=512, seed=i)
        genomes = self.root / "genomes"
        genomes.mkdir()
        for i in range(5):
            fixtures.write_fasta(genomes / f"STRAIN{i}__GCF_{i}.fna", f"STRAIN{i}", length=512, seed=i)
        dest = self.root / "work3"
        with mock.patch.object(real_datasets, "MORI_CLOSED", closed), \
             mock.patch.object(real_datasets, "MORI_GENOMES", genomes):
            code = real_datasets.main(
                ["--datasets", "s_aureus_cc", "--work-dir", str(dest),
                 "--max-per-cc", "2", "--max-extra-genomes", "3"]
            )
        self.assertEqual(code, 0)
        doc = json.loads((dest / "s_aureus_cc" / "dataset.json").read_text(encoding="utf-8"))
        self.assertEqual(doc["genome_count"], 2 + 2 + 3)
        clusters = doc["truth"]["expected_clusters"]["clusters"]
        self.assertEqual({k: len(v) for k, v in clusters.items()}, {"CC1": 2, "CC2": 2})
        self.assertEqual(len(doc["truth"]["expected_clusters"]["unassigned"]), 3)
        self.assertTrue(
            any("far apart" in c.upper() or "FAR APART" in c for c in doc["caveats"]),
            "the 'different clonal complexes' caveat must be recorded",
        )
        # running again is idempotent and re-verifies rather than re-copying blindly
        with mock.patch.object(real_datasets, "MORI_CLOSED", closed), \
             mock.patch.object(real_datasets, "MORI_GENOMES", genomes):
            self.assertEqual(
                real_datasets.main(
                    ["--datasets", "s_aureus_cc", "--work-dir", str(dest),
                     "--max-per-cc", "2", "--max-extra-genomes", "3"]
                ),
                0,
            )

    def test_list_reports_source_availability(self) -> None:
        self.assertEqual(real_datasets.main(["--list"]), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
