#!/usr/bin/env python3
"""Timing honesty: repeated runs, medians with ranges, and what is not judged.

Many matrix cells finish in well under a second, where a 20% difference is
noise. So every configuration runs several times, the section 11 rule is
applied to medians, a time whose baseline median is under 2 s is "too short to
judge", and the 3-genome trio (an exactness check) is not part of the
performance verdict at all.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List

TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent
for path in (str(SCRIPTS_DIR), str(TESTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

import common  # noqa: E402
import fixtures  # noqa: E402
import make_report  # noqa: E402
import run_matrix  # noqa: E402


def cell_with_repeats(dataset: str, config: str, exact: List[float], memory: List[int]) -> Dict[str, Any]:
    repeats = [
        {"repeat": i + 1, "status": "ok", "wall_s": e + 0.05, "first_tree_s": e / 3,
         "exact_tree_s": e, "total_s": e + 0.01, "peak_memory_bytes": m,
         "peak_memory_is_floor": False}
        for i, (e, m) in enumerate(zip(exact, memory))
    ]
    return {
        "dataset": dataset, "config": config, "status": "ok", "out_dir": None,
        "wall_s": repeats[0]["wall_s"], "exact_tree_s": exact[0],
        "peak_memory_bytes": memory[0], "repeats": repeats,
        "run": {"timings_s": {"exact_tree_s": exact[0]}, "genomes": 10},
    }


class VerdictOnMediansTests(unittest.TestCase):
    def measure(self, exact: List[float], memory: List[int]) -> Dict[str, Any]:
        return make_report.measurements_for_cell(cell_with_repeats("d", "c", exact, memory))

    def test_the_median_and_range_come_from_every_usable_repeat(self) -> None:
        m = self.measure([3.0, 9.0, 4.0], [100, 300, 200])
        self.assertEqual(m["exact_tree_s"], 4.0, "median, not the first or the mean")
        self.assertEqual(m["exact_tree_s_stats"]["min"], 3.0)
        self.assertEqual(m["exact_tree_s_stats"]["max"], 9.0)
        self.assertEqual(m["exact_tree_s_stats"]["n"], 3)
        self.assertEqual(m["peak_memory_bytes"], 200)
        self.assertEqual(m["repeats_usable"], 3)
        self.assertIn("[", make_report.fmt_stat(m, "exact_tree_s"))

    def test_a_failed_repeat_is_left_out_of_the_median(self) -> None:
        cell = cell_with_repeats("d", "c", [3.0, 9.0, 4.0], [100, 300, 200])
        cell["repeats"][1]["status"] = "failed"
        m = make_report.measurements_for_cell(cell)
        self.assertEqual(m["exact_tree_s_stats"]["n"], 2)
        self.assertEqual(m["exact_tree_s"], 3.5)

    def test_an_old_single_run_summary_is_n_equals_one(self) -> None:
        cell = {"dataset": "d", "config": "c", "status": "ok", "wall_s": 1.5,
                "peak_memory_bytes": 7, "run": {"timings_s": {"exact_tree_s": 1.2}}}
        m = make_report.measurements_for_cell(cell)
        self.assertEqual(m["exact_tree_s"], 1.2)
        self.assertEqual(m["exact_tree_s_stats"]["n"], 1)
        self.assertNotIn("[", make_report.fmt_stat(m, "exact_tree_s"))

    def verdict(self, base: Dict[str, Dict[str, Any]], off: Dict[str, Dict[str, Any]], **kw: Any) -> common.Verdict:
        return common.evaluate_verdict("flat", "x", base, off, dataset_names=sorted(base), **kw)

    def test_a_short_baseline_time_is_too_short_to_judge(self) -> None:
        base = {"a": self.measure([1.0, 1.1, 0.9], [100] * 3),
                "b": self.measure([1.0, 1.0, 1.0], [100] * 3)}
        off = {"a": self.measure([1.5, 1.6, 1.4], [100] * 3),
               "b": self.measure([1.5, 1.5, 1.5], [100] * 3)}
        judged = self.verdict(base, off)
        self.assertEqual(judged.state, "pays", "without the floor +50% on two datasets pays")
        floored = self.verdict(base, off, min_baseline_s=2.0)
        time_rows = [m for m in floored.metrics if m.metric == "time to exact tree"]
        self.assertTrue(all(m.state == "too_short" for m in time_rows))
        self.assertTrue(all(m.worse is None for m in time_rows))
        self.assertEqual(floored.datasets_worse, [])
        # memory is still judged, so both datasets remain measurable
        self.assertEqual(floored.datasets_measurable, ["a", "b"])
        self.assertEqual(floored.state, "does_not_pay")
        self.assertIn("too short to judge", floored.reason)

    def test_memory_still_carries_the_verdict_on_a_short_run(self) -> None:
        base = {d: self.measure([0.5] * 3, [100] * 3) for d in ("a", "b")}
        off = {d: self.measure([0.5] * 3, [200] * 3) for d in ("a", "b")}
        v = self.verdict(base, off, min_baseline_s=2.0)
        self.assertEqual(v.state, "pays")

    def test_a_long_enough_baseline_is_judged_on_medians(self) -> None:
        base = {d: self.measure([3.0, 30.0, 3.1], [100] * 3) for d in ("a", "b")}
        off = {d: self.measure([3.9, 3.8, 40.0], [100] * 3) for d in ("a", "b")}
        v = self.verdict(base, off, min_baseline_s=2.0)
        time_rows = [m for m in v.metrics if m.metric == "time to exact tree"]
        # medians 3.1 vs 3.9: +26%, although one outlier each way
        self.assertTrue(all(m.worse for m in time_rows))
        self.assertTrue(all(m.ranges_overlap for m in time_rows))
        self.assertEqual(v.state, "pays")

    def test_an_excluded_dataset_is_named_and_left_out_of_the_denominator(self) -> None:
        base = {d: self.measure([3.0] * 3, [100] * 3) for d in ("a", "b", "trio")}
        off = dict(base)
        off["trio"] = self.measure([30.0] * 3, [1000] * 3)
        v = self.verdict(base, off, excluded={"trio": "3 genomes: an exactness check"})
        self.assertEqual(v.datasets_total, 2)
        self.assertNotIn("trio", v.datasets_worse)
        self.assertNotIn("trio", v.datasets_measurable)
        self.assertEqual(v.datasets_excluded, {"trio": "3 genomes: an exactness check"})
        self.assertTrue(all(m.state == "excluded" for m in v.metrics if m.dataset == "trio"))
        self.assertIn("left out: trio", v.reason)

    def test_single_runs_never_claim_separated_ranges(self) -> None:
        base = {d: self.measure([3.0], [100]) for d in ("a", "b")}
        off = {d: self.measure([9.0], [100]) for d in ("a", "b")}
        v = self.verdict(base, off)
        self.assertTrue(all(m.ranges_overlap is None for m in v.metrics))

    def test_the_trio_is_excluded_by_default_and_a_dataset_can_opt_out(self) -> None:
        excluded = common.performance_exclusions(
            [{"name": "salmonella_trio"}, {"name": "x", "performance_verdict": False},
             {"name": "sim_100"}]
        )
        self.assertEqual(set(excluded), {"salmonella_trio", "x"})


class RepeatsInTheMatrixTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="infinity-rep-")
        self.root = Path(self._tmp.name)
        self.work = self.root / "work"
        self.matrix = self.work / "matrix"
        self.infinity = fixtures.write_fake_infinity(self.root / "infinity")
        fixtures.make_work_dir(self.work, dataset_names=("fix_small",), genomes_per_dataset=8)
        self._env = dict(os.environ)
        os.environ["FAKE_INFINITY_SLEEP"] = "0.01"

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self._env)
        self._tmp.cleanup()

    def run_matrix(self, *extra: str) -> int:
        return run_matrix.main(
            ["--work-dir", str(self.work), "--matrix-dir", str(self.matrix),
             "--infinity", str(self.infinity), "--datasets", "fix_small", *extra]
        )

    def cells(self) -> Dict[str, Dict[str, Any]]:
        summary = json.loads((self.matrix / "summary.json").read_text(encoding="utf-8"))
        return {c["config"]: c for c in summary["cells"]}

    def test_every_repeat_is_recorded_and_summarised(self) -> None:
        self.assertEqual(self.run_matrix("--configs", "baseline", "flat", "--repeats", "3"), 0)
        cells = self.cells()
        for config in ("baseline", "flat"):
            cell = cells[config]
            self.assertEqual([r["repeat"] for r in cell["repeats"]], [1, 2, 3])
            self.assertTrue(all(r["status"] == "ok" for r in cell["repeats"]))
            self.assertTrue(all(r["wall_s"] for r in cell["repeats"]))
            self.assertEqual(cell["timing"]["exact_tree_s"]["n"], 3)
            self.assertEqual(cell["timing"]["wall_s"]["n"], 3)
            self.assertEqual(cell["repeats_usable"], 3)
        summary = json.loads((self.matrix / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["request"]["repeats"], 3)

    def test_extra_repeats_keep_only_run_json_and_logs(self) -> None:
        self.run_matrix("--configs", "baseline", "--repeats", "2")
        r2 = self.matrix / "fix_small" / "baseline" / "repeats" / "r2"
        self.assertEqual(sorted(p.name for p in r2.iterdir()), sorted(run_matrix.REPEAT_KEEP))
        canonical = self.matrix / "fix_small" / "baseline"
        self.assertTrue((canonical / "events.jsonl").exists(), "repeat 1 keeps everything")
        self.run_matrix("--configs", "baseline", "--repeats", "2", "--force",
                        "--keep-repeat-outputs")
        self.assertTrue((r2 / "events.jsonl").exists())

    def test_repeats_are_resumable(self) -> None:
        self.run_matrix("--configs", "baseline", "--repeats", "2")
        log = self.matrix / "fix_small" / "baseline" / "repeats" / "r2" / "infinity.run.log"
        before = log.stat().st_mtime_ns
        self.run_matrix("--configs", "baseline", "--repeats", "3")
        cell = self.cells()["baseline"]
        self.assertEqual([r["status"] for r in cell["repeats"]], ["cached", "cached", "ok"])
        self.assertEqual(log.stat().st_mtime_ns, before, "a cached repeat must not re-run")

    def test_a_failed_first_run_is_not_repeated(self) -> None:
        os.environ["FAKE_INFINITY_FAIL"] = "1"
        self.run_matrix("--configs", "baseline", "--repeats", "3")
        cell = self.cells()["baseline"]
        self.assertEqual(cell["status"], "failed")
        self.assertEqual(len(cell["repeats"]), 1)
        self.assertTrue(any("not repeated" in p for p in cell["problems"]))

    def test_a_failed_later_repeat_is_reported_not_dropped(self) -> None:
        self.run_matrix("--configs", "baseline", "--repeats", "1")
        os.environ["FAKE_INFINITY_FAIL"] = "1"
        code = self.run_matrix("--configs", "baseline", "--repeats", "2")
        self.assertEqual(code, 1, "a failed repeat is a failure of the matrix")
        cell = self.cells()["baseline"]
        self.assertEqual(cell["status"], "cached")
        self.assertEqual(cell["repeats"][1]["status"], "failed")
        self.assertEqual(cell["timing"]["exact_tree_s"]["n"], 1, "the failure adds no number")
        summary = json.loads((self.matrix / "summary.json").read_text(encoding="utf-8"))
        self.assertTrue(any("repeat 2" in f for f in summary["failures"]))

    def test_dry_run_mentions_the_repeats(self) -> None:
        import contextlib
        import io

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            self.run_matrix("--configs", "baseline", "--repeats", "3", "--dry-run")
        self.assertIn("repeats 2..3", buffer.getvalue())
        self.assertFalse((self.matrix / "summary.json").exists())

    def test_zero_repeats_is_refused(self) -> None:
        with self.assertRaises(common.PrereqError):
            self.run_matrix("--repeats", "0")

    def test_the_report_shows_medians_ranges_and_the_rules(self) -> None:
        self.run_matrix("--repeats", "2")
        out = self.root / "REPORT.md"
        make_report.main(["--work-dir", str(self.work), "--matrix-dir", str(self.matrix),
                          "--out", str(out), "--no-figures"])
        text = out.read_text(encoding="utf-8")
        self.assertIn("Medians of repeated runs", text)
        self.assertIn("too short to judge", text, "the fake runs take well under 2 s")
        self.assertIn("baseline median [min-max]", text)
        self.assertIn("2/2", text, "usable/total runs per cell")
        self.assertIn("**simulated** genomes", text)
        self.assertNotIn("the real ones", text)


if __name__ == "__main__":
    unittest.main()
