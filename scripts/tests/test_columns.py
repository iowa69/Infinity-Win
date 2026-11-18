#!/usr/bin/env python3
"""Forward compatibility with the output columns the engine is adding.

``cluster_<id>.dist.tsv`` gains ``raw_snps`` and ``raw_comparable_bp`` after the
existing columns (``snps`` / ``comparable_bp`` then become the filtered
core-genome values), ``cluster_<id>.variants.tsv`` gains a final ``filter``
column, and ``genomes.tsv`` may gain columns. Every reader must find columns by
header name - never by position - and the reports must show both distance
bases against the references when both exist.
"""

from __future__ import annotations

import json
import os
import random
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

NAMES = fixtures.genome_names(6, prefix="g")


def shuffle_columns(path: Path, seed: int = 3) -> None:
    """Rewrite a TSV with its columns in a different order (same content)."""
    table = common.read_table(path)
    header = list(table.header)
    random.Random(seed).shuffle(header)
    lines = ["\t".join(header)] + ["\t".join(row[h] for h in header) for row in table.rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class ReadersTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="infinity-cols-")
        self.root = Path(self._tmp.name)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def folder(self, new_columns: bool) -> Path:
        return fixtures.write_output_folder(
            self.root / ("new" if new_columns else "old"), names=NAMES, clusters=1,
            new_columns=new_columns,
        )

    def test_raw_columns_are_read_when_present_and_none_when_absent(self) -> None:
        new_rows, problems = common.load_cluster_distances(self.folder(True))
        self.assertEqual(problems, [])
        self.assertTrue(new_rows)
        for row in new_rows:
            self.assertEqual(row["raw_snps"], row["snps"] + 1)
            self.assertEqual(row["raw_comparable_bp"], 2048)
            self.assertEqual(row["comparable_bp"], 1900)
            self.assertTrue(row["a"].startswith("g"), "names, not the numeric id columns")
        old_rows, _ = common.load_cluster_distances(self.folder(False))
        self.assertTrue(all(r["raw_snps"] is None for r in old_rows))
        self.assertEqual(common.distance_basis(new_rows)["snps_label"], common.SNPS_LABEL_FILTERED)
        self.assertEqual(common.distance_basis(old_rows)["snps_label"], common.SNPS_LABEL_ONLY)

    def test_column_order_never_matters(self) -> None:
        folder = self.folder(True)
        before, _ = common.load_cluster_distances(folder)
        variants_before = make_report.variants_fingerprint(folder)
        recomputed_before = common.recompute_pairwise_snps(folder).pairs
        for path in list(folder.glob("cluster_*.tsv")) + [folder / "genomes.tsv"]:
            shuffle_columns(path)
        after, problems = common.load_cluster_distances(folder)
        self.assertEqual(problems, [])
        key = lambda r: (r["a"], r["b"])  # noqa: E731
        strip = lambda r: {k: v for k, v in r.items() if k != "_source"}  # noqa: E731
        self.assertEqual(
            [strip(r) for r in sorted(before, key=key)], [strip(r) for r in sorted(after, key=key)]
        )
        self.assertEqual(common.recompute_pairwise_snps(folder).pairs, recomputed_before)
        # the fingerprint is per column order, so it differs in text but not in size
        self.assertEqual(len(make_report.variants_fingerprint(folder)["rows"]),
                         len(variants_before["rows"]))

    def test_a_filter_column_does_not_disturb_the_variant_readers(self) -> None:
        folder = self.folder(True)
        table, _ = common.load_variants_tsv(folder, "0")
        self.assertEqual(table.header[-1], "filter")
        result = common.recompute_pairwise_snps(folder)
        self.assertTrue(result.pairs)
        self.assertEqual(result.problems, [])
        # one PASS row per non-pivot member plus one filtered row
        self.assertEqual(result.variants_total, len(NAMES))

    def test_extra_genomes_columns_are_ignored_by_name_lookups(self) -> None:
        folder = self.folder(True)
        members = common._cluster_members(folder)
        self.assertEqual(sorted(members["0"]), NAMES)
        table, problems = common.load_genomes_tsv(folder)
        self.assertEqual(problems, [])
        self.assertIn("core_bp", table.header)


class ReportsBothBasesTests(unittest.TestCase):
    """End to end: the fake engine writes the new layout, the report shows both."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="infinity-cols-e2e-")
        self.root = Path(self._tmp.name)
        self.work = self.root / "work"
        self.matrix = self.work / "matrix"
        self.infinity = fixtures.write_fake_infinity(self.root / "infinity")
        self._env = dict(os.environ)
        os.environ["FAKE_INFINITY_SLEEP"] = "0.01"
        os.environ["FAKE_INFINITY_NEW_COLUMNS"] = "1"

    def tearDown(self) -> None:
        os.environ.clear()
        os.environ.update(self._env)
        self._tmp.cleanup()

    def run_cells(self, dataset: str, *extra: str) -> None:
        run_matrix.main(
            ["--work-dir", str(self.work), "--matrix-dir", str(self.matrix),
             "--infinity", str(self.infinity), "--datasets", dataset, "--repeats", "1", *extra]
        )

    def report(self) -> str:
        out = self.root / "REPORT.md"
        make_report.main(["--work-dir", str(self.work), "--matrix-dir", str(self.matrix),
                          "--out", str(out), "--no-figures"])
        return out.read_text(encoding="utf-8")

    def test_switch_invariance_compares_raw_distances_too(self) -> None:
        fixtures.make_work_dir(self.work, dataset_names=("fix_small",), genomes_per_dataset=8)
        self.run_cells("fix_small", "--configs", "baseline", "flat")
        summary = json.loads((self.matrix / "summary.json").read_text(encoding="utf-8"))
        doc = make_report.check_switch_invariance(make_report.collect(summary))
        flat = doc["fix_small"]["configs"]["flat"]
        self.assertEqual(flat["status"], "identical")
        self.assertEqual(flat["raw_pairs_differing"], 0)
        # change only a raw count: property 6 must still notice
        dist = self.matrix / "fix_small" / "flat" / "cluster_0.dist.tsv"
        table = common.read_table(dist)
        table.rows[0]["raw_snps"] = str(int(table.rows[0]["raw_snps"]) + 5)
        common.write_tsv(dist, table.header, [[r[h] for h in table.header] for r in table.rows])
        doc = make_report.check_switch_invariance(make_report.collect(summary))
        flat = doc["fix_small"]["configs"]["flat"]
        self.assertEqual(flat["pairs_differing"], 0, "the filtered counts still agree")
        self.assertEqual(flat["raw_pairs_differing"], 1)
        self.assertEqual(flat["status"], "DEFECT")
        self.assertIn("raw pairs differing", self.report())

    def test_the_acceptance_check_shows_raw_and_core_counts(self) -> None:
        names = ["isolate-A", "isolate-B", "isolate-C"]
        distances = [(names[0], names[1], 2.0), (names[0], names[2], 2.0), (names[1], names[2], 0.0)]
        fixtures.write_dataset(self.work, "salmonella_trio", names,
                               truth_distances=distances, kind="real", simulator_truth=False)
        dist_json = self.root / "trio.json"
        dist_json.write_text(json.dumps([list(d) for d in distances]), encoding="utf-8")
        os.environ["FAKE_INFINITY_DIST_JSON"] = str(dist_json)
        self.run_cells("salmonella_trio", "--configs", "baseline")
        expected_dir = self.root / "expected"
        expected_dir.mkdir()
        (expected_dir / "trio.json").write_text(json.dumps({
            "dataset": "salmonella_trio", "abs_tolerance_snps": 0,
            "pairs": [{"genome_a": a, "genome_b": b, "snps": d} for a, b, d in distances],
        }), encoding="utf-8")
        summary = json.loads((self.matrix / "summary.json").read_text(encoding="utf-8"))
        entries, _ = make_report.run_expected_checks(
            self.work, summary, make_report.collect(summary), expected_dir=expected_dir
        )
        info = entries[0]["checks_by_config"]["baseline"]
        self.assertEqual(info["snps_label"], common.SNPS_LABEL_FILTERED)
        self.assertEqual(info["status"], "PASS", "the core counts equal the published values")
        self.assertEqual(info["status_raw"], "FAIL", "raw is one SNP over, tolerance 0")
        text = make_report.render_expected(entries)
        self.assertIn("engine core (`snps`)", text)
        self.assertIn("engine raw (`raw_snps`)", text)


if __name__ == "__main__":
    unittest.main()
