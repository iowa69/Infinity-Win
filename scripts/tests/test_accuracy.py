#!/usr/bin/env python3
"""Unit tests for the full-leaf-set accuracy curve and the key figure.

The curve must score the tree the program was *showing* at every instant, over
every genome: clusters not yet refined count as unresolved stars, refined ones
as their emitted trees, joined by the backbone. These tests build tiny event
streams by hand so each rule can be checked in isolation.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Optional

TESTS_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = TESTS_DIR.parent
for path in (str(SCRIPTS_DIR), str(TESTS_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

import accuracy_over_time as aot  # noqa: E402
import common  # noqa: E402

CLUSTERS = {
    0: ["a1", "a2", "a3", "a4"],
    1: ["b1", "b2", "b3", "b4"],
    2: ["c1", "c2", "c3", "c4"],
}


def cluster_tree(names: List[str], inner: float = 2.0) -> str:
    """((x1,x2),(x3,x4)) with SNP-like integer branch lengths."""
    a, b, c, d = names
    return f"(({a}:1,{b}:1):{inner},({c}:1,{d}:1):{inner});"


def write_truth(root: Path) -> Path:
    truth = root / "truth"
    truth.mkdir(parents=True)
    for cid, names in CLUSTERS.items():
        (truth / f"tree_c{cid:02d}.nwk").write_text(cluster_tree(names) + "\n", encoding="utf-8")
    (truth / "clusters.tsv").write_text(
        "genome\tcluster_id\n"
        + "".join(f"{n}\t{cid}\n" for cid, names in CLUSTERS.items() for n in names),
        encoding="utf-8",
    )
    return truth


class EventStream:
    """Builds an events.jsonl record list in the real envelope shape."""

    def __init__(self) -> None:
        self.events: List[Dict[str, Any]] = []
        self.t = 0.0

    def emit(self, name: str, data: Dict[str, Any], dt: float = 10.0) -> None:
        self.t += dt
        self.events.append(
            {"seq": len(self.events), "t_ms": self.t, "wall_ms": 0,
             "event": {"event": name, "data": data}}
        )

    def standard_start(self, backbone: str = "((a1:0.02,b1:0.02):0.001,c1:0.02);") -> None:
        gid = 0
        for names in CLUSTERS.values():
            for name in names:
                self.emit("genome_scanned", {"id": gid, "name": name}, dt=1.0)
                gid += 1
        self.emit(
            "clustered",
            {"clusters": [
                {"id": cid, "pivot_name": names[0], "member_names": list(names)}
                for cid, names in CLUSTERS.items()
            ]},
        )
        self.emit("backbone_tree", {"newick": backbone, "pivots": [n[0] for n in CLUSTERS.values()]})

    def focus(self, cid: int, focus: int, newick: str) -> None:
        self.emit(
            "cluster_focus",
            {"cluster": cid, "focus": focus, "approximate": focus < 2, "newick": newick,
             "members": [{"name": n} for n in CLUSTERS[cid]]},
        )

    def log(self, root: Path) -> common.EventLog:
        path = root / "events.jsonl"
        path.write_text("".join(json.dumps(e) + "\n" for e in self.events), encoding="utf-8")
        return common.read_event_log(path)


class FullLeafSetTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory(prefix="infinity-acc-")
        self.root = Path(self._tmp.name)
        self.truth_dir = write_truth(self.root)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def truth(self) -> aot.FullTruth:
        truth, problems = aot.load_full_truth(
            self.truth_dir, None, collapse_snps=0.5, collapse_mode="auto",
            collapse_length=None, comparable_bp=None,
        )
        self.assertIsNotNone(truth, problems)
        return truth

    def score(self, stream: EventStream, comparable_bp: Optional[int] = 1000) -> List[Dict[str, Any]]:
        thresholds = aot.run_thresholds(
            collapse_snps=0.5, collapse_mode="auto", collapse_length=None,
            comparable_bp=comparable_bp,
        )
        rows, _problems = aot.score_full_leaf_set(
            stream.log(self.root), self.truth(), thresholds, "toy", "baseline"
        )
        return rows

    def test_simulator_truth_is_the_cluster_trees_joined_at_a_star(self) -> None:
        truth = self.truth()
        self.assertEqual(truth.kind, "simulator")
        self.assertEqual(len(truth.leaves), 12)
        self.assertEqual(len(truth.root.children), 3)
        index = common.leaf_index(truth.leaves)
        _, masks = common.clade_masks(truth.root, index, truth.threshold)
        # 3 cluster clades + 2 cherries in each of the 3 clusters
        self.assertEqual(len(common.canonical_splits(masks, (1 << 12) - 1)), 9)

    def test_unrefined_clusters_are_stars_and_the_curve_falls_to_zero(self) -> None:
        stream = EventStream()
        stream.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        for cid, names in CLUSTERS.items():
            stream.focus(cid, 2, cluster_tree(names))
        rows = self.score(stream)
        self.assertEqual(len(rows), 4)
        first = rows[0]
        self.assertEqual(first["event"], "backbone_tree")
        # stars: only the 3 cluster clades are claimed, 6 within-cluster splits missing
        self.assertEqual(first["splits_run"], 3)
        self.assertEqual(first["rf"], 6)
        self.assertAlmostEqual(first["rf_normalized"], 6 / 12)
        self.assertEqual(first["leaves_run"], 12)
        self.assertEqual(first["leaves_exact"], 0)
        self.assertEqual([r["rf"] for r in rows], [6, 4, 2, 0])
        self.assertEqual(rows[-1]["leaves_exact"], 12)
        self.assertEqual(rows[-1]["clusters_exact"], 3)

    def test_a_resolved_backbone_over_a_star_truth_counts_between_clusters(self) -> None:
        # four clusters would be needed for a resolved backbone over three to
        # differ; add a fourth cluster on the fly
        CLUSTERS[3] = ["d1", "d2", "d3", "d4"]
        try:
            (self.truth_dir / "tree_c03.nwk").write_text(cluster_tree(CLUSTERS[3]) + "\n", encoding="utf-8")
            stream = EventStream()
            stream.standard_start(backbone="((a1:0.02,b1:0.02):0.001,(c1:0.02,d1:0.02):0.001);")
            for cid, names in CLUSTERS.items():
                stream.focus(cid, 2, cluster_tree(names))
            rows = self.score(stream)
            last = rows[-1]
            # regression: genome "c1" (pivot of cluster 2) must not be read as the
            # generic alias "c1" of cluster 1
            self.assertEqual(last["note"], "", "every pivot must land on its own cluster")
            self.assertEqual(last["rf"], 1, "ab|cd is claimed by the backbone, not by the truth")
            self.assertEqual(last["rf_between_clusters"], 1)
            self.assertEqual(last["rf_within_clusters"], 0)
        finally:
            del CLUSTERS[3]

    def test_a_cluster_tree_rooted_on_the_wrong_edge_counts_within_clusters(self) -> None:
        stream = EventStream()
        stream.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        stream.focus(0, 2, "(a1:1,(a2:1,(a3:1,a4:1):2):1);")  # same unrooted tree, other root
        stream.focus(1, 2, cluster_tree(CLUSTERS[1]))
        stream.focus(2, 2, cluster_tree(CLUSTERS[2]))
        last = self.score(stream)[-1]
        self.assertGreater(last["rf"], 0)
        self.assertEqual(last["rf_between_clusters"], 0)
        self.assertEqual(last["rf_within_clusters"], last["rf"])

    def test_each_focus_level_is_collapsed_in_its_own_units(self) -> None:
        stream = EventStream()
        stream.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        # Focus 1 with 0.3 difference pieces inside: under half a piece, collapsed
        stream.focus(0, 1, cluster_tree(CLUSTERS[0], inner=0.3))
        # Focus 2 with 3 SNPs: kept
        stream.focus(1, 2, cluster_tree(CLUSTERS[1], inner=3.0))
        rows = self.score(stream)
        after_focus1, after_focus2 = rows[1], rows[2]
        self.assertEqual(after_focus1["splits_run"], 3, "a0.3-piece branch carries no evidence")
        self.assertEqual(after_focus2["splits_run"], 5)
        # A Focus 0 tree in subs/site: 0.0004 is 0.4 SNP over 1000 bp -> collapsed
        stream2 = EventStream()
        stream2.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        stream2.focus(0, 0, cluster_tree(CLUSTERS[0], inner=0.0004))
        self.assertEqual(self.score(stream2)[-1]["splits_run"], 3)
        stream3 = EventStream()
        stream3.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        stream3.focus(0, 0, cluster_tree(CLUSTERS[0], inner=0.004))
        self.assertEqual(self.score(stream3)[-1]["splits_run"], 5)

    def test_genomes_the_truth_does_not_know_are_left_out(self) -> None:
        stream = EventStream()
        stream.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        stream.emit("genome_scanned", {"id": 99, "name": "stranger"}, dt=1.0)
        rows = self.score(stream)
        self.assertEqual(rows[0]["leaves_shared"], 12)

    def test_fewer_than_four_scored_genomes_leaves_rf_undefined(self) -> None:
        tiny = self.root / "tiny"
        tiny.mkdir()
        (tiny / "tree.nwk").write_text("(x:1,y:1,z:1);\n", encoding="utf-8")
        truth, _ = aot.load_full_truth(
            None, tiny / "tree.nwk", collapse_snps=0.5, collapse_mode="none",
            collapse_length=None, comparable_bp=None,
        )
        stream = EventStream()
        for i, name in enumerate("xyz"):
            stream.emit("genome_scanned", {"id": i, "name": name}, dt=1.0)
        stream.emit("clustered", {"clusters": [{"id": 0, "pivot_name": "x", "member_names": list("xyz")}]})
        stream.emit("backbone_tree", {"newick": "x;"})
        thresholds = aot.run_thresholds(collapse_snps=0.5, collapse_mode="none",
                                        collapse_length=None, comparable_bp=None)
        rows, problems = aot.score_full_leaf_set(stream.log(self.root), truth, thresholds, "t", "c")
        self.assertTrue(rows)
        self.assertIsNone(rows[0]["rf"])
        self.assertIsNone(rows[0]["rf_normalized"])
        self.assertTrue(any("fewer than 4 taxa" in p for p in problems))

    def test_events_beyond_max_points_are_replayed_but_not_scored(self) -> None:
        stream = EventStream()
        stream.standard_start(backbone="(a1:0.02,b1:0.02,c1:0.02);")
        for cid, names in CLUSTERS.items():
            stream.focus(cid, 2, cluster_tree(names))
        thresholds = aot.run_thresholds(collapse_snps=0.5, collapse_mode="auto",
                                        collapse_length=None, comparable_bp=1000)
        rows, problems = aot.score_full_leaf_set(
            stream.log(self.root), self.truth(), thresholds, "t", "c", max_points=2
        )
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[-1]["rf"], 0, "the state is replayed in full even when skipped")
        self.assertTrue(any("not scored" in p for p in problems))


class FigureTests(unittest.TestCase):
    def test_every_configuration_has_one_fixed_distinct_style(self) -> None:
        styles = [aot.config_style(c.name) for c in common.CONFIGS]
        self.assertEqual(len({s[0] for s in styles}), len(styles), "colours must differ")
        self.assertEqual(len({s[1] for s in styles}), len(styles), "markers must differ")
        # the style follows the configuration, never its position in a list
        self.assertEqual(aot.config_style("flat"), aot.config_style("flat", ["zzz", "flat"]))
        extra = aot.config_style("custom_a", ["custom_b", "custom_a"])
        self.assertNotIn(extra[0], {s[0] for s in styles})

    def test_small_multiples_draw_empty_panels_with_their_reason(self) -> None:
        plt, _ = common.try_import_matplotlib()
        if plt is None:
            self.skipTest("matplotlib not importable")
        rows = [
            {"t_ms": 100.0, "rf_normalized": 0.9},
            {"t_ms": 300.0, "rf_normalized": 0.2},
            {"t_ms": 400.0, "rf_normalized": None},
        ]
        series = {"sim/baseline": rows, "sim/flat": rows, "real/baseline": rows}
        with tempfile.TemporaryDirectory() as tmp:
            png = Path(tmp) / "fig.png"
            written, reason = aot.plot_curves(
                series, png, title="t", empty={"trio": "3 genomes: RF undefined"},
                labels={"sim": "sim (10 simulated genomes)"},
            )
            self.assertTrue(written, reason)
            self.assertGreater(png.stat().st_size, 1000)
        panels = aot.panels_from_series(series, empty={"trio": "why"})
        self.assertEqual([p["dataset"] for p in panels], ["real", "sim", "trio"])
        self.assertEqual(panels[2]["empty_reason"], "why")
        self.assertEqual(set(panels[1]["series"]), {"baseline", "flat"})

    def test_nothing_to_plot_is_reported_not_drawn(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            written, reason = aot.plot_curves({"x/baseline": []}, Path(tmp) / "f.png", title="t")
        self.assertFalse(written)
        self.assertTrue(reason)

    def test_dataset_labels_say_what_the_dataset_is(self) -> None:
        summary = {
            "datasets": [{"name": "sim_100", "kind": "simulated", "genome_count": 100}],
            "cells": [{"dataset": "sim_100", "config": "baseline",
                       "run": {"total_length": 100_000_000, "genomes": 100}}],
        }
        label = aot.labels_from_summary(summary)["sim_100"]
        self.assertIn("100 simulated genomes", label)
        self.assertIn("1.0 Mb", label)


if __name__ == "__main__":
    unittest.main()
