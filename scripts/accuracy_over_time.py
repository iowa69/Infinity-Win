#!/usr/bin/env python3
"""Tree accuracy versus elapsed time, read out of ``events.jsonl``.

PROJECT.md section 11 names this the key figure: "tree accuracy plotted against
elapsed time, one curve per configuration. It shows the 'photo coming into focus'
as a measurement."

Accuracy at a given instant is the **Robinson-Foulds distance** between the tree
Infinity-Win had emitted by then and the truth tree, *after collapsing branches shorter
than half a SNP* (a branch that short carries no evidence, so the split it implies
must not be counted).

Where the trees come from
    ``events.jsonl`` carries one JSON object per line::

        {"seq": 12, "t_ms": 4310.5, "wall_ms": 1763..., 
         "event": {"event": "cluster_focus", "data": {...}}}

    The tree-emitting events are ``backbone_tree`` (leaves are clusters) and
    ``cluster_focus`` (one cluster at a focus level). This script replays them in
    ``t_ms`` order and, after each one, re-assembles the *full* tree exactly the
    way section 6.9 describes: take the backbone tree and replace each pivot leaf
    with its cluster tree. That assembled tree is what is compared with truth, so
    the curve shows the whole picture sharpening, not one cluster at a time.

Where the RF values come from (``--rf-mode``)
    ``auto``     use a ``bench.json`` timeline when one exists, else compute RF here
    ``bench``    require RF values from a bench JSON (``--rf-from``)
    ``internal`` always compute RF in pure Python from the emitted Newick
    ``infinity``    write each emitted tree to a temporary output folder and ask
                 ``infinity bench`` to score it (works only once the CLI exists)

The numbers are always written as CSV. The PNG is a bonus: when matplotlib is not
importable the script says so and still leaves the CSV behind.

Usage::

    scripts/accuracy_over_time.py --run-dir bench-work/matrix/ncbi_pd_cluster/baseline \
        --truth-tree bench-work/ncbi_pd_cluster/truth/cluster.treefile

    scripts/accuracy_over_time.py --matrix-dir bench-work/matrix --dataset ncbi_pd_cluster

    scripts/accuracy_over_time.py --run-dir <dir> --rf-from <dir>/bench.json
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

CSV_HEADER = [
    "dataset",
    "config",
    "t_ms",
    "t_s",
    "seq",
    "event",
    "cluster",
    "focus",
    "rf_source",
    "leaves_run",
    "leaves_truth",
    "leaves_shared",
    "splits_run",
    "splits_truth",
    "splits_shared",
    "rf",
    "rf_normalized",
    "collapse_threshold",
    "note",
]

#: Colour/marker cycle for the figure, in the order configurations are drawn.
STYLE_CYCLE = [
    ("#1f77b4", "o"),
    ("#d62728", "s"),
    ("#2ca02c", "^"),
    ("#ff7f0e", "D"),
    ("#9467bd", "v"),
    ("#8c564b", "P"),
    ("#e377c2", "X"),
    ("#7f7f7f", "*"),
]


# ---------------------------------------------------------------------------
# Replaying the event stream into full trees
# ---------------------------------------------------------------------------


@dataclass
class ClusterState:
    """What is known about one cluster while the stream is replayed."""

    cluster_id: str
    pivot_name: Optional[str] = None
    member_names: List[str] = field(default_factory=list)
    latest_newick: Optional[str] = None
    latest_focus: Optional[int] = None

    def aliases(self) -> List[str]:
        """Every name the backbone tree might use for this cluster's leaf."""
        out: List[str] = []
        if self.pivot_name:
            out.append(self.pivot_name)
        cid = str(self.cluster_id)
        out.extend([cid, f"c{cid}", f"cluster_{cid}", f"cluster{cid}", f"Cluster{cid}"])
        return out


def _as_name_list(value: Any) -> List[str]:
    if not isinstance(value, (list, tuple)):
        return []
    return [str(v) for v in value if v is not None]


def collect_clusters(log: common.EventLog) -> Dict[str, ClusterState]:
    """Cluster ids, pivots and member names from ``clustered``/``genome_scanned``."""
    clusters: Dict[str, ClusterState] = {}
    clustered = log.first("clustered")
    if clustered is not None:
        data = common.event_data(clustered)
        for entry in data.get("clusters") or []:
            if not isinstance(entry, dict):
                continue
            cid = entry.get("id", entry.get("cluster"))
            if cid is None:
                continue
            state = clusters.setdefault(str(cid), ClusterState(str(cid)))
            state.pivot_name = entry.get("pivot_name") or (
                str(entry.get("pivot")) if entry.get("pivot") is not None else None
            )
            names = _as_name_list(entry.get("member_names"))
            if names:
                state.member_names = names
    # genome_scanned gives an id -> name map, used when a tree is labelled by id.
    id_to_name: Dict[str, str] = {}
    for ev in log.of_type("genome_scanned"):
        data = common.event_data(ev)
        gid = data.get("id", data.get("genome_id", data.get("genome")))
        name = data.get("name", data.get("genome_name"))
        if gid is not None and name:
            id_to_name[str(gid)] = str(name)
    for state in clusters.values():
        if state.pivot_name in id_to_name:
            state.pivot_name = id_to_name[state.pivot_name]
        state.member_names = [id_to_name.get(m, m) for m in state.member_names]
    return clusters


def remap_leaf_names(tree: common.NewickNode, mapping: Dict[str, str]) -> int:
    """Rename leaves through ``mapping`` (genome id -> name). Returns the count."""
    changed = 0
    for node in tree.iter_nodes():
        if node.is_leaf() and node.name in mapping:
            node.name = mapping[node.name]
            changed += 1
    return changed


def graft(
    backbone: common.NewickNode,
    clusters: Dict[str, ClusterState],
    id_to_name: Optional[Dict[str, str]] = None,
) -> Tuple[common.NewickNode, List[str]]:
    """Replace each cluster leaf of ``backbone`` with that cluster's tree (6.9)."""
    notes: List[str] = []
    id_to_name = id_to_name or {}
    alias_to_cluster: Dict[str, str] = {}
    for cid, state in clusters.items():
        if state.latest_newick is None:
            continue
        for alias in state.aliases():
            alias_to_cluster.setdefault(alias.lower(), cid)
    if not alias_to_cluster:
        return backbone.copy(), ["no cluster tree had been emitted yet"]
    used: set = set()

    def parse_cluster(cid: str) -> Optional[common.NewickNode]:
        state = clusters[cid]
        try:
            sub = common.parse_newick(state.latest_newick or "")
        except common.PrereqError as exc:
            notes.append(f"cluster {cid}: unparsable Newick ({exc})")
            return None
        # A tree may be labelled by genome id rather than name; the truth tree is
        # labelled by name, so ids are translated before any comparison.
        remap_leaf_names(sub, id_to_name)
        return sub

    def build(src: common.NewickNode) -> common.NewickNode:
        if src.is_leaf():
            cid = alias_to_cluster.get((src.name or "").lower())
            if cid is not None and cid not in used:
                sub = parse_cluster(cid)
                if sub is not None:
                    used.add(cid)
                    if sub.length is None:
                        sub.length = src.length
                    return sub
                return common.NewickNode(src.name, src.length)
            return common.NewickNode(src.name, src.length)
        out = common.NewickNode("" if src is backbone else src.name, src.length)
        for child in src.children:
            out.add(build(child))
        return out

    grafted = build(backbone)
    # Clusters that the backbone did not mention (a single-cluster run, or a
    # backbone whose leaves are not cluster names) are attached at the root so
    # their genomes still count towards the accuracy of this instant.
    leftovers = [cid for cid in clusters if cid not in used and clusters[cid].latest_newick]
    if leftovers:
        notes.append(
            f"{len(leftovers)} cluster tree(s) could not be placed on the backbone "
            f"and were attached at the root: {sorted(leftovers)[:5]}"
        )
        for cid in leftovers:
            sub = parse_cluster(cid)
            if sub is not None:
                grafted.add(sub)
    return grafted, notes


def assemble(
    backbone_newick: Optional[str],
    clusters: Dict[str, ClusterState],
    id_to_name: Dict[str, str],
) -> Tuple[Optional[common.NewickNode], List[str]]:
    """The full tree as it looks right now, or ``(None, reasons)``."""
    notes: List[str] = []
    emitted = {cid: s for cid, s in clusters.items() if s.latest_newick}
    if backbone_newick:
        try:
            backbone = common.parse_newick(backbone_newick)
        except common.PrereqError as exc:
            return None, [f"backbone_tree Newick is unparsable: {exc}"]
        remap_leaf_names(backbone, id_to_name)
        tree, graft_notes = graft(backbone, clusters, id_to_name)
        notes.extend(graft_notes)
        return tree, notes
    if len(emitted) == 1:
        only = next(iter(emitted.values()))
        try:
            tree = common.parse_newick(only.latest_newick or "")
        except common.PrereqError as exc:
            return None, [f"cluster {only.cluster_id} Newick is unparsable: {exc}"]
        remap_leaf_names(tree, id_to_name)
        notes.append("no backbone tree yet: the single emitted cluster tree is used")
        return tree, notes
    if emitted:
        root = common.NewickNode()
        for cid in sorted(emitted):
            try:
                root.add(common.parse_newick(emitted[cid].latest_newick or ""))
            except common.PrereqError as exc:
                notes.append(f"cluster {cid} Newick is unparsable: {exc}")
        remap_leaf_names(root, id_to_name)
        notes.append(
            "no backbone tree yet: the emitted cluster trees are joined at a "
            "temporary root, so the topology between clusters is not meaningful"
        )
        return root, notes
    return None, ["no tree had been emitted yet"]


def tree_points(
    log: common.EventLog,
    *,
    max_points: Optional[int] = None,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """One entry per tree-emitting event, in ``t_ms`` order."""
    clusters = collect_clusters(log)
    id_to_name: Dict[str, str] = {}
    for ev in log.of_type("genome_scanned"):
        data = common.event_data(ev)
        gid = data.get("id", data.get("genome_id", data.get("genome")))
        name = data.get("name", data.get("genome_name"))
        if gid is not None and name:
            id_to_name[str(gid)] = str(name)

    emitting = [e for e in log.events if common.event_name(e) in ("backbone_tree", "cluster_focus")]
    emitting.sort(key=lambda e: (common.event_t_ms(e) if common.event_t_ms(e) is not None else 0.0,
                                 e.get("seq", 0)))
    if max_points and len(emitting) > max_points:
        # Keep the first and last, sample the middle evenly, so the shape of the
        # curve survives subsampling.
        step = (len(emitting) - 1) / (max_points - 1)
        picked = sorted({round(i * step) for i in range(max_points)})
        dropped = len(emitting) - len(picked)
        emitting = [emitting[i] for i in picked]
        note = f"{dropped} tree events were skipped to stay within --max-points {max_points}"
    else:
        note = ""

    backbone_newick: Optional[str] = None
    points: List[Dict[str, Any]] = []
    problems: List[str] = []
    if note:
        problems.append(note)
    if not emitting:
        problems.append(
            "events.jsonl contains no backbone_tree or cluster_focus event, so there "
            "is no tree to score"
        )
    for ev in emitting:
        name = common.event_name(ev)
        data = common.event_data(ev)
        t_ms = common.event_t_ms(ev)
        newick = data.get("newick")
        cluster_id = data.get("cluster", data.get("cluster_id"))
        focus = data.get("focus")
        if name == "backbone_tree":
            if not isinstance(newick, str) or not newick.strip():
                problems.append(f"seq {ev.get('seq')}: backbone_tree has no newick")
                continue
            backbone_newick = newick
        else:
            if cluster_id is None:
                problems.append(f"seq {ev.get('seq')}: cluster_focus has no cluster id")
                continue
            cid = str(cluster_id)
            state = clusters.setdefault(cid, ClusterState(cid))
            if not isinstance(newick, str) or not newick.strip():
                problems.append(f"seq {ev.get('seq')}: cluster_focus {cid} has no newick")
                continue
            if state.pivot_name is None:
                pivots = _as_name_list(data.get("pivots"))
                if pivots:
                    state.pivot_name = pivots[0]
            state.latest_newick = newick
            state.latest_focus = focus if isinstance(focus, int) else state.latest_focus
            names = data.get("members")
            if isinstance(names, list) and not state.member_names:
                state.member_names = [
                    str(m.get("name")) for m in names if isinstance(m, dict) and m.get("name")
                ]
        tree, notes = assemble(backbone_newick, clusters, id_to_name)
        if tree is None:
            problems.extend(notes)
            continue
        points.append(
            {
                "t_ms": t_ms,
                "seq": ev.get("seq"),
                "event": name,
                "cluster": None if cluster_id is None else str(cluster_id),
                "focus": focus if isinstance(focus, int) else None,
                "tree": tree,
                # Snapshot of every cluster tree emitted so far, as Newick text.
                # Scoring against a per-cluster truth (what `infinity simulate`
                # writes) needs the individual cluster trees, not just the
                # assembled full tree.
                "cluster_newicks": {
                    cid: state.latest_newick
                    for cid, state in clusters.items()
                    if state.latest_newick
                },
                "notes": notes,
            }
        )
    return points, problems


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def score_internal(
    points: List[Dict[str, Any]],
    truth: common.NewickNode,
    collapse: float,
    dataset: str,
    config: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for point in points:
        result = common.rf_distance(point["tree"], truth, collapse=collapse)
        rows.append(
            {
                "dataset": dataset,
                "config": config,
                **{k: point[k] for k in ("t_ms", "seq", "event", "cluster", "focus")},
                "rf_source": "internal RF",
                "rf": result.rf,
                "rf_normalized": result.rf_normalized,
                "leaves_run": result.leaves_a,
                "leaves_truth": result.leaves_b,
                "leaves_shared": result.leaves_shared,
                "splits_run": result.splits_a,
                "splits_truth": result.splits_b,
                "splits_shared": result.shared_splits,
                "collapse_threshold": collapse or None,
                "note": "; ".join(point.get("notes", []) + result.problems),
            }
        )
    return rows


def load_truth_trees(truth_dir: Path) -> Tuple[List[Tuple[str, common.NewickNode]], List[str]]:
    """Parse the per-cluster truth trees of a ``infinity simulate`` truth folder.

    The simulator writes one truth tree per cluster (``tree_c00.nwk``,
    ``tree_c01.nwk``, ...) in SNP units and **no** single full truth tree. The
    engine's own ``bench`` scores each emitted cluster tree against its matching
    truth tree, so this script does the same: comparing an assembled full tree
    against a made-up join of the per-cluster truths would invent inter-cluster
    splits that no truth file ever claimed.
    """
    problems: List[str] = []
    trees: List[Tuple[str, common.NewickNode]] = []
    if not truth_dir or not Path(truth_dir).is_dir():
        return trees, [f"no truth folder at {truth_dir}"]
    paths = sorted(Path(truth_dir).glob("tree_*.nwk"))
    if not paths:
        return trees, [f"no tree_*.nwk truth trees in {truth_dir}"]
    for path in paths:
        try:
            trees.append((path.stem, common.read_newick_file(path)))
        except common.PrereqError as exc:
            problems.append(str(exc))
    return trees, problems


def match_truth_tree(
    leaves: set, truth_trees: Sequence[Tuple[str, common.NewickNode]]
) -> Tuple[Optional[str], Optional[common.NewickNode], int]:
    """The truth tree sharing the most leaves with ``leaves``."""
    best_name = None
    best_tree = None
    best_overlap = 0
    for name, tree in truth_trees:
        overlap = len(leaves & set(tree.leaf_names()))
        if overlap > best_overlap:
            best_name, best_tree, best_overlap = name, tree, overlap
    return best_name, best_tree, best_overlap


def score_internal_per_cluster(
    points: List[Dict[str, Any]],
    truth_trees: Sequence[Tuple[str, common.NewickNode]],
    collapse: float,
    dataset: str,
    config: str,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Aggregate RF over the clusters emitted so far, at each instant.

    Mirrors what ``infinity bench`` reports in ``trees[]`` (one RF per cluster
    against that cluster's truth tree), but evaluated at *every* point of the
    event stream so the curve shows the picture sharpening. ``rf`` is the sum
    over clusters and ``rf_normalized`` their mean, both labelled as such in the
    CSV's ``note`` column.
    """
    rows: List[Dict[str, Any]] = []
    problems: List[str] = []
    for point in points:
        newicks = point.get("cluster_newicks") or {}
        rf_sum = 0
        norms: List[float] = []
        leaves_shared = 0
        splits_run = 0
        splits_truth = 0
        scored = 0
        unscored: List[str] = []
        for cid in sorted(newicks):
            try:
                tree = common.parse_newick(newicks[cid])
            except common.PrereqError as exc:
                unscored.append(f"cluster {cid}: unparsable ({exc})")
                continue
            leaves = set(tree.leaf_names())
            name, truth, overlap = match_truth_tree(leaves, truth_trees)
            if truth is None or overlap < 4:
                unscored.append(
                    f"cluster {cid}: only {overlap} shared leaves with the closest "
                    f"truth tree ({name}); RF needs at least 4"
                )
                continue
            result = common.rf_distance(tree, truth, collapse=collapse)
            if not result.ok:
                unscored.append(f"cluster {cid}: {result.problems[0] if result.problems else 'RF undefined'}")
                continue
            scored += 1
            rf_sum += result.rf or 0
            if result.rf_normalized is not None:
                norms.append(result.rf_normalized)
            leaves_shared += result.leaves_shared or 0
            splits_run += result.splits_a or 0
            splits_truth += result.splits_b or 0
        notes = list(point.get("notes") or [])
        notes.append(
            f"rf is the sum over {scored} cluster tree(s) and rf_normalized their mean"
        )
        if point.get("event") == "backbone_tree" and not newicks:
            notes.append("backbone only: no cluster resolved yet, so accuracy is undefined")
        for item in dict.fromkeys(unscored):
            notes.append(item)
        rows.append(
            {
                "dataset": dataset,
                "config": config,
                **{k: point[k] for k in ("t_ms", "seq", "event", "cluster", "focus")},
                "rf_source": "internal RF (per cluster)",
                "rf": rf_sum if scored else None,
                "rf_normalized": (sum(norms) / len(norms)) if norms else None,
                "leaves_run": None,
                "leaves_truth": None,
                "leaves_shared": leaves_shared if scored else None,
                "splits_run": splits_run if scored else None,
                "splits_truth": splits_truth if scored else None,
                "splits_shared": None,
                "collapse_threshold": collapse or None,
                "note": "; ".join(notes),
            }
        )
    if not any(r["rf"] is not None for r in rows):
        problems.append(
            "no cluster tree could be scored against the per-cluster truth trees; "
            "check that genome names match between the run and truth/tree_*.nwk"
        )
    return rows, problems


def score_with_infinity(
    run_dir: Path,
    truth_dir: Path,
    binary: Path,
    collapse_below: float,
    dataset: str,
    config: str,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Ask the CLI's own ``infinity bench`` for the RF timeline.

    ``infinity bench --output <run dir> --truth <truth dir>`` already replays the
    event stream and scores every emitted tree against the simulator truth, so
    this mode runs bench once and reads its ``timeline`` - no temporary folders.
    The curve then comes from the engine's own Robinson-Foulds code rather than
    from the Python reimplementation in this directory, which is the point of
    offering both.

    It requires a ``infinity simulate`` truth folder. Real datasets have a reference
    tree and a SKA2 table instead, which ``bench --truth`` does not accept; for
    those use ``--rf-mode internal`` with ``--truth-tree`` (one full reference
    tree) or ``--rf-mode bench`` against an existing ``bench.json``.
    """
    problems: List[str] = []
    if not common.is_simulator_truth_dir(truth_dir):
        raise common.PrereqError(
            f"--rf-mode infinity needs a `infinity simulate` truth folder; {truth_dir} is not "
            "one (expected tree_*.nwk, or pairwise_snps.tsv together with clusters.tsv)\n"
            "       for a real dataset use --rf-mode internal with --truth-tree <newick>, "
            "or --rf-mode bench with an existing bench.json"
        )
    out_json = run_dir / "bench.json"
    cmd = [
        str(binary), "bench",
        "--output", str(run_dir),
        "--truth", str(truth_dir),
        "--json", str(out_json),
        "--collapse-below", str(collapse_below),
        "--quiet",
    ]
    res = common.run_child(
        cmd,
        stdout_path=run_dir / "infinity.bench.log",
        stderr_path=run_dir / "infinity.bench.err",
        timeout=3600,
        poll_memory=False,
    )
    bench = common.load_bench_json(out_json, required=False)
    if bench is None:
        problems.append(
            f"`{' '.join(cmd)}` exited {res.returncode} and wrote no bench JSON: "
            f"{res.stderr_tail(300)}"
        )
        return [], problems
    rows, bench_problems = rows_from_bench(bench, dataset, config)
    problems.extend(bench_problems)
    for row in rows:
        row["rf_source"] = "infinity bench"
    return rows, problems


def rows_from_bench(bench: common.BenchJson, dataset: str, config: str) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Turn the RF timeline of a bench JSON into curve rows.

    Handles the schema the CLI actually writes, where each entry looks like::

        {"source": "cluster 0 (live)", "truth": "", "leaves": 5, "rf": 0,
         "splits_emitted": 0, "splits_true": 0, "rf_normalised": null,
         "t_ms": 24.9, "approximate": true, "focus": 0, "units": "subs_per_site"}

    ``rf_normalised`` is null for the earliest points: with fewer than four taxa
    in the emitted tree the RF distance is undefined, and the engine says so
    rather than reporting zero. That null is passed through untouched.
    """
    series = common.rf_values_from_bench(bench)
    rows: List[Dict[str, Any]] = []
    problems: List[str] = []
    if not series:
        problems.append(
            f"{bench.path.name} carries no RF values (looked for timeline / trees / "
            "tree.rf); the accuracy curve cannot be drawn from it"
        )
        return rows, problems
    for kind, entries in series.items():
        for entry in entries:
            t_ms = entry.get("t_ms", entry.get("elapsed_ms"))
            source = entry.get("source") or kind
            cluster = entry.get("cluster")
            if cluster is None and isinstance(source, str):
                match = re.search(r"cluster[_ ](\d+)", source)
                if match:
                    cluster = match.group(1)
            rf_norm = entry.get("rf_normalised", entry.get("rf_normalized", entry.get("rf_norm")))
            notes = []
            if entry.get("approximate"):
                notes.append("approximate (Focus 1 estimate)")
            if entry.get("units"):
                notes.append(f"branch units: {entry['units']}")
            if entry.get("truth"):
                notes.append(f"vs {entry['truth']}")
            if rf_norm is None and entry.get("rf") is not None:
                notes.append(
                    f"normalised RF undefined ({entry.get('leaves')} leaves in the "
                    "emitted tree); raw RF is kept"
                )
            rows.append(
                {
                    "dataset": dataset,
                    "config": config,
                    "t_ms": None if t_ms is None else float(t_ms),
                    "seq": entry.get("seq"),
                    "event": source,
                    "cluster": None if cluster is None else str(cluster),
                    "focus": entry.get("focus"),
                    "rf_source": f"bench.json ({kind})",
                    "rf": common._to_float(entry.get("rf")),
                    "rf_normalized": common._to_float(rf_norm),
                    "leaves_run": entry.get("leaves", entry.get("leaves_run")),
                    "leaves_truth": entry.get("leaves_truth"),
                    "leaves_shared": entry.get("leaves_shared", entry.get("leaves")),
                    "splits_run": entry.get("splits_emitted", entry.get("splits_run")),
                    "splits_truth": entry.get("splits_true", entry.get("splits_truth")),
                    "splits_shared": entry.get("splits_shared"),
                    "collapse_threshold": entry.get("collapse_threshold"),
                    "note": "; ".join(notes),
                }
            )
    if not rows:
        problems.append(f"{bench.path.name}: RF series was present but empty")
    return rows, problems


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_rows(csv_path: Path, rows: Sequence[Dict[str, Any]]) -> Path:
    """Write the accuracy curve as CSV. Each row carries its own dataset/config."""
    common.write_csv(
        csv_path,
        CSV_HEADER,
        [
            [
                r.get("dataset"),
                r.get("config"),
                r.get("t_ms"),
                None if r.get("t_ms") is None else round(float(r["t_ms"]) / 1000.0, 6),
                r.get("seq"),
                r.get("event"),
                r.get("cluster"),
                r.get("focus"),
                r.get("rf_source"),
                r.get("leaves_run"),
                r.get("leaves_truth"),
                r.get("leaves_shared"),
                r.get("splits_run"),
                r.get("splits_truth"),
                r.get("splits_shared"),
                r.get("rf"),
                r.get("rf_normalized"),
                r.get("collapse_threshold"),
                r.get("note"),
            ]
            for r in rows
        ],
    )
    return csv_path


def plot_curves(
    series: Dict[str, List[Dict[str, Any]]],
    png_path: Path,
    *,
    title: str,
    log_x: bool = True,
) -> Tuple[bool, Optional[str]]:
    """One curve per configuration. Returns ``(written, reason)``."""
    plt, reason = common.try_import_matplotlib()
    if plt is None:
        return False, reason
    drawable = {k: v for k, v in series.items() if v}
    if not drawable:
        return False, "no data points to plot"
    try:
        fig, ax = plt.subplots(figsize=(8.2, 5.0), dpi=150)
        for i, (label, rows) in enumerate(sorted(drawable.items())):
            colour, marker = STYLE_CYCLE[i % len(STYLE_CYCLE)]
            pts = sorted(
                ((float(r["t_ms"]) / 1000.0, r.get("rf_normalized")) for r in rows
                 if r.get("t_ms") is not None and r.get("rf_normalized") is not None),
                key=lambda p: p[0],
            )
            if not pts:
                continue
            xs = [max(p[0], 1e-3) for p in pts]
            ys = [p[1] for p in pts]
            ax.plot(xs, ys, color=colour, marker=marker, markersize=3.5, linewidth=1.4,
                    label=f"{label} ({len(pts)} pts)")
        ax.set_xscale("log" if log_x else "linear")
        ax.set_xlabel("elapsed time (s)")
        ax.set_ylabel("normalised Robinson-Foulds distance to truth")
        ax.set_ylim(-0.02, 1.05)
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.25, linewidth=0.6)
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
        fig.tight_layout()
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path)
        plt.close(fig)
    except Exception as exc:  # noqa: BLE001 - plotting must never lose the numbers
        return False, f"plotting failed ({type(exc).__name__}: {exc})"
    return True, None


# ---------------------------------------------------------------------------
# One run directory
# ---------------------------------------------------------------------------


def find_truth_tree(run_dir: Path, explicit: Optional[Path], dataset_dir: Optional[Path]) -> Tuple[Optional[Path], List[str]]:
    notes: List[str] = []
    if explicit:
        p = Path(explicit)
        if p.exists():
            return p, notes
        notes.append(f"--truth-tree {p} does not exist")
        return None, notes
    candidates: List[Path] = []
    if dataset_dir:
        candidates.append(Path(dataset_dir) / "truth" / "cluster.treefile")
        candidates.append(Path(dataset_dir) / "truth" / "truth.nwk")
    # run_dir is <matrix>/<dataset>/<config>: walk up to the work dir
    for parent in [run_dir.parent.parent, run_dir.parent]:
        candidates.append(parent / "truth" / "cluster.treefile")
        candidates.append(parent / "truth" / "truth.nwk")
        candidates.append(parent / "truth.nwk")
    for cand in candidates:
        if cand.exists():
            notes.append(f"truth tree found by convention at {cand}")
            return cand, notes
    notes.append(
        "no truth tree found; pass --truth-tree <newick> (or a dataset folder with "
        "truth/ inside bench-work/<dataset>/)"
    )
    return None, notes


def analyse_run_dir(
    run_dir: Path,
    *,
    dataset: str,
    config: str,
    out_dir: Path,
    truth_tree: Optional[Path],
    truth_dir: Optional[Path] = None,
    rf_from: Optional[Path],
    rf_mode: str,
    collapse_snps: float,
    collapse_mode: str,
    collapse_length: Optional[float],
    comparable_bp: Optional[int],
    max_points: Optional[int],
    infinity: Optional[Path],
) -> Dict[str, Any]:
    """Produce the accuracy-over-time rows for one output folder."""
    result: Dict[str, Any] = {
        "dataset": dataset,
        "config": config,
        "run_dir": str(run_dir),
        "status": "pending",
        "problems": [],
        "points": 0,
        "csv": None,
    }
    events_path = run_dir / "events.jsonl"
    if not events_path.exists():
        result["status"] = "failed"
        result["problems"].append(
            f"no events.jsonl in {run_dir}: this configuration produced no event "
            "stream, so accuracy over time cannot be measured"
        )
        return result
    log = common.read_event_log(events_path)
    result["problems"].extend(log.problems)
    result["events"] = {"read": log.lines_read, "bad": log.lines_bad, "used": len(log.events)}

    run_json = common.load_run_json(run_dir, required=False)

    # --- RF from a bench JSON ------------------------------------------
    bench_rows: List[Dict[str, Any]] = []
    if rf_from is not None or rf_mode in ("auto", "bench"):
        bench_path = Path(rf_from) if rf_from else run_dir / "bench.json"
        if bench_path.is_dir():
            bench_path = bench_path / "bench.json"
        if bench_path.exists():
            bench = common.load_bench_json(bench_path, required=False)
            if bench is None:
                result["problems"].append(f"{bench_path} is unreadable")
            else:
                rows, problems = rows_from_bench(bench, dataset, config)
                bench_rows = rows
                result["problems"].extend(problems)
                result["bench_json"] = str(bench_path)
        elif rf_mode == "bench":
            result["status"] = "failed"
            result["problems"].append(f"--rf-mode bench needs a bench JSON; none at {bench_path}")
            return result

    # --- collapse threshold --------------------------------------------
    if comparable_bp is None:
        comparable_bp = common.comparable_bp_for(run_dir, run_json)
    truth_node: Optional[common.NewickNode] = None
    if truth_tree and Path(truth_tree).exists():
        truth_node = common.read_newick_file(Path(truth_tree))
    threshold, explanation = common.guess_collapse_threshold(
        truth_node,
        collapse_snps=collapse_snps,
        comparable_bp=comparable_bp,
        mode=collapse_mode,
        explicit_length=collapse_length,
    )
    result["collapse"] = {
        "threshold": threshold,
        "explanation": explanation,
        "comparable_bp": comparable_bp,
    }

    # --- rows -----------------------------------------------------------
    # Which RF source gives the best accuracy-over-time curve:
    #
    #   an explicit --rf-mode   that mode, and nothing else
    #   auto                    1. a bench.json series that has RF at more than
    #                               one instant - the engine's own numbers
    #                            2. RF recomputed here at every emitted tree,
    #                               against per-cluster truth (simulated data) or
    #                               one reference tree (real data). Today this is
    #                               the only source with real time resolution: the
    #                               engine scores its live trees against nothing,
    #                               and its final per-cluster trees all share one
    #                               timestamp.
    #                            3. the bench.json series anyway, labelled,
    #                               because one point beats none
    rows: List[Dict[str, Any]] = []
    source: Optional[str] = None

    def defined_instants(series: Sequence[Dict[str, Any]]) -> int:
        return len(
            {
                round(float(r["t_ms"]), 6)
                for r in series
                if r.get("rf_normalized") is not None and r.get("t_ms") is not None
            }
        )

    # --- internal (pure Python) RF, computed once and then selected ------
    internal_rows: List[Dict[str, Any]] = []
    internal_source: Optional[str] = None
    cluster_threshold = threshold
    if rf_mode in ("auto", "internal"):
        points, problems = tree_points(log, max_points=max_points)
        result["problems"].extend(problems)
        # Simulated truth is a set of per-cluster trees in SNP units; real truth
        # is one reference tree in substitutions per site. Which one exists
        # decides how RF is computed and what "half a SNP" means.
        per_cluster_truth: List[Tuple[str, common.NewickNode]] = []
        if truth_dir and common.is_simulator_truth_dir(Path(truth_dir)):
            per_cluster_truth, tc_problems = load_truth_trees(Path(truth_dir))
            result["problems"].extend(tc_problems)
        if per_cluster_truth:
            cluster_threshold, cluster_why = common.guess_collapse_threshold(
                per_cluster_truth[0][1],
                collapse_snps=collapse_snps,
                comparable_bp=None,
                mode=collapse_mode,
                explicit_length=collapse_length,
            )
            result["collapse"]["per_cluster_threshold"] = cluster_threshold
            result["collapse"]["per_cluster_explanation"] = cluster_why
        internal_reason: Optional[str] = None
        if not points:
            internal_reason = "no tree-emitting events were found in the stream"
        elif per_cluster_truth:
            internal_rows, problems = score_internal_per_cluster(
                points, per_cluster_truth, cluster_threshold, dataset, config
            )
            result["problems"].extend(problems)
            internal_source = "internal RF (per cluster)"
        elif truth_node is not None:
            internal_rows = score_internal(points, truth_node, threshold, dataset, config)
            internal_source = "internal RF"
        else:
            internal_reason = (
                "no truth tree (--truth-tree) and no simulator truth folder "
                "(--truth-dir), so RF cannot be recomputed here"
            )
        if internal_reason:
            result["problems"].append(internal_reason)

    if rf_mode == "bench":
        rows, source = bench_rows, "bench.json"
    elif rf_mode == "infinity":
        binary = infinity or common.resolve_binary(None, required=False)
        if binary is None:
            raise common.PrereqError(
                "--rf-mode infinity needs the infinity binary; none was found. "
                "Use --rf-mode internal (pure Python RF) or --rf-from <bench.json>."
            )
        if not (truth_dir and common.is_simulator_truth_dir(Path(truth_dir))):
            raise common.PrereqError(
                "--rf-mode infinity needs a `infinity simulate` truth folder; got "
                f"{truth_dir or 'nothing'}. Pass --truth-dir <truth folder>, or use "
                "--rf-mode internal with --truth-tree <newick>."
            )
        rows, problems = score_with_infinity(
            run_dir, Path(truth_dir), Path(binary), collapse_snps, dataset, config
        )
        result["problems"].extend(problems)
        result["bench_json"] = str(run_dir / "bench.json")
        source = "infinity bench" if rows else None
    elif rf_mode == "internal":
        rows, source = internal_rows, internal_source
    else:  # auto
        instants = defined_instants(bench_rows)
        if bench_rows and instants >= 2:
            rows, source = bench_rows, "bench.json"
            if internal_rows:
                result["problems"].append(
                    f"RF was also recomputed here ({len(internal_rows)} points); the "
                    "engine's own bench values were preferred because they cover "
                    f"{instants} distinct instants. Use --rf-mode internal for the "
                    "recomputed curve."
                )
        elif internal_rows:
            rows, source = internal_rows, internal_source
            if bench_rows:
                result["problems"].append(
                    f"the bench.json RF series has {len(bench_rows)} entries but only "
                    f"{instants} distinct timestamp(s) with a defined normalised RF, "
                    "which cannot show accuracy growing over time; RF was recomputed "
                    "here from the emitted Newick instead. Use --rf-mode bench to "
                    "force the bench values."
                )
        elif bench_rows:
            rows, source = bench_rows, "bench.json"
            result["problems"].append(
                f"only the bench.json series was available and it has {instants} "
                "distinct timestamp(s) with a defined normalised RF, so this curve has "
                "little or no time resolution"
            )

    result["points"] = len(rows)
    result["rf_source"] = source
    csv_path = out_dir / f"accuracy_over_time__{dataset}__{config}.csv"
    write_rows(csv_path, rows)
    result["csv"] = str(csv_path)
    result["rows"] = rows
    if rows:
        result["status"] = "ok"
    else:
        # An explicitly requested RF source that produced nothing is a failure;
        # `auto` producing nothing is an empty result with the reasons recorded.
        result["status"] = "failed" if rf_mode in ("bench", "internal", "infinity") else "empty"
        result["problems"].append(
            "no accuracy points could be produced; the CSV was written with a header only"
        )
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="accuracy_over_time.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--run-dir", type=Path, help="one infinity output folder (has events.jsonl)")
    src.add_argument("--matrix-dir", type=Path, help="a matrix folder: <dataset>/<config>/")
    p.add_argument("--dataset", default=None, help="dataset name (default: inferred/all)")
    p.add_argument("--configs", nargs="+", default=None, help="restrict to these configurations")
    p.add_argument("--truth-tree", type=Path, default=None, help="one full truth Newick tree (for --rf-mode internal)")
    p.add_argument(
        "--truth-dir",
        type=Path,
        default=None,
        help="a `infinity simulate` truth folder (tree_*.nwk, pairwise_snps.tsv); "
        "lets --rf-mode infinity/auto ask the CLI's own bench for the RF timeline "
        "(default: the dataset's truth.dir from dataset.json)",
    )
    p.add_argument("--work-dir", type=Path, default=None, help="bench-work folder (for truth lookup)")
    p.add_argument("--out-dir", type=Path, default=None, help="where CSV/PNG go (default: matrix dir or run dir)")
    p.add_argument(
        "--rf-from",
        type=Path,
        default=None,
        help="read RF values from this bench JSON (file or folder) instead of recomputing",
    )
    p.add_argument(
        "--rf-mode",
        choices=("auto", "bench", "internal", "infinity"),
        default="auto",
        help="where RF values come from (default: auto = bench.json if present, else internal)",
    )
    p.add_argument("--infinity", default=None, help="infinity binary, for --rf-mode infinity")
    p.add_argument(
        "--collapse-mode",
        choices=("auto", "snps", "length", "none"),
        default="auto",
        help="how to convert 'half a SNP' into a branch-length threshold (default auto)",
    )
    p.add_argument("--collapse-snps", type=float, default=0.5, help="SNP threshold (default 0.5)")
    p.add_argument("--collapse-length", type=float, default=None, help="absolute length threshold")
    p.add_argument(
        "--comparable-bp",
        type=int,
        default=None,
        help="comparable sites, to convert SNPs into substitutions/site "
        "(default: from cluster_*.dist.tsv or run.json total_length)",
    )
    p.add_argument("--max-points", type=int, default=400, help="cap scored points per run (default 400)")
    p.add_argument("--no-plot", action="store_true", help="write CSV only")
    p.add_argument("--linear-x", action="store_true", help="linear instead of log time axis")
    p.add_argument("--json", type=Path, default=None, help="write a status JSON here")
    p.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    return p


def discover_run_dirs(matrix_dir: Path, dataset: Optional[str], configs: Optional[Sequence[str]]) -> List[Tuple[str, str, Path]]:
    """``(dataset, config, run_dir)`` for every folder holding an events.jsonl."""
    if not matrix_dir.is_dir():
        raise common.PrereqError(f"matrix directory does not exist: {matrix_dir}")
    out: List[Tuple[str, str, Path]] = []
    for dataset_dir in sorted(p for p in matrix_dir.iterdir() if p.is_dir()):
        dname = dataset_dir.name
        if dataset and dname != dataset:
            continue
        for config_dir in sorted(p for p in dataset_dir.iterdir() if p.is_dir()):
            if configs and config_dir.name not in configs:
                continue
            if (config_dir / "events.jsonl").exists() or (config_dir / "run.json").exists():
                out.append((dname, config_dir.name, config_dir))
    if not out:
        raise common.PrereqError(
            f"no run folders with events.jsonl under {matrix_dir}"
            + (f" for dataset {dataset!r}" if dataset else "")
            + (f" and configs {list(configs)}" if configs else "")
            + "\n       run scripts/run_matrix.py first"
        )
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()

    targets: List[Tuple[str, str, Path]] = []
    if args.run_dir:
        run_dir = Path(args.run_dir)
        if not run_dir.is_dir():
            raise common.PrereqError(f"--run-dir is not a directory: {run_dir}")
        dname = args.dataset or run_dir.parent.name
        cname = run_dir.name
        if run_dir.name in {c.name for c in common.CONFIGS}:
            dname = args.dataset or run_dir.parent.name
        targets.append((dname, cname, run_dir))
        out_dir = Path(args.out_dir) if args.out_dir else run_dir
    else:
        matrix_dir = Path(args.matrix_dir)
        targets = discover_run_dirs(matrix_dir, args.dataset, args.configs)
        out_dir = Path(args.out_dir) if args.out_dir else matrix_dir / "figures"
    common.ensure_dir(out_dir)

    results: List[Dict[str, Any]] = []
    all_series: Dict[str, List[Dict[str, Any]]] = {}
    all_rows: List[Dict[str, Any]] = []
    for dataset, config, run_dir in targets:
        dataset_dir = work / dataset
        truth_tree = args.truth_tree
        truth_dir = args.truth_dir
        truth_notes: List[str] = []
        # The dataset folder knows both kinds of truth: a simulator truth folder
        # (per-cluster trees + pairwise SNPs) and, for real data, a reference tree.
        if dataset_dir.is_dir():
            obj, _ = common.load_json_optional(dataset_dir / "dataset.json")
            if isinstance(obj, dict):
                truth = obj.get("truth") or {}
                if truth_dir is None and truth.get("dir"):
                    truth_dir = Path(truth["dir"])
        if truth_tree is None and dataset_dir.is_dir():
            truth_tree, truth_notes = find_truth_tree(
                run_dir, None, dataset_dir if dataset_dir.is_dir() else None
            )
        elif truth_tree is None:
            truth_tree, truth_notes = find_truth_tree(run_dir, None, None)
        print(f"[{dataset}/{config}] {run_dir}")
        for note in truth_notes:
            print(f"    note: {note}")
        try:
            res = analyse_run_dir(
                run_dir,
                dataset=dataset,
                config=config,
                out_dir=out_dir,
                truth_tree=Path(truth_tree) if truth_tree else None,
                truth_dir=Path(truth_dir) if truth_dir else None,
                rf_from=args.rf_from,
                rf_mode=args.rf_mode,
                collapse_snps=args.collapse_snps,
                collapse_mode=args.collapse_mode,
                collapse_length=args.collapse_length,
                comparable_bp=args.comparable_bp,
                max_points=args.max_points,
                infinity=Path(args.infinity) if args.infinity else None,
            )
        except common.PrereqError as exc:
            res = {
                "dataset": dataset,
                "config": config,
                "run_dir": str(run_dir),
                "status": "failed",
                "problems": [str(exc)],
                "points": 0,
                "rows": [],
            }
        rows = res.pop("rows", [])
        res["problems"] = list(dict.fromkeys(res.get("problems", [])))
        results.append(res)
        print(f"    status={res['status']} points={res.get('points', 0)} "
              f"rf_source={res.get('rf_source')} collapse={common.fmt((res.get('collapse') or {}).get('threshold'), '.3g')}")
        for problem in res["problems"][:6]:
            print(f"    NOTE: {problem}")
        if len(res["problems"]) > 6:
            print(f"    ... and {len(res['problems']) - 6} more note(s)")
        if res.get("csv"):
            print(f"    csv: {res['csv']}")
        if rows:
            all_series[f"{dataset}/{config}" if len(targets) > 1 else config] = rows
            all_rows.extend(
                [{**r, "dataset": dataset, "config": config} for r in rows]
            )

    # combined CSV + figure. The CSV is written first and always: a missing
    # matplotlib must never cost the numbers.
    combined_csv = out_dir / "accuracy_over_time.csv"
    write_rows(combined_csv, all_rows)
    print(f"\ncombined CSV: {combined_csv} ({len(all_rows)} rows)")

    png_written: Dict[str, Any] = {}
    if not args.no_plot and all_series:
        title = "Tree accuracy over time"
        if len(targets) == 1:
            title += f" - {targets[0][0]}/{targets[0][1]}"
        png = out_dir / "accuracy_over_time.png"
        written, reason = plot_curves(all_series, png, title=title, log_x=not args.linear_x)
        if written:
            print(f"figure: {png}")
            png_written["combined"] = str(png)
        else:
            print(f"NOTE: no figure written - {reason}; the CSV above holds every number")
            png_written["error"] = reason
        # one figure per dataset when several were processed
        by_dataset: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
        for label, rows in all_series.items():
            d = label.split("/")[0] if "/" in label else (targets[0][0] if targets else "dataset")
            c = label.split("/", 1)[1] if "/" in label else label
            by_dataset.setdefault(d, {})[c] = rows
        if len(by_dataset) > 1:
            for d, series in sorted(by_dataset.items()):
                png = out_dir / f"accuracy_over_time__{d}.png"
                written, reason = plot_curves(
                    series, png, title=f"Tree accuracy over time - {d}", log_x=not args.linear_x
                )
                if written:
                    print(f"figure: {png}")
                    png_written[d] = str(png)
                else:
                    print(f"NOTE: no figure for {d} - {reason}")
                    png_written[d] = f"error: {reason}"

    if args.json:
        common.dump_json(
            args.json,
            {
                "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "out_dir": str(out_dir),
                "combined_csv": str(combined_csv),
                "rows": len(all_rows),
                "figures": png_written,
                "runs": results,
            },
        )
        print(f"status JSON: {args.json}")

    ok = sum(1 for r in results if r.get("points"))
    if not ok:
        print(
            "\nno accuracy points were produced for any run; see the notes above "
            "(a missing truth tree or a missing events.jsonl are the usual causes)",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    common.guard(main)()
