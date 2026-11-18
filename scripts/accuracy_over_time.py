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
    ``auto``     compute RF here whenever a truth exists, else use a ``bench.json``
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
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

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
    "leaves_exact",
    "clusters_exact",
    "rf",
    "rf_normalized",
    "rf_between_clusters",
    "rf_within_clusters",
    "collapse_threshold",
    "note",
]

#: One fixed style per configuration: (colour, marker, line style). The same
#: configuration has the same look in every panel and every figure, so one
#: legend serves them all. Colours are slots 1-6 of the dataviz reference
#: categorical palette (validated for adjacent-pair colour-vision-deficiency
#: separation on a light surface); marker and dash pattern are a second,
#: colour-free encoding of the same identity.
CONFIG_STYLES: Dict[str, Tuple[str, str, Any]] = {
    "baseline": ("#2a78d6", "o", "-"),
    "no_block_skip": ("#eb6834", "s", (0, (5, 2))),
    "flat": ("#1baf7a", "^", (0, (1, 1.5))),
    "no_piece_cache": ("#eda100", "D", (0, (6, 2, 1, 2))),
    "no_profile_grouping": ("#e87ba4", "v", (0, (3, 1.5, 1, 1.5, 1, 1.5))),
    "single_cluster": ("#008300", "P", (0, (8, 2))),
}
#: Configurations outside the six of section 11, in sorted-name order. Past
#: these, a configuration is drawn grey rather than given a generated hue.
EXTRA_STYLES: List[Tuple[str, str, Any]] = [
    ("#6250d6", "X", "-"),
    ("#e34948", "*", (0, (5, 2))),
]
FALLBACK_STYLE: Tuple[str, str, Any] = ("#8a8984", ".", (0, (1, 1)))
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"

#: How the figure defines accuracy; printed under the figure and in REPORT.md.
ACCURACY_DEFINITION = (
    "Each point scores the whole tree on screen at that instant against the truth, "
    "over every genome of the dataset: clusters not yet refined are unresolved stars "
    "of their members, refined clusters are their emitted trees, all joined by the "
    "backbone. Normalised Robinson-Foulds = (splits in one tree but not the other) / "
    "(splits of the two trees added together); 1 = no correct split, 0 = identical "
    "topology, after collapsing branches shorter than half a SNP. Simulated truth = the "
    "per-cluster truth trees joined at one star, because every cluster root is "
    "derived independently from the reference."
)


def config_style(config: str, known: Sequence[str] = ()) -> Tuple[str, str, Any]:
    """The fixed (colour, marker, line style) of one configuration."""
    if config in CONFIG_STYLES:
        return CONFIG_STYLES[config]
    extras = sorted({c for c in known if c not in CONFIG_STYLES} | {config})
    position = extras.index(config)
    return EXTRA_STYLES[position] if position < len(EXTRA_STYLES) else FALLBACK_STYLE


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


def alias_map(clusters: Dict[str, "ClusterState"], only_emitted: bool = False) -> Dict[str, str]:
    """``lower-cased backbone leaf name -> cluster id``.

    Pivot names are claimed first; the generic spellings (``3``, ``c3``,
    ``cluster_3``...) only fill in afterwards, so a genome that happens to be
    called ``c1`` is never taken for cluster 1 when it is another cluster's pivot.
    """
    out: Dict[str, str] = {}
    chosen = {
        cid: state for cid, state in clusters.items()
        if not only_emitted or state.latest_newick is not None
    }
    for cid, state in chosen.items():
        if state.pivot_name:
            out.setdefault(state.pivot_name.lower(), cid)
    for cid, state in chosen.items():
        for alias in state.aliases():
            out.setdefault(alias.lower(), cid)
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
    alias_to_cluster = alias_map(clusters, only_emitted=True)
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


FULL_LEAF_SET_SOURCE = "internal RF (full leaf set)"


@dataclass
class FullTruth:
    """The truth tree over every genome, ready to be scored against."""

    root: common.NewickNode
    leaves: frozenset
    kind: str  # "simulator" or "reference"
    threshold: float
    explanation: str
    source: str


@dataclass
class RunThresholds:
    """'Half a SNP' in the units of each part of the tree the program shows.

    The backbone and Focus 0 cluster trees are in substitutions per site, Focus 1
    cluster trees in difference pieces (PROJECT.md 6.9: the Focus 1 estimate
    counts differing pieces, roughly one per SNP) and Focus 2 cluster trees in
    SNPs, so one number cannot collapse all of them.
    """

    per_site: float
    pieces: float
    snps: float
    explanation: str

    def for_focus(self, focus: Optional[int]) -> float:
        if focus == 2:
            return self.snps
        if focus == 1:
            return self.pieces
        return self.per_site


def load_full_truth(
    truth_dir: Optional[Path],
    truth_tree: Optional[Path],
    *,
    collapse_snps: float,
    collapse_mode: str,
    collapse_length: Optional[float],
    comparable_bp: Optional[int],
) -> Tuple[Optional[FullTruth], List[str]]:
    """One truth tree over the full leaf set, or ``(None, reasons)``.

    A simulator truth folder holds one tree per cluster (in SNPs, rooted at the
    cluster root) and no full tree. Every cluster root is the reference with its
    own independent substitutions (``infinity-sim``: "each is the reference with
    `divergence` of its positions substituted"), so the full truth is those
    trees joined at one star: each cluster is a clade, and nothing is claimed
    about the order in which clusters split from the reference. A real dataset
    has one reference tree instead, used as it is.
    """
    problems: List[str] = []
    if truth_dir and common.is_simulator_truth_dir(Path(truth_dir)):
        trees, tree_problems = load_truth_trees(Path(truth_dir))
        problems.extend(tree_problems)
        if trees:
            root = common.NewickNode()
            for _, tree in trees:
                tree.length = None  # the long branch from the reference to the cluster root
                root.add(tree)
            if collapse_mode == "none":
                threshold, why = 0.0, "collapsing disabled (--collapse-mode none)"
            elif collapse_mode == "length":
                threshold, why = common.guess_collapse_threshold(
                    None, mode="length", explicit_length=collapse_length
                )
            else:
                threshold = float(collapse_snps)
                why = (
                    f"simulator truth trees are in SNPs; collapsing branches shorter "
                    f"than {collapse_snps} SNP"
                )
            return (
                FullTruth(
                    root=root,
                    leaves=frozenset(root.leaf_names()),
                    kind="simulator",
                    threshold=threshold,
                    explanation=why,
                    source=f"{len(trees)} per-cluster truth tree(s) in {truth_dir}, joined at a star",
                ),
                problems,
            )
    if truth_tree and Path(truth_tree).exists():
        node = common.read_newick_file(Path(truth_tree))
        threshold, why = common.guess_collapse_threshold(
            node,
            collapse_snps=collapse_snps,
            comparable_bp=comparable_bp,
            mode=collapse_mode,
            explicit_length=collapse_length,
        )
        return (
            FullTruth(
                root=node,
                leaves=frozenset(node.leaf_names()),
                kind="reference",
                threshold=threshold,
                explanation=why,
                source=str(truth_tree),
            ),
            problems,
        )
    problems.append(
        "no truth tree (--truth-tree) and no simulator truth folder (--truth-dir), "
        "so RF cannot be recomputed here"
    )
    return None, problems


def run_thresholds(
    *,
    collapse_snps: float,
    collapse_mode: str,
    collapse_length: Optional[float],
    comparable_bp: Optional[int],
) -> RunThresholds:
    if collapse_mode == "none":
        return RunThresholds(0.0, 0.0, 0.0, "collapsing disabled (--collapse-mode none)")
    if collapse_mode == "length":
        if collapse_length is None:
            raise common.PrereqError("--collapse-mode length needs --collapse-length <float>")
        x = float(collapse_length)
        return RunThresholds(x, x, x, f"every branch shorter than {x} collapsed")
    snp = float(collapse_snps)
    if comparable_bp:
        per_site = snp / float(comparable_bp)
        site_note = f"{per_site:.3g} subs/site ({snp} SNP over {comparable_bp} comparable bp)"
    else:
        per_site = 0.0
        site_note = "not collapsed (no comparable-site count to convert SNPs into subs/site)"
    return RunThresholds(
        per_site=per_site,
        pieces=snp,
        snps=snp,
        explanation=(
            f"Focus 2 cluster trees: < {snp} SNP; Focus 1: < {snp} difference piece; "
            f"backbone and Focus 0: {site_note}"
        ),
    )


def _genome_universe(log: common.EventLog, id_to_name: Dict[str, str]) -> List[str]:
    """Every genome name the run knows about, from the whole event stream."""
    names: List[str] = list(id_to_name.values())
    for ev in log.of_type("clustered"):
        for entry in common.event_data(ev).get("clusters") or []:
            if isinstance(entry, dict):
                names.extend(id_to_name.get(str(m), str(m)) for m in _as_name_list(entry.get("member_names")))
    for ev in log.of_type("cluster_focus"):
        for member in common.event_data(ev).get("members") or []:
            if isinstance(member, dict) and member.get("name"):
                names.append(str(member["name"]))
    return sorted(set(names))


def score_full_leaf_set(
    log: common.EventLog,
    truth: FullTruth,
    thresholds: RunThresholds,
    dataset: str,
    config: str,
    *,
    max_points: Optional[int] = None,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """RF of the whole displayed tree against the truth, at every tree event.

    At each ``backbone_tree`` / ``cluster_focus`` event the tree the program is
    showing is rebuilt the way section 6.9 describes - the backbone, with each
    pivot leaf replaced by its cluster - where a cluster that has emitted no tree
    yet is an unresolved star of its members and every other cluster is its
    latest emitted tree. That tree is scored against the truth over the *full*
    leaf set: every genome of the run that the truth tree knows. A genome not
    yet placed anywhere counts as unresolved, never as missing. So the curve
    starts near 1 (a star resolves nothing) and falls as clusters become exact.

    Splits are bitmasks over the scored genomes; each cluster tree's splits are
    computed once when it is emitted, so a point costs one set union.
    """
    problems: List[str] = []
    rows: List[Dict[str, Any]] = []
    id_to_name: Dict[str, str] = {}
    for ev in log.of_type("genome_scanned"):
        data = common.event_data(ev)
        gid = data.get("id", data.get("genome_id", data.get("genome")))
        name = data.get("name", data.get("genome_name"))
        if gid is not None and name:
            id_to_name[str(gid)] = str(name)
    universe_names = _genome_universe(log, id_to_name)
    if not universe_names:
        problems.append("the event stream names no genomes, so there is nothing to score")
        return rows, problems
    scored_names = sorted(set(universe_names) & truth.leaves)
    not_in_truth = len(universe_names) - len(scored_names)
    if not_in_truth:
        problems.append(
            f"{not_in_truth} of {len(universe_names)} genomes are not in the truth tree "
            "and are left out of the score"
        )
    extra_truth = len(truth.leaves - set(universe_names))
    if extra_truth:
        problems.append(
            f"the truth tree has {extra_truth} leaves that are not genomes of this run "
            "(e.g. reference/outgroup genomes); they are left out of the score"
        )
    index = common.leaf_index(scored_names)
    universe = (1 << len(index)) - 1
    n_scored = len(index)
    defined = n_scored >= 4
    if not defined:
        problems.append(
            f"only {n_scored} genomes are both in the run and in the truth tree: "
            "Robinson-Foulds is undefined for fewer than 4 taxa"
        )
    _, truth_masks = common.clade_masks(truth.root, index, truth.threshold)
    truth_splits = common.canonical_splits(truth_masks, universe)
    # For simulated truth, each differing split is classed as "between clusters"
    # (it separates whole truth clusters: the backbone resolving what the truth
    # leaves as a star) or "within" (it cuts a cluster: a wrong cluster tree, or
    # a cluster tree rooted somewhere other than its true root, which changes
    # where the cluster hangs off the rest of the tree).
    truth_cluster_masks: List[int] = []
    if truth.kind == "simulator":
        for child in truth.root.children:
            mask, _ = common.clade_masks(child, index)
            if mask:
                truth_cluster_masks.append(mask)

    def between_clusters(split: int) -> bool:
        return all((split & cm) in (0, cm) for cm in truth_cluster_masks)

    clusters = collect_clusters(log)
    cache: Dict[str, Dict[str, Any]] = {}

    def member_mask(cid: str) -> int:
        entry = cache.get(cid)
        if entry is not None:
            return entry["leaf_mask"]
        mask = 0
        for name in clusters[cid].member_names:
            mask |= index.get(name, 0)
        return mask

    def emit_cluster_tree(cid: str, newick: str, focus: Optional[int]) -> None:
        try:
            tree = common.parse_newick(newick)
        except common.PrereqError as exc:
            problems.append(f"cluster {cid}: unparsable Newick ({exc})")
            return
        remap_leaf_names(tree, id_to_name)
        leaf_mask, masks = common.clade_masks(tree, index, thresholds.for_focus(focus))
        cache[cid] = {
            "leaf_mask": leaf_mask,
            "splits": common.canonical_splits(masks, universe),
            "focus": focus,
        }
        state = clusters.setdefault(cid, ClusterState(cid))
        if not state.member_names:
            state.member_names = sorted(n for n in tree.leaf_names() if n)

    backbone: Optional[common.NewickNode] = None

    def displayed_splits() -> Tuple[set, int, List[str]]:
        notes: List[str] = []
        clade_edges: List[int] = []
        present = 0
        used: set = set()
        alias_to_cluster = alias_map(clusters)
        if backbone is not None:
            # post-order over the backbone; a pivot leaf expands to its cluster
            value: Dict[int, int] = {}
            stack: List[Tuple[common.NewickNode, bool]] = [(backbone, False)]
            while stack:
                node, expanded = stack.pop()
                short = bool(
                    thresholds.per_site and node.length is not None
                    and node.length < thresholds.per_site
                )
                if not node.children:
                    cid = alias_to_cluster.get((node.name or "").lower())
                    if cid is not None and cid not in used:
                        used.add(cid)
                        mask = member_mask(cid)
                        if not short:
                            clade_edges.append(mask)  # the cluster as a clade
                    else:
                        mask = index.get(node.name, 0) if node.name else 0
                    value[id(node)] = mask
                    continue
                if not expanded:
                    stack.append((node, True))
                    stack.extend((child, False) for child in node.children)
                    continue
                mask = 0
                for child in node.children:
                    mask |= value.pop(id(child), 0)
                value[id(node)] = mask
                if node is not backbone and not short:
                    clade_edges.append(mask)
            present |= value.get(id(backbone), 0)
        else:
            notes.append(
                "no backbone tree yet: clusters are joined at one star, so no split "
                "between clusters is claimed"
            )
        leftovers = [cid for cid in sorted(clusters) if cid not in used]
        if backbone is not None and leftovers:
            notes.append(
                f"{len(leftovers)} cluster(s) are not on the backbone and were attached "
                "at the root"
            )
        for cid in leftovers:
            mask = member_mask(cid)
            clade_edges.append(mask)
            present |= mask
        found = common.canonical_splits(clade_edges, universe)
        for cid in used | set(leftovers):
            entry = cache.get(cid)
            if entry is not None:
                found |= entry["splits"]
        return found, present, notes

    emitting = [e for e in log.events if common.event_name(e) in ("backbone_tree", "cluster_focus")]
    emitting.sort(
        key=lambda e: (
            common.event_t_ms(e) if common.event_t_ms(e) is not None else 0.0,
            e.get("seq", 0),
        )
    )
    if not emitting:
        problems.append(
            "events.jsonl contains no backbone_tree or cluster_focus event, so there "
            "is no tree to score"
        )
        return rows, problems
    picked = set(range(len(emitting)))
    if max_points and len(emitting) > max_points:
        step = (len(emitting) - 1) / (max_points - 1)
        picked = {round(i * step) for i in range(max_points)}
        problems.append(
            f"{len(emitting) - len(picked)} of {len(emitting)} tree events were replayed "
            f"but not scored, to stay within --max-points {max_points}"
        )
    for position, ev in enumerate(emitting):
        name = common.event_name(ev)
        data = common.event_data(ev)
        newick = data.get("newick")
        cluster_id = data.get("cluster", data.get("cluster_id"))
        focus = data.get("focus") if isinstance(data.get("focus"), int) else None
        if name == "backbone_tree":
            if not isinstance(newick, str) or not newick.strip():
                problems.append(f"seq {ev.get('seq')}: backbone_tree has no newick")
                continue
            try:
                backbone = common.parse_newick(newick)
            except common.PrereqError as exc:
                problems.append(f"seq {ev.get('seq')}: backbone Newick is unparsable ({exc})")
                continue
            remap_leaf_names(backbone, id_to_name)
        else:
            if cluster_id is None or not isinstance(newick, str) or not newick.strip():
                problems.append(f"seq {ev.get('seq')}: cluster_focus without cluster id or newick")
                continue
            cid = str(cluster_id)
            state = clusters.setdefault(cid, ClusterState(cid))
            if not state.member_names:
                state.member_names = [
                    str(m.get("name")) for m in data.get("members") or []
                    if isinstance(m, dict) and m.get("name")
                ]
            emit_cluster_tree(cid, newick, focus)
        if position not in picked:
            continue
        run_splits, present, notes = displayed_splits()
        exact_mask = 0
        clusters_exact = 0
        for cid, entry in cache.items():
            if entry.get("focus") == 2:
                exact_mask |= entry["leaf_mask"]
                clusters_exact += 1
        rf_between: Optional[int] = None
        rf_within: Optional[int] = None
        if defined:
            result = common.rf_from_split_sets(run_splits, truth_splits)
            if truth_cluster_masks:
                differing = run_splits ^ truth_splits
                rf_between = sum(1 for split in differing if between_clusters(split))
                rf_within = len(differing) - rf_between
        else:
            result = common.RfResult()
        rows.append(
            {
                "dataset": dataset,
                "config": config,
                "t_ms": common.event_t_ms(ev),
                "seq": ev.get("seq"),
                "event": name,
                "cluster": None if cluster_id is None else str(cluster_id),
                "focus": focus,
                "rf_source": FULL_LEAF_SET_SOURCE,
                "leaves_run": common.popcount(present & universe),
                "leaves_truth": len(truth.leaves),
                "leaves_shared": n_scored,
                "splits_run": result.splits_a,
                "splits_truth": result.splits_b,
                "splits_shared": result.shared_splits,
                "leaves_exact": common.popcount(exact_mask & universe),
                "clusters_exact": clusters_exact,
                "rf": result.rf,
                "rf_normalized": result.rf_normalized,
                "rf_between_clusters": rf_between,
                "rf_within_clusters": rf_within,
                "collapse_threshold": truth.threshold or None,
                "note": "; ".join(notes),
            }
        )
    return rows, problems


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
    def cell(r: Dict[str, Any], column: str) -> Any:
        if column == "t_s":
            return None if r.get("t_ms") is None else round(float(r["t_ms"]) / 1000.0, 6)
        return r.get(column)

    common.write_csv(csv_path, CSV_HEADER, [[cell(r, c) for c in CSV_HEADER] for r in rows])
    return csv_path


def _curve_points(rows: Sequence[Dict[str, Any]]) -> List[Tuple[float, float]]:
    """``(seconds, normalised RF)`` in time order; undefined RF is left out."""
    pts = []
    for r in rows:
        if r.get("t_ms") is None or r.get("rf_normalized") in (None, ""):
            continue
        try:
            pts.append((max(float(r["t_ms"]) / 1000.0, 1e-3), float(r["rf_normalized"])))
        except (TypeError, ValueError):
            continue
    return sorted(pts, key=lambda p: p[0])


def _ordered_configs(names: Iterable[str]) -> List[str]:
    names = set(names)
    known = [c for c in CONFIG_STYLES if c in names]
    return known + sorted(names - set(known))


def plot_panels(
    panels: Sequence[Dict[str, Any]],
    png_path: Path,
    *,
    title: str,
    log_x: Optional[bool] = None,
    footnote: Optional[str] = ACCURACY_DEFINITION,
) -> Tuple[bool, Optional[str]]:
    """Small multiples: one panel per dataset, one curve per configuration.

    ``panels`` is a list of ``{"dataset", "label", "series": {config: rows},
    "empty_reason"}``. A panel without data is still drawn, with its reason, so
    a dataset that could not be scored is visible rather than silently absent.
    Each configuration keeps its fixed colour, marker and dash pattern
    (:data:`CONFIG_STYLES`) in every panel, so one legend serves the figure.
    The curve is drawn as steps: the tree on screen stays the same until the
    next event replaces it. ``log_x=None`` picks a log time axis per panel when
    its times span more than a factor of five.
    """
    plt, reason = common.try_import_matplotlib()
    if plt is None:
        return False, reason
    panels = [p for p in panels if p]
    if not any(_curve_points(rows) for p in panels for rows in (p.get("series") or {}).values()):
        return False, "no data points to plot"
    try:
        import textwrap

        from matplotlib.lines import Line2D

        n = len(panels)
        ncols = 1 if n == 1 else (2 if n in (2, 4) else 3)
        nrows = (n + ncols - 1) // ncols
        panel_w, panel_h = (7.2, 3.6) if n == 1 else (4.4, 3.0)
        all_configs = _ordered_configs(c for p in panels for c in (p.get("series") or {}))
        width = panel_w * ncols + 0.9
        wrap_chars = int(width * 15)
        foot = textwrap.wrap(footnote, wrap_chars) if footnote else []
        # about 2.2 in per legend entry (dash sample + the longest config name)
        legend_cols = max(1, min(len(all_configs), int((width - 0.9) // 2.2)))
        legend_rows = (len(all_configs) + legend_cols - 1) // legend_cols if len(all_configs) > 1 else 0
        top_in = 0.6 + 0.26 * legend_rows + (0.2 if legend_rows else 0.0)
        bottom_in = 0.62 + 0.15 * len(foot)
        height = panel_h * nrows + top_in + bottom_in + 0.35 * (nrows - 1)
        fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), dpi=150, squeeze=False)
        fig.patch.set_facecolor(SURFACE)
        fig.subplots_adjust(
            left=0.75 / width,
            right=1 - 0.2 / width,
            top=1 - top_in / height,
            bottom=bottom_in / height,
            hspace=0.6 * nrows / max(1.0, panel_h * nrows),
            wspace=0.22,
        )
        for slot, ax in enumerate(axes.flat):
            if slot >= n:
                ax.set_visible(False)
                continue
            panel = panels[slot]
            ax.set_facecolor(SURFACE)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            for side in ("left", "bottom"):
                ax.spines[side].set_color("#c3c2b7")
            ax.tick_params(colors=INK_MUTED, labelsize=8)
            ax.set_title(panel.get("label") or panel.get("dataset") or "", loc="left",
                         fontsize=9.5, color=INK)
            series = panel.get("series") or {}
            drawn = 0
            xs_all: List[float] = []
            for config in _ordered_configs(series):
                pts = _curve_points(series[config])
                if not pts:
                    continue
                colour, marker, dash = config_style(config, all_configs)
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                xs_all.extend(xs)
                ax.plot(
                    xs, ys, drawstyle="steps-post", color=colour, linestyle=dash,
                    linewidth=1.8, marker=marker, markersize=5.5,
                    markevery=max(1, len(xs) // 10), markerfacecolor=colour,
                    markeredgecolor=SURFACE, markeredgewidth=0.9, zorder=3,
                )
                drawn += 1
            if not drawn:
                ax.set_xticks([])
                ax.set_yticks([])
                why = panel.get("empty_reason") or "no accuracy points"
                ax.text(0.5, 0.5, "\n".join(textwrap.wrap(why, 44)), ha="center",
                        va="center", fontsize=8, color=INK_MUTED, transform=ax.transAxes)
                continue
            use_log = log_x
            if use_log is None:
                use_log = bool(xs_all) and max(xs_all) / max(min(xs_all), 1e-3) > 5.0
            ax.set_xscale("log" if use_log else "linear")
            if use_log:
                from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

                ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _pos: f"{v:g}"))
                ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_ylim(-0.04, 1.04)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.grid(True, which="major", color="#e4e3df", linewidth=0.6, zorder=0)
            row, col = divmod(slot, ncols)
            if col == 0:
                ax.set_ylabel("normalised RF to truth\n(1 = nothing right, 0 = exact)",
                              fontsize=8.5, color=INK)
            last_in_column = row == nrows - 1 or (slot + ncols) >= n
            if last_in_column:
                ax.set_xlabel(
                    "time since the run started (s" + (", log scale)" if use_log else ")"),
                    fontsize=8.5, color=INK,
                )
        fig.suptitle(title, x=0.75 / width, ha="left", y=1 - 0.12 / height,
                     fontsize=12, color=INK)
        if len(all_configs) > 1:
            handles = []
            for config in all_configs:
                colour, marker, dash = config_style(config, all_configs)
                handles.append(
                    Line2D([0], [0], color=colour, linestyle=dash, linewidth=1.8,
                           marker=marker, markersize=5.5, markerfacecolor=colour,
                           markeredgecolor=SURFACE, label=config)
                )
            fig.legend(
                handles=handles, loc="upper left", ncol=legend_cols, frameon=False,
                fontsize=8.5, bbox_to_anchor=(0.75 / width - 0.01, 1 - 0.48 / height),
                handlelength=3.2, columnspacing=1.4,
            )
        if foot:
            fig.text(0.75 / width, 0.12 / height, "\n".join(foot), fontsize=7.2,
                     color=INK_MUTED, ha="left", va="bottom")
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path, facecolor=SURFACE)
        plt.close(fig)
    except Exception as exc:  # noqa: BLE001 - plotting must never lose the numbers
        return False, f"plotting failed ({type(exc).__name__}: {exc})"
    return True, None


def dataset_label(
    name: str,
    meta: Optional[Dict[str, Any]] = None,
    mean_genome_bp: Optional[float] = None,
) -> str:
    """A panel title that says what the dataset is: size, and real or simulated."""
    meta = meta or {}
    bits: List[str] = []
    count = meta.get("genome_count")
    kind = meta.get("kind")
    if count:
        noun = f"{kind} genomes" if kind else "genomes"
        bits.append(f"{count} {noun}")
    elif kind:
        bits.append(str(kind))
    if mean_genome_bp:
        bits.append(f"~{mean_genome_bp / 1e6:.1f} Mb each")
    return f"{name} ({', '.join(bits)})" if bits else name


def labels_from_summary(summary: Optional[Dict[str, Any]]) -> Dict[str, str]:
    """Panel titles for every dataset a matrix summary knows."""
    lengths: Dict[str, float] = {}
    for cell in (summary or {}).get("cells") or []:
        run = cell.get("run") or {}
        total, genomes = run.get("total_length"), run.get("genomes")
        if cell.get("config") == "baseline" and total and genomes:
            lengths[str(cell.get("dataset"))] = float(total) / float(genomes)
    return {
        str(entry.get("name")): dataset_label(
            str(entry.get("name")), entry, lengths.get(str(entry.get("name")))
        )
        for entry in (summary or {}).get("datasets") or []
        if isinstance(entry, dict) and entry.get("name")
    }


def panels_from_series(
    series: Dict[str, List[Dict[str, Any]]],
    *,
    default_dataset: str = "dataset",
    labels: Optional[Dict[str, str]] = None,
    empty: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """Group ``{"dataset/config": rows}`` (or ``{config: rows}``) into panels."""
    grouped: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
    for key, rows in series.items():
        dataset, _, config = key.partition("/") if "/" in key else (default_dataset, "", key)
        grouped.setdefault(dataset, {})[config] = rows
    for dataset in (empty or {}):
        grouped.setdefault(dataset, {})
    return [
        {
            "dataset": dataset,
            "label": (labels or {}).get(dataset, dataset),
            "series": grouped[dataset],
            "empty_reason": (empty or {}).get(dataset),
        }
        for dataset in sorted(grouped)
    ]


def plot_curves(
    series: Dict[str, List[Dict[str, Any]]],
    png_path: Path,
    *,
    title: str,
    log_x: Optional[bool] = None,
    labels: Optional[Dict[str, str]] = None,
    empty: Optional[Dict[str, str]] = None,
    default_dataset: str = "dataset",
) -> Tuple[bool, Optional[str]]:
    """Figure for ``{"dataset/config": rows}``: one panel per dataset."""
    panels = panels_from_series(
        series, default_dataset=default_dataset, labels=labels, empty=empty
    )
    return plot_panels(panels, png_path, title=title, log_x=log_x)


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
    scoring: str = "full",
) -> Dict[str, Any]:
    """Produce the accuracy-over-time rows for one output folder.

    ``scoring="full"`` (the default, and what the key figure shows) scores the
    whole displayed tree over the full leaf set; ``"per-cluster"`` keeps the
    older measure - the mean RF of the cluster trees emitted so far, each
    against its own truth tree - for comparison with ``infinity bench``.
    """
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

    # --- collapse thresholds -------------------------------------------
    if comparable_bp is None:
        comparable_bp = common.comparable_bp_for(run_dir, run_json)
    result["collapse"] = {"threshold": None, "explanation": None, "comparable_bp": comparable_bp}

    # --- rows -----------------------------------------------------------
    # Which RF source gives the accuracy-over-time curve:
    #
    #   an explicit --rf-mode   that mode, and nothing else
    #   auto                    1. RF recomputed here at every emitted tree, over
    #                               the full leaf set, whenever a truth exists - the
    #                               one definition that is the same for simulated
    #                               and real data and has real time resolution
    #                            2. otherwise a bench.json series (labelled), because
    #                               one point beats none
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
    if rf_mode in ("auto", "internal"):
        truth, truth_problems = load_full_truth(
            truth_dir,
            truth_tree,
            collapse_snps=collapse_snps,
            collapse_mode=collapse_mode,
            collapse_length=collapse_length,
            comparable_bp=comparable_bp,
        )
        result["problems"].extend(truth_problems)
        if truth is not None:
            result["truth"] = {"kind": truth.kind, "source": truth.source, "leaves": len(truth.leaves)}
            result["collapse"]["threshold"] = truth.threshold
            result["collapse"]["explanation"] = f"truth: {truth.explanation}"
        if truth is not None and scoring == "per-cluster":
            per_cluster_truth: List[Tuple[str, common.NewickNode]] = []
            if truth.kind == "simulator" and truth_dir:
                per_cluster_truth, tc_problems = load_truth_trees(Path(truth_dir))
                result["problems"].extend(tc_problems)
            if not per_cluster_truth:
                result["problems"].append(
                    "--scoring per-cluster needs a simulator truth folder with tree_*.nwk"
                )
            else:
                points, problems = tree_points(log, max_points=max_points)
                result["problems"].extend(problems)
                internal_rows, problems = score_internal_per_cluster(
                    points, per_cluster_truth, truth.threshold, dataset, config
                )
                result["problems"].extend(problems)
                internal_source = "internal RF (per cluster)"
        elif truth is not None:
            thresholds = run_thresholds(
                collapse_snps=collapse_snps,
                collapse_mode=collapse_mode,
                collapse_length=collapse_length,
                comparable_bp=comparable_bp,
            )
            result["collapse"]["run"] = thresholds.explanation
            internal_rows, problems = score_full_leaf_set(
                log, truth, thresholds, dataset, config, max_points=max_points
            )
            result["problems"].extend(problems)
            internal_source = FULL_LEAF_SET_SOURCE

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
        if internal_rows:
            rows, source = internal_rows, internal_source
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
        help="where RF values come from (default: auto = RF recomputed here over the "
        "full leaf set whenever a truth exists, else a bench.json series)",
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
    p.add_argument(
        "--scoring",
        choices=("full", "per-cluster"),
        default="full",
        help="full (default): RF of the whole displayed tree over every genome, "
        "unrefined clusters as stars; per-cluster: mean RF of the cluster trees "
        "emitted so far, each against its own truth tree (simulated data only)",
    )
    p.add_argument("--max-points", type=int, default=400, help="cap scored points per run (default 400)")
    p.add_argument("--no-plot", action="store_true", help="write CSV only")
    p.add_argument(
        "--linear-x",
        action="store_true",
        help="linear time axis everywhere (default: log per panel when its times span > 5x)",
    )
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

    # A matrix summary remembers where every dataset's truth is, including
    # simulated datasets whose dataset.json lives outside the work folder.
    summary: Optional[Dict[str, Any]] = None
    summary_path = (
        Path(args.matrix_dir) / "summary.json" if args.matrix_dir
        else Path(args.run_dir).parent.parent / "summary.json"
    )
    obj, _ = common.load_json_optional(summary_path)
    if isinstance(obj, dict):
        summary = obj

    results: List[Dict[str, Any]] = []
    all_series: Dict[str, List[Dict[str, Any]]] = {}
    all_rows: List[Dict[str, Any]] = []
    empty_reasons: Dict[str, str] = {}
    for dataset, config, run_dir in targets:
        dataset_dir = work / dataset
        truth_tree = args.truth_tree
        truth_dir = args.truth_dir
        truth_notes: List[str] = []
        # The dataset descriptor knows both kinds of truth: a simulator truth
        # folder (per-cluster trees + pairwise SNPs) and, for real data, a
        # reference tree.
        # An explicit --truth-tree / --truth-dir is the truth, full stop: a
        # discovered simulator folder must not silently override it.
        explicit = truth_tree is not None or truth_dir is not None
        known_tree, known_dir, origin = (
            (None, None, "") if explicit else common.dataset_truth(work, dataset, summary)
        )
        if known_dir is not None:
            truth_dir = known_dir
            truth_notes.append(f"truth folder from {origin}: {known_dir}")
        if known_tree is not None:
            truth_tree = known_tree
            truth_notes.append(f"reference tree from {origin}: {known_tree}")
        if not explicit and truth_tree is None and truth_dir is None:
            truth_tree, more = find_truth_tree(
                run_dir, None, dataset_dir if dataset_dir.is_dir() else None
            )
            truth_notes.extend(more)
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
                scoring=args.scoring,
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
            all_series[f"{dataset}/{config}"] = rows
            all_rows.extend(
                [{**r, "dataset": dataset, "config": config} for r in rows]
            )
        elif res.get("problems"):
            empty_reasons.setdefault(dataset, res["problems"][0])

    # combined CSV + figure. The CSV is written first and always: a missing
    # matplotlib must never cost the numbers.
    combined_csv = out_dir / "accuracy_over_time.csv"
    write_rows(combined_csv, all_rows)
    print(f"\ncombined CSV: {combined_csv} ({len(all_rows)} rows)")

    png_written: Dict[str, Any] = {}
    log_x = False if args.linear_x else None
    with_points = {key.split("/", 1)[0] for key in all_series}
    empty = {d: why for d, why in empty_reasons.items() if d not in with_points}
    labels = labels_from_summary(summary)
    if not args.no_plot and all_series:
        title = "Tree accuracy over time"
        if len(targets) == 1:
            title += f" - {targets[0][0]}/{targets[0][1]}"
        png = out_dir / "accuracy_over_time.png"
        written, reason = plot_curves(
            all_series, png, title=title, log_x=log_x, empty=empty, labels=labels
        )
        if written:
            print(f"figure: {png}")
            png_written["combined"] = str(png)
        else:
            print(f"NOTE: no figure written - {reason}; the CSV above holds every number")
            png_written["error"] = reason
        # one figure per dataset when several were processed
        by_dataset: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
        for label, rows in all_series.items():
            d, _, c = label.partition("/")
            by_dataset.setdefault(d, {})[c] = rows
        if len(by_dataset) > 1:
            for d, series in sorted(by_dataset.items()):
                png = out_dir / f"accuracy_over_time__{d}.png"
                written, reason = plot_curves(
                    series, png, title=f"Tree accuracy over time - {d}", log_x=log_x,
                    default_dataset=d, labels=labels,
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
