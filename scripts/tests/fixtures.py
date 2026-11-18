#!/usr/bin/env python3
"""Synthetic fixtures for testing the M-Tato benchmark scripts.

The real ``mtato`` binary does not exist yet, so these fixtures build everything
the scripts consume -- a tiny dataset, a truth tree, ``run.json``,
``events.jsonl``, ``cluster_*.dist.tsv`` -- plus a **fake ``mtato`` executable**
that speaks the documented CLI contract. That makes it possible to test
``run_matrix.py``, ``accuracy_over_time.py``, ``baselines.py`` and
``make_report.py`` end to end, including their failure paths, before the Rust
side is finished.

Everything is written into a directory the caller supplies (the tests pass a
``tempfile.mkdtemp()``, i.e. somewhere under ``/tmp``), so nothing here touches
the repository or the real datasets.
"""

from __future__ import annotations

import json
import os
import stat
import sys
import textwrap
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

SCRIPTS_DIR = Path(__file__).resolve().parent.parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")


# ---------------------------------------------------------------------------
# Genomes, trees, distances
# ---------------------------------------------------------------------------


def genome_names(n: int, prefix: str = "genome") -> List[str]:
    return [f"{prefix}{i:03d}" for i in range(n)]


def write_fasta(path: Path, name: str, length: int = 4096, seed: int = 0) -> Path:
    """A small deterministic FASTA file. Sequence content is irrelevant here."""
    bases = "ACGT"
    path.parent.mkdir(parents=True, exist_ok=True)
    chunks = []
    state = (seed * 1103515245 + 12345) & 0x7FFFFFFF
    for i in range(length):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        chunks.append(bases[(state >> 16) & 3])
    seq = "".join(chunks)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i : i + 70] + "\n")
    return path


def write_simulator_truth(
    truth_dir: Path,
    names: Sequence[str],
    *,
    clusters: int = 2,
    distances: Optional[Sequence[Tuple[str, str, float]]] = None,
) -> Dict[str, Any]:
    """A truth folder in the layout ``mtato simulate`` writes and ``bench`` reads.

    ``mtato bench --truth`` requires ``tree_*.nwk`` plus ``pairwise_snps.tsv`` and
    ``clusters.tsv``. The per-cluster truth trees are built with the same split
    the fake engine uses (first half / second half), so a correct fake run scores
    RF 0 and a wrong one does not.
    """
    truth_dir.mkdir(parents=True, exist_ok=True)
    items = list(names)
    per = max(1, (len(items) + clusters - 1) // clusters)
    groups = [g for g in (items[i * per : (i + 1) * per] for i in range(clusters)) if g]
    cluster_of: Dict[str, int] = {}
    for ci, group in enumerate(groups):
        (truth_dir / f"tree_c{ci:02d}.nwk").write_text(
            balanced_tree(group) + "\n", encoding="utf-8"
        )
        for gname in group:
            cluster_of[gname] = ci
    with open(truth_dir / "clusters.tsv", "w", encoding="utf-8") as fh:
        fh.write("genome\tcluster_id\n")
        for gname in items:
            fh.write(f"{gname}\t{cluster_of[gname]}\n")
    dist = {(a, b) if a <= b else (b, a): d for a, b, d in (distances or [])}
    with open(truth_dir / "pairwise_snps.tsv", "w", encoding="utf-8") as fh:
        fh.write("genome_a\tgenome_b\tsnps\n")
        for i, a in enumerate(items):
            for b in items[i + 1 :]:
                key = (a, b) if a <= b else (b, a)
                snps = dist.get(key, float(abs(i - items.index(b))))
                fh.write(f"{a}\t{b}\t{int(snps)}\n")
    (truth_dir / "manifest.json").write_text(
        json.dumps({"format_version": 1, "source": "scripts/tests/fixtures.py"}),
        encoding="utf-8",
    )
    return {"clusters": len(groups), "groups": groups}


def write_dataset(
    root: Path,
    name: str,
    names: Sequence[str],
    *,
    truth_tree: Optional[str] = None,
    truth_distances: Optional[Sequence[Tuple[str, str, float]]] = None,
    kind: str = "simulated",
    expected_values: Optional[Dict[str, Any]] = None,
    simulator_truth: bool = True,
    truth_clusters: int = 2,
) -> Path:
    """A ``bench-work``-style dataset folder: ``asm/`` plus ``dataset.json``."""
    dest = root / name
    asm = dest / "asm"
    asm.mkdir(parents=True, exist_ok=True)
    genomes = []
    for i, gname in enumerate(names):
        path = write_fasta(asm / f"{gname}.fasta", gname, length=2048, seed=i + 1)
        genomes.append(
            {
                "name": gname,
                "file": path.name,
                "source": str(path),
                "bytes": path.stat().st_size,
                "sha256": "0" * 64,
                "copy_verified": True,
            }
        )
    truth: Dict[str, Any] = {"kind": "none", "dir": None, "tree": None, "distances": None}
    if truth_tree or truth_distances:
        truth_dir = dest / "truth"
        truth_dir.mkdir(parents=True, exist_ok=True)
        truth["kind"] = "simulated"
        truth["dir"] = str(truth_dir)
        if truth_tree:
            tree_path = truth_dir / "truth.nwk"
            tree_path.write_text(truth_tree + "\n", encoding="utf-8")
            truth["tree"] = str(tree_path)
        if truth_distances:
            dist_path = truth_dir / "ska2_distance.tsv"
            lines = ["Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count"]
            for a, b, d in truth_distances:
                lines.append(f"{a}\t{b}\t{d:.2f}\t{d / 2000:.5f}\t2000\t{int(d) * 10}")
            dist_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            truth["distances"] = {"format": "ska2", "path": str(dist_path), "rows": len(truth_distances)}
        if simulator_truth:
            # Also lay the folder out the way `mtato simulate` does, so the real
            # `mtato bench --truth` (and the fake one) can consume it.
            meta = write_simulator_truth(
                truth_dir,
                names,
                clusters=truth_clusters,
                distances=truth_distances,
            )
            truth["simulator_layout"] = meta
    if expected_values:
        truth["expected_values"] = str(expected_values)
    doc = {
        "schema_version": 1,
        "name": name,
        "kind": kind,
        "organism": "Syntheticus testis",
        "description": f"synthetic fixture dataset with {len(names)} genomes",
        "asm_dir": str(asm),
        "genome_count": len(genomes),
        "total_bytes": sum(g["bytes"] for g in genomes),
        "copy_verified": True,
        "genomes": genomes,
        "source": {"paths": [str(asm)], "exists": {str(asm): True}},
        "licence": "fixture, no licence",
        "truth": truth,
        "notes": ["synthetic fixture created by scripts/tests/fixtures.py"],
        "caveats": [],
        "problems": [],
    }
    (dest / "dataset.json").write_text(json.dumps(doc, indent=2), encoding="utf-8")
    return dest


def balanced_tree(names: Sequence[str]) -> str:
    """A deterministic balanced Newick tree with SNP-ish integer branch lengths."""
    items = list(names)
    if not items:
        return "();"
    if len(items) == 1:
        return f"{items[0]};"

    def build(group: List[str], depth: int) -> str:
        if len(group) == 1:
            return f"{group[0]}:{depth + 1}"
        mid = len(group) // 2
        return f"({build(group[:mid], depth + 1)},{build(group[mid:], depth + 1)}):{depth + 1}"

    return build(items, 0) + ";"


def pair_distances(names: Sequence[str], base: int = 5) -> List[Tuple[str, str, float]]:
    """Deterministic 'true' pairwise SNP counts for every pair."""
    out = []
    items = list(names)
    for i, a in enumerate(items):
        for b in items[i + 1 :]:
            out.append((a, b, float(base + abs(i - items.index(b)))))
    return out


def perturb_tree(newick: str, swaps: int = 1) -> str:
    """Return a tree that differs from ``newick`` by ``swaps`` leaf exchanges.

    Used to make a fake run tree that is *not* identical to truth, so the RF code
    has something to measure.
    """
    import re

    tokens = re.findall(r"[A-Za-z0-9_]+", newick)
    leaves = [t for t in tokens if not t.isdigit()]
    if len(leaves) < 4 or swaps < 1:
        return newick
    out = newick
    for i in range(swaps):
        a = leaves[i]
        b = leaves[-(i + 1)]
        if a == b:
            continue
        out = out.replace(f"{a}:", f"@@{i}@@:").replace(f"{b}:", f"{a}:")
        out = out.replace(f"@@{i}@@:", f"{b}:")
    return out


# ---------------------------------------------------------------------------
# run.json / events.jsonl / output folder
# ---------------------------------------------------------------------------


def make_run_json(
    *,
    names: Sequence[str],
    clusters: int = 2,
    params: Optional[Dict[str, Any]] = None,
    first_tree_ms: Optional[int] = 120,
    exact_tree_ms: Optional[int] = 900,
    total_ms: Optional[int] = 1000,
    peak_memory_bytes: Optional[int] = 64 * 1024 * 1024,
    cancelled: bool = False,
    counters: Optional[Dict[str, Any]] = None,
    version: str = "0.0.0-fake",
) -> Dict[str, Any]:
    default_params = {
        "k": 31,
        "densities": [16384, 2048, 256],
        "sketch_level": 2,
        "cluster_dist": 0.005,
        "max_span_ratio": 10.0,
        "max_span_extra": 50000,
        "dp_max_len": 5000,
        "min_identity": 0.70,
        "xdrop": 15,
        "threads": 0,
        "block_skip": True,
        "hierarchical": True,
        "piece_cache": True,
        "profile_grouping": True,
        "single_cluster": False,
        "memory_budget_bytes": 0,
    }
    default_params.update(params or {})
    default_counters = {
        "blocks_looked_up": 4000,
        "blocks_matched": [1200, 900, 1800],
        "pieces_compared": 120,
        "cache_hits": 30,
        "alignments_run": 90,
        "bases_aligned": 45000,
        "kmers_scanned": len(names) * 2048,
        "genomes_scanned": len(names),
        "end_extensions": len(names) * 2,
        "stray_points_dropped": 3,
        "breakpoints": 1,
        "blocks_total": 4200,
    }
    default_counters.update(counters or {})
    n = len(names)
    per_cluster = max(1, (n + clusters - 1) // clusters)
    focus_counts = [0, 0, n] if exact_tree_ms is not None else [n, 0, 0]
    return {
        "version": version,
        "params": default_params,
        "notes": [],
        "genomes": n,
        "clusters": clusters,
        "total_length": n * 2048,
        "focus_counts": focus_counts,
        "timings": {
            "scan_ms": 50,
            "cluster_ms": 20,
            "backbone_ms": 10,
            "refine_ms": (exact_tree_ms or 0) - 80 if exact_tree_ms else None,
            "output_ms": 15,
            "total_ms": total_ms,
            "first_tree_ms": first_tree_ms,
            "exact_tree_ms": exact_tree_ms,
        },
        "counters": default_counters,
        "peak_memory_bytes": peak_memory_bytes,
        "cancelled": cancelled,
        "_per_cluster": per_cluster,
    }


def make_events(
    *,
    names: Sequence[str],
    clusters: int = 2,
    tree: str,
    cluster_trees: Sequence[str],
    start_wall_ms: int = 1_700_000_000_000,
) -> List[Dict[str, Any]]:
    """An ``events.jsonl`` record list in the documented envelope shape.

    The ``backbone_tree`` event carries a tree whose leaves are *clusters*
    (represented by their pivot names), as PROJECT.md section 6.4 requires; the
    ``tree`` argument is the full tree that ends up in ``tree.nwk``. Replacing
    each pivot leaf with its ``cluster_focus`` subtree must reproduce it, which
    is what makes the fixture a faithful double of the real engine.
    """
    items = list(names)
    per = max(1, (len(items) + clusters - 1) // clusters)
    groups = [items[i * per : (i + 1) * per] for i in range(clusters)]
    groups = [g for g in groups if g]
    events: List[Dict[str, Any]] = []
    seq = 0
    t = 0.0

    def emit(name: str, data: Any, dt: float = 10.0) -> None:
        nonlocal seq, t
        t += dt
        events.append(
            {
                "seq": seq,
                "t_ms": round(t, 3),
                "wall_ms": start_wall_ms + int(t),
                "event": {"event": name, "data": data},
            }
        )
        seq += 1

    for i, gname in enumerate(items):
        emit(
            "genome_scanned",
            {
                "id": i,
                "name": gname,
                "length": 2048,
                "contigs": 1,
                "n50": 2048,
                "gc": 0.5,
                "sketch_size": 64,
            },
            dt=1.0,
        )
    emit(
        "clustered",
        {
            "threshold": 0.005,
            "clusters": [
                {
                    "id": ci,
                    "pivot": ci * per,
                    "pivot_name": group[0],
                    "members": list(range(ci * per, ci * per + len(group))),
                    "member_names": list(group),
                    "focus": 0,
                    "cluster_dist": 0.005,
                }
                for ci, group in enumerate(groups)
            ],
        },
        dt=5.0,
    )
    pivot_names = [g[0] for g in groups]
    # Leaves are clusters, represented by their pivot (section 6.4).
    emit("backbone_tree", {"newick": balanced_tree(pivot_names), "pivots": pivot_names}, dt=5.0)
    for ci, group in enumerate(groups):
        for focus in (0, 1, 2):
            emit(
                "cluster_focus",
                {
                    "cluster": ci,
                    "focus": focus,
                    "approximate": focus < 2,
                    "newick": cluster_trees[ci] if ci < len(cluster_trees) else balanced_tree(group),
                    "members": [
                        {
                            "id": ci * per + j,
                            "name": gname,
                            "diff_pieces": 4 - focus,
                            "snps": focus * 2,
                            "indels": 0,
                            "uncovered_bp": 0,
                            "novel_bp": 0,
                        }
                        for j, gname in enumerate(group)
                    ],
                },
                dt=20.0,
            )
    emit(
        "counters",
        {
            "blocks_looked_up": 4000,
            "blocks_matched": [1200, 900, 1800],
            "pieces_compared": 120,
            "cache_hits": 30,
            "alignments_run": 90,
            "bases_aligned": 45000,
            "kmers_scanned": len(items) * 2048,
            "genomes_scanned": len(items),
            "end_extensions": len(items) * 2,
            "stray_points_dropped": 3,
            "breakpoints": 1,
            "blocks_total": 4200,
        },
        dt=1.0,
    )
    emit("progress", {"stage": "output", "done": 1, "total": 1, "elapsed_ms": int(t)}, dt=1.0)
    return events


def write_output_folder(
    outdir: Path,
    *,
    names: Sequence[str],
    clusters: int = 2,
    tree: Optional[str] = None,
    cluster_trees: Optional[Sequence[str]] = None,
    distances: Optional[Sequence[Tuple[str, str, float]]] = None,
    run_json: Optional[Dict[str, Any]] = None,
    events: Optional[Sequence[Dict[str, Any]]] = None,
    write_run_json: bool = True,
    write_events: bool = True,
    comparable_bp: int = 2048,
    only_listed_pairs: bool = False,
) -> Path:
    """A complete-looking ``mtato run`` output folder."""
    outdir.mkdir(parents=True, exist_ok=True)
    items = list(names)
    tree = tree or balanced_tree(items)
    per = max(1, (len(items) + clusters - 1) // clusters)
    groups = [items[i * per : (i + 1) * per] for i in range(clusters)]
    groups = [g for g in groups if g]
    if cluster_trees is None:
        cluster_trees = [balanced_tree(g) for g in groups]

    (outdir / "tree.nwk").write_text(tree + "\n", encoding="utf-8")
    with open(outdir / "genomes.tsv", "w", encoding="utf-8") as fh:
        fh.write("id\tname\tpath\tlength\tcontigs\tn50\tgc\tcluster\tis_pivot\tfocus\t"
                 "snps_vs_pivot\tindels_vs_pivot\tuncovered_bp\tnovel_bp\n")
        for i, gname in enumerate(items):
            cid = min(i // per, len(groups) - 1)
            fh.write(
                f"{i}\t{gname}\t/tmp/{gname}.fasta\t2048\t1\t2048\t0.5\t{cid}\t"
                f"{'true' if i % per == 0 else 'false'}\t2\t0\t0\t0\t0\n"
            )
    with open(outdir / "clusters.tsv", "w", encoding="utf-8") as fh:
        fh.write("cluster\tpivot\tsize\tfocus\n")
        for ci, group in enumerate(groups):
            fh.write(f"{ci}\t{group[0]}\t{len(group)}\t2\n")

    dist_index = {(a, b) if a <= b else (b, a): d for a, b, d in (distances or [])}
    for ci, group in enumerate(groups):
        (outdir / f"cluster_{ci}.nwk").write_text(cluster_trees[ci] + "\n", encoding="utf-8")
        with open(outdir / f"cluster_{ci}.dist.tsv", "w", encoding="utf-8") as fh:
            fh.write("genome_a\tgenome_b\tsnps\tcomparable_bp\tdistance\tfocus\n")
            for i, a in enumerate(group):
                for b in group[i + 1 :]:
                    key = (a, b) if a <= b else (b, a)
                    if key in dist_index:
                        snps = dist_index[key]
                    elif only_listed_pairs:
                        # A pair M-Tato never reported must stay absent, so an
                        # acceptance check can see it as "not measured".
                        continue
                    else:
                        snps = float(abs(items.index(a) - items.index(b)))
                    fh.write(f"{a}\t{b}\t{int(snps)}\t{comparable_bp}\t{snps / comparable_bp:.9f}\t2\n")
        with open(outdir / f"cluster_{ci}.variants.tsv", "w", encoding="utf-8") as fh:
            fh.write("genome\tpivot_contig\tpivot_position\tref\talt\ttype\t"
                     "genome_contig\tgenome_position\tstrand\tflags\n")
            for j, gname in enumerate(group):
                if j == 0:
                    continue
                fh.write(f"{gname}\tcontig1\t{100 + j}\tA\tG\tSNP\tcontig1\t{100 + j}\t+\t\n")

    if write_run_json:
        doc = run_json if run_json is not None else make_run_json(names=items, clusters=len(groups))
        doc = {k: v for k, v in doc.items() if not k.startswith("_")}
        (outdir / "run.json").write_text(json.dumps(doc, indent=2), encoding="utf-8")
    if write_events:
        records = (
            list(events)
            if events is not None
            else make_events(names=items, clusters=len(groups), tree=tree, cluster_trees=cluster_trees)
        )
        with open(outdir / "events.jsonl", "w", encoding="utf-8") as fh:
            for record in records:
                fh.write(json.dumps(record) + "\n")
    return outdir


def write_ska2_table(path: Path, rows: Sequence[Tuple[str, str, float]], *,
                     match_count: int = 2048, header: Optional[str] = None) -> Path:
    """A SKA2 ``ska distance``-shaped TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        header
        or "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count"
    ]
    for a, b, d in rows:
        lines.append(f"{a}\t{b}\t{d:.2f}\t{d / match_count:.5f}\t{match_count}\t{int(round(d * 10))}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The fake mtato executable
# ---------------------------------------------------------------------------

FAKE_MTATO_TEMPLATE = r'''#!/usr/bin/env python3
"""Fake `mtato` used only by the benchmark-script tests. Not part of the product.

Speaks the CLI contract from PROJECT.md section 9 and writes the documented
output files. Behaviour is steered with environment variables so the failure
paths of the benchmark scripts can be exercised:

  FAKE_MTATO_FAIL=1          exit non-zero after writing a partial folder
  FAKE_MTATO_NO_RUN_JSON=1   exit 0 without writing run.json
  FAKE_MTATO_NO_TIMINGS=1    write run.json with null timings
  FAKE_MTATO_NO_PEAKMEM=1    write run.json with peak_memory_bytes = null
  FAKE_MTATO_IGNORE_SWITCH=1 accept the switch flags but echo default params
  FAKE_MTATO_BENCH_FAIL=1    `bench` exits non-zero and writes nothing
  FAKE_MTATO_DIST_JSON=<p>   take the pairwise SNP counts from a JSON list of
                             [genome_a, genome_b, snps] triples
  FAKE_MTATO_SLEEP=<float>   sleep this many seconds (to make wall time visible)
"""
import json
import os
import sys
import time
from pathlib import Path

SCRIPTS_DIR = Path(__r__{scripts_dir!r}__r__)
sys.path.insert(0, str(SCRIPTS_DIR))
sys.path.insert(0, str(SCRIPTS_DIR / "tests"))

import common
import fixtures

FLAG_TO_PARAM = {{
    "--no-block-skip": ("block_skip", False),
    "--flat": ("hierarchical", False),
    "--no-piece-cache": ("piece_cache", False),
    "--no-profile-grouping": ("profile_grouping", False),
    "--single-cluster": ("single_cluster", True),
}}


def env_flag(name):
    return os.environ.get(name, "") not in ("", "0", "false", "False")


def collect_inputs(argv):
    names = []
    for token in argv:
        p = Path(token)
        if p.is_dir():
            for ext in fixtures.FASTA_EXTENSIONS:
                names.extend(sorted(f.stem for f in p.glob("*" + ext)))
        elif p.is_file():
            names.append(p.stem)
    seen = []
    for n in names:
        if n not in seen:
            seen.append(n)
    return sorted(seen)


def cmd_run(argv):
    outdir = None
    inputs = []
    params = {{}}
    skip = 0
    i = 0
    while i < len(argv):
        token = argv[i]
        if token == "-o":
            i += 1
            outdir = Path(argv[i])
        elif token in ("--threads", "--k", "--densities", "--cluster-dist"):
            i += 1
        elif token in FLAG_TO_PARAM:
            key, value = FLAG_TO_PARAM[token]
            params[key] = value
        elif token.startswith("-"):
            pass
        else:
            inputs.append(token)
        i += 1
    if outdir is None:
        print("fake mtato run: -o <outdir> is required", file=sys.stderr)
        return 2
    time.sleep(float(os.environ.get("FAKE_MTATO_SLEEP", "0.02")))
    names = collect_inputs(inputs)
    if not names:
        print(f"fake mtato run: no genomes found in {{inputs}}", file=sys.stderr)
        return 3
    clusters = 1 if params.get("single_cluster") else (2 if len(names) > 3 else 1)
    if env_flag("FAKE_MTATO_IGNORE_SWITCH"):
        params = {{}}
    # FAKE_MTATO_DIST_JSON lets a test dictate the exact pairwise SNP counts the
    # fake reports, so acceptance checks against published values can be tested.
    dist_json = os.environ.get("FAKE_MTATO_DIST_JSON")
    if dist_json:
        raw = json.loads(Path(dist_json).read_text())
        distances = [(r[0], r[1], float(r[2])) for r in raw]
    else:
        distances = fixtures.pair_distances(names)
    tree = fixtures.balanced_tree(names)
    per = max(1, (len(names) + clusters - 1) // clusters)
    groups = [g for g in (names[i2 * per:(i2 + 1) * per] for i2 in range(clusters)) if g]
    cluster_trees = [fixtures.balanced_tree(g) for g in groups]
    if env_flag("FAKE_MTATO_FAIL"):
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "tree.nwk").write_text(tree + "\n")
        print("fake mtato: simulated crash", file=sys.stderr)
        return 17
    fixtures.write_output_folder(
        outdir,
        names=names,
        clusters=len(groups),
        tree=tree,
        cluster_trees=cluster_trees,
        distances=distances,
        only_listed_pairs=bool(dist_json),
        write_run_json=not env_flag("FAKE_MTATO_NO_RUN_JSON"),
    )
    run_json_path = outdir / "run.json"
    if run_json_path.exists():
        doc = json.loads(run_json_path.read_text())
        doc["params"].update(params)
        if env_flag("FAKE_MTATO_NO_TIMINGS"):
            for key in list(doc["timings"]):
                doc["timings"][key] = None
        if env_flag("FAKE_MTATO_NO_PEAKMEM"):
            doc["peak_memory_bytes"] = None
        run_json_path.write_text(json.dumps(doc, indent=2))
    return 0


def cmd_bench(argv):
    outdir = None
    truth_tree = None
    truth_dir = None
    truth_distances = None
    out_json = None
    dataset = None
    config = None
    collapse_snps = 0.5
    comparable_bp = None
    i = 0
    positional = []
    while i < len(argv):
        token = argv[i]
        if token == "--truth-tree":
            i += 1
            truth_tree = Path(argv[i])
        elif token == "--truth":
            i += 1
            truth_dir = Path(argv[i])
        elif token == "--truth-distances":
            i += 1
            truth_distances = Path(argv[i])
        elif token == "--json":
            i += 1
            out_json = Path(argv[i])
        elif token == "--dataset":
            i += 1
            dataset = argv[i]
        elif token == "--config":
            i += 1
            config = argv[i]
        elif token == "--collapse-snps":
            i += 1
            collapse_snps = float(argv[i])
        elif token == "--comparable-bp":
            i += 1
            comparable_bp = int(argv[i])
        elif token.startswith("--truth-distances-format"):
            pass
        elif token.startswith("-"):
            pass
        else:
            positional.append(token)
        i += 1
    if positional:
        outdir = Path(positional[0])
    if outdir is None or not outdir.is_dir():
        print(f"fake mtato bench: no such output folder: {{outdir}}", file=sys.stderr)
        return 2
    if env_flag("FAKE_MTATO_BENCH_FAIL"):
        print("fake mtato bench: simulated failure", file=sys.stderr)
        return 19
    if truth_tree is None and truth_dir is not None:
        for candidate in ("truth.nwk", "cluster.treefile", "tree.nwk"):
            if (truth_dir / candidate).exists():
                truth_tree = truth_dir / candidate
                break
    if truth_tree is None or not truth_tree.exists():
        print("fake mtato bench: no truth tree", file=sys.stderr)
        return 4
    run_tree_path = outdir / "tree.nwk"
    if not run_tree_path.exists():
        print("fake mtato bench: output folder has no tree.nwk", file=sys.stderr)
        return 5
    run_tree = common.read_newick_file(run_tree_path)
    truth = common.read_newick_file(truth_tree)
    if comparable_bp is None:
        comparable_bp = common.comparable_bp_for(outdir)
    threshold, _ = common.guess_collapse_threshold(
        truth, collapse_snps=collapse_snps, comparable_bp=comparable_bp, mode="auto"
    )
    result = common.rf_distance(run_tree, truth, collapse=threshold)
    doc = {{
        "schema_version": 1,
        "output_dir": str(outdir),
        "truth_tree": str(truth_tree),
        "dataset": dataset,
        "config": config,
        "collapse_threshold": threshold,
        "tree": result.as_dict(),
    }}
    # distance comparison against a SKA2-shaped reference table
    if truth_distances and truth_distances.exists():
        ska = common.parse_ska2_distance_tsv(truth_distances)
        ska_index = {{common.pair_key(p.sample1, p.sample2): p.distance for p in ska.pairs}}
        rows, _problems = common.load_cluster_distances(outdir)
        diffs = []
        exact = 0
        for row in rows:
            key = common.pair_key(row["a"], row["b"])
            if key in ska_index and ska_index[key] is not None and row.get("snps") is not None:
                diff = abs(float(row["snps"]) - float(ska_index[key]))
                diffs.append(diff)
                exact += int(diff == 0)
        if diffs:
            doc["distances"] = {{
                "reference": str(truth_distances),
                "pairs_compared": len(diffs),
                "pairs_exact": exact,
                "mean_abs_error_snps": sum(diffs) / len(diffs),
                "max_abs_error_snps": max(diffs),
                "correlation": common.pearson(
                    [float(r["snps"]) for r in rows if r.get("snps") is not None][: len(diffs)],
                    diffs[: len(diffs)],
                ),
            }}
    # per-emitted-tree RF timeline, replayed from events.jsonl
    events_path = outdir / "events.jsonl"
    if events_path.exists():
        try:
            log = common.read_event_log(events_path)
            points, _ = _aot().tree_points(log)
            timeline = []
            for point in points:
                rf = common.rf_distance(point["tree"], truth, collapse=threshold)
                timeline.append({{
                    "t_ms": point["t_ms"],
                    "seq": point["seq"],
                    "event": point["event"],
                    "cluster": point["cluster"],
                    "focus": point["focus"],
                    "rf": rf.rf,
                    "rf_normalized": rf.rf_normalized,
                    "leaves_shared": rf.leaves_shared,
                    "splits_run": rf.splits_a,
                    "splits_truth": rf.splits_b,
                }})
            doc["rf_timeline"] = timeline
        except Exception as exc:
            doc["rf_timeline_error"] = f"{{type(exc).__name__}}: {{exc}}"
    out_json = out_json or (outdir / "bench.json")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(doc, indent=2))
    return 0


def _aot():
    import accuracy_over_time
    return accuracy_over_time


def main():
    argv = sys.argv[1:]
    if not argv:
        print(__doc__)
        return 2
    if argv[0] in ("--version", "-V"):
        print("mtato 0.0.0-fake")
        return 0
    if argv[0] == "run":
        return cmd_run(argv[1:])
    if argv[0] == "bench":
        return cmd_bench(argv[1:])
    print(f"fake mtato: unknown subcommand {{argv[0]}}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
'''


def write_fake_mtato(dest: Path, scripts_dir: Path = SCRIPTS_DIR) -> Path:
    """Write an executable fake ``mtato`` and return its path."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    source = FAKE_MTATO_TEMPLATE.replace(
        "__r__{scripts_dir!r}__r__", repr(str(scripts_dir))
    ).replace("{{", "{").replace("}}", "}")
    dest.write_text(source, encoding="utf-8")
    dest.chmod(dest.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return dest


# ---------------------------------------------------------------------------
# A whole synthetic bench-work tree
# ---------------------------------------------------------------------------


def make_work_dir(
    root: Path,
    *,
    dataset_names: Sequence[str] = ("fixture_small", "fixture_pairs"),
    genomes_per_dataset: int = 8,
    with_truth: bool = True,
) -> Path:
    """``root`` laid out like ``bench-work/`` with two fixture datasets."""
    root.mkdir(parents=True, exist_ok=True)
    for i, name in enumerate(dataset_names):
        n = genomes_per_dataset + i
        names = genome_names(n, prefix=f"{name[:4]}")
        tree = balanced_tree(names)
        distances = pair_distances(names)
        write_dataset(
            root,
            name,
            names,
            truth_tree=tree if with_truth else None,
            truth_distances=distances if with_truth else None,
            kind="simulated",
        )
    return root
