#!/usr/bin/env python3
"""Outside baselines for PROJECT.md section 11 (Linux only, NOT part of M-Tato).

Two references are compared against M-Tato:

**Mash + neighbour-joining** -- the fast, blurry reference. Mash reduces each
genome to a small sketch and forgets position, which is exactly the trade-off
M-Tato claims to improve on. The neighbour-joining tree is built by this script
(pure Python, in ``common.py``), not by an external phylogenetics package.

**SKA2** -- the exact, alignment-free reference. Split k-mers see every SNP that
lies far enough from another one.

Neither tool is installed on the machine this was written on, so **every mode
detects a missing tool and degrades to consuming a pre-computed result instead of
failing**. That is not a shortcut: a pre-computed SKA2 distance table for the real
32-genome cluster already exists and is the reference this project is judged
against::

    .../analysis/ncbi_pd/ska/cluster_reads_distance.tsv
    header: Sample1  Sample2  Distance  Mismatches (proportion)  Match count  Mismatch count
    496 data rows for 32 genomes

``real_datasets.py`` copies that table (and the IQ-TREE tree built from the SKA2
core alignment) into ``bench-work/<dataset>/truth/``, and this script reads them
from there by default.

Sub-commands::

    status                       what is installed, what pre-computed data exists
    ska2   --dataset NAME        parse/normalise SKA2 distances  (or run ska2)
    mash   --dataset NAME        sketch + NJ tree                (or read a mash table)
    compare --dataset NAME --run-dir DIR
                                 join SKA2 against M-Tato's cluster_*.dist.tsv
    all    --dataset NAME --run-dir DIR

Outputs land in ``bench-work/baselines/<dataset>/`` and are always accompanied by
a TSV/JSON with the numbers, whether or not a figure could be drawn.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

TOOL_MASH = "mash"
TOOL_SKA2 = "ska2"


# ---------------------------------------------------------------------------
# Tool availability
# ---------------------------------------------------------------------------


def tool_status(name: str) -> Dict[str, Any]:
    """Where a tool is, what version it reports, and why it may be unusable."""
    info: Dict[str, Any] = {"tool": name, "found": False, "path": None, "version": None}
    path = shutil.which(name)
    if not path:
        info["reason"] = (
            f"{name} is not installed on this machine (nothing named '{name}' on PATH)"
        )
        return info
    info["found"] = True
    info["path"] = path
    if not sys.platform.startswith("linux"):
        info["reason"] = (
            f"the outside baselines are Linux only (PROJECT.md section 11); "
            f"this platform is {platform.system()}"
        )
        return info
    res = common.run_child([path, "--version"], capture=True, timeout=30, poll_memory=False)
    out, err = common.captured(res)
    text = (out or err or "").strip()
    info["version"] = text.splitlines()[0] if text else None
    if res.returncode != 0:
        info["reason"] = f"`{path} --version` exited {res.returncode}: {text[:200]}"
    return info


def linux_only_or_precomputed(action: str) -> Optional[str]:
    """A note when the platform cannot run outside tools at all."""
    if sys.platform.startswith("linux"):
        return None
    return (
        f"{action} needs Linux (PROJECT.md section 11 restricts the outside "
        f"baselines to Linux); this is {platform.system()}, so only "
        "pre-computed results can be used"
    )


# ---------------------------------------------------------------------------
# Mash
# ---------------------------------------------------------------------------

MASH_COLUMN_ALIASES = {
    "a": ("id", "query", "id1", "sample1", "#id", "ref"),
    "b": ("id2", "subject", "id-2", "sample2", "query"),
    "distance": ("distance", "dist", "ani-distance", "mash-distance"),
    "pvalue": ("p-value", "pvalue", "p"),
    "matching": ("matching-hashes", "matching_hashes", "matching", "shared-hashes"),
}


def parse_mash_table(path: Path) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Parse a ``mash dist -t`` table (or any two-sample distance table).

    ``mash dist -t`` writes **no header** at all - the columns are positionally
    ``ID  ID2  distance  p-value  matching-hashes`` - while a hand-made table may
    well have one. So the header is first looked up by name and, failing that, the
    file is re-read as headerless and the columns are taken positionally. Which of
    the two happened is recorded in the returned problems.
    """
    problems: List[str] = []
    table = common.read_table(path)
    cols: Dict[str, Optional[str]] = {
        field: table.find_column(*aliases) for field, aliases in MASH_COLUMN_ALIASES.items()
    }
    if cols["a"] is None or cols["b"] is None or cols["distance"] is None:
        # Re-read without consuming the first row as a header: `mash dist -t`
        # output starts straight away with data, and sample IDs are not numeric,
        # so the auto-detection above mistakes the first row for a header.
        table = common.read_table(path, has_header=False)
        if len(table.header) < 3:
            raise common.PrereqError(
                f"{path}: cannot find the two sample columns of a mash distance "
                f"table; the widest row has only {len(table.header)} field(s). "
                "Expected `mash dist -t` output: ID ID2 distance p-value matching-hashes."
            )
        cols = {
            "a": table.header[0],
            "b": table.header[1],
            "distance": table.header[2],
            "pvalue": table.header[3] if len(table.header) > 3 else None,
            "matching": table.header[4] if len(table.header) > 4 else None,
        }
        problems.append(
            f"{path.name} has no recognisable header; columns were taken positionally "
            "as ID, ID2, distance[, p-value, matching-hashes] (the `mash dist -t` layout)"
        )
    problems.extend(table.problems)
    rows: List[Dict[str, Any]] = []
    for row in table.rows:
        a = (row.get(cols["a"]) or "").strip()
        b = (row.get(cols["b"]) or "").strip()
        if not a or not b or a == b:
            continue
        rows.append(
            {
                "a": a,
                "b": b,
                "distance": common._to_float(row.get(cols["distance"])) if cols["distance"] else None,
                "pvalue": common._to_float(row.get(cols["pvalue"])) if cols.get("pvalue") else None,
                "matching": common._to_int(row.get(cols["matching"])) if cols.get("matching") else None,
            }
        )
    if not rows:
        problems.append(f"no usable rows in {path}")
    return rows, problems


def run_mash(
    genomes: Sequence[Path], out_dir: Path, *, k: int, sketch_size: int, threads: int
) -> Tuple[List[Dict[str, Any]], List[str], Dict[str, Any]]:
    """Sketch the genomes and produce an all-vs-all mash distance table."""
    notes: List[str] = []
    meta: Dict[str, Any] = {"k": k, "sketch_size": sketch_size, "threads": threads}
    mash = tool_status(TOOL_MASH)
    if not mash["found"]:
        notes.append(mash.get("reason") or f"{TOOL_MASH} not found")
        return [], notes, {"tool": mash, "ran": False}
    platform_note = linux_only_or_precomputed("running mash")
    if platform_note:
        notes.append(platform_note)
        return [], notes, {"tool": mash, "ran": False}
    common.ensure_dir(out_dir)
    list_file = out_dir / "mash_input.txt"
    list_file.write_text("\n".join(str(g) for g in genomes) + "\n", encoding="utf-8")
    sketch = out_dir / "mash"
    cmds = [
        [
            mash["path"], "sketch", "-o", str(sketch), "-l", str(list_file),
            "-k", str(k), "-s", str(sketch_size),
        ],
        [mash["path"], "dist", "-t", f"{sketch}.msh", f"{sketch}.msh"],
    ]
    for cmd in cmds:
        res = common.run_child(
            cmd,
            stdout_path=out_dir / "mash.log",
            stderr_path=out_dir / "mash.err",
            timeout=3600,
            poll_memory=False,
        )
        if not res.ok:
            notes.append(
                f"`{' '.join(cmd)}` failed with code {res.returncode}: {res.stderr_tail(300)}"
            )
            return [], notes, {"tool": mash, "ran": True, "failed": True, "commands": cmds}
    dist_table = out_dir / "mash.dist.tab"
    # `mash dist -t` writes to stdout; re-run it into the table file.
    res = common.run_child(
        [mash["path"], "dist", "-t", f"{sketch}.msh", f"{sketch}.msh"],
        stdout_path=dist_table,
        stderr_path=out_dir / "mash.err",
        timeout=3600,
        poll_memory=False,
    )
    if not res.ok:
        notes.append(f"writing {dist_table} failed: {res.stderr_tail(300)}")
        return [], notes, {"tool": mash, "ran": True, "failed": True}
    rows, problems = parse_mash_table(dist_table)
    notes.extend(problems)
    meta.update({"ran": True, "commands": cmds, "table": str(dist_table), "rows": len(rows)})
    return rows, notes, {"tool": mash, **meta}


def nj_tree_from_rows(rows: Sequence[Dict[str, Any]]) -> Tuple[str, List[str]]:
    """Neighbour-joining tree from a symmetric distance table."""
    notes: List[str] = []
    table: Dict[Tuple[str, str], float] = {}
    names: List[str] = []
    for row in rows:
        a, b = row["a"], row["b"]
        d = row.get("distance")
        if d is None:
            continue
        table[common.pair_key(a, b)] = float(d)
        for n in (a, b):
            if n not in names:
                names.append(n)
    if len(names) < 3:
        return "", [f"only {len(names)} samples have distances; NJ needs at least 3"]
    missing = 0

    def distance(a: str, b: str) -> float:
        nonlocal missing
        if a == b:
            return 0.0
        value = table.get(common.pair_key(a, b))
        if value is None:
            missing += 1
            # No mash distance for this pair: use the largest observed one so the
            # pair is placed far apart rather than falsely close. Recorded below.
            return max(table.values()) if table else 1.0
        return value

    newick = common.neighbor_joining(sorted(names), distance)
    if missing:
        notes.append(
            f"{missing} pair(s) had no distance and were given the maximum observed "
            "distance; the NJ topology near those samples is not meaningful"
        )
    return newick, notes


# ---------------------------------------------------------------------------
# SKA2
# ---------------------------------------------------------------------------


def run_ska2(genomes: Sequence[Path], out_dir: Path, *, k: int, threads: int) -> Tuple[Optional[Path], List[str], Dict[str, Any]]:
    """Build a SKA2 split-k-mer file and write its distance table."""
    notes: List[str] = []
    ska = tool_status(TOOL_SKA2)
    if not ska["found"]:
        notes.append(ska.get("reason") or f"{TOOL_SKA2} not found")
        return None, notes, {"tool": ska, "ran": False}
    platform_note = linux_only_or_precomputed("running ska2")
    if platform_note:
        notes.append(platform_note)
        return None, notes, {"tool": ska, "ran": False}
    common.ensure_dir(out_dir)
    list_file = out_dir / "ska_input.txt"
    list_file.write_text(
        "\n".join(f"{g.stem}\t{g}" for g in genomes) + "\n", encoding="utf-8"
    )
    skf = out_dir / "baseline.skf"
    dist = out_dir / "ska2_distance.tsv"
    cmds = [
        [ska["path"], "build", "-k", str(k), "-t", str(threads), "-o", str(skf), str(list_file)],
        [ska["path"], "distance", "--threads", str(threads), "-o", str(dist), str(skf)],
    ]
    for cmd in cmds:
        res = common.run_child(
            cmd,
            stdout_path=out_dir / "ska2.log",
            stderr_path=out_dir / "ska2.err",
            timeout=7200,
            poll_memory=False,
        )
        if not res.ok:
            notes.append(
                f"`{' '.join(cmd)}` failed with code {res.returncode}: {res.stderr_tail(300)}"
            )
            return None, notes, {"tool": ska, "ran": True, "failed": True, "commands": cmds}
    if not dist.exists():
        notes.append(f"ska2 finished but wrote no {dist}")
        return None, notes, {"tool": ska, "ran": True, "failed": True}
    return dist, notes, {"tool": ska, "ran": True, "commands": cmds, "table": str(dist)}


def normalise_ska(ska: common.SkaDistances, out_dir: Path) -> Path:
    """Write the parsed SKA2 table out in a stable, documented layout."""
    common.ensure_dir(out_dir)
    path = out_dir / "ska2.dist.tsv"
    common.write_tsv(
        path,
        ["sample1", "sample2", "snps", "mismatch_proportion", "match_count", "mismatch_count"],
        [
            [p.sample1, p.sample2, p.distance, p.mismatch_proportion, p.match_count, p.mismatch_count]
            for p in sorted(ska.pairs, key=lambda x: x.key())
        ],
    )
    return path


# ---------------------------------------------------------------------------
# The comparison that matters
# ---------------------------------------------------------------------------


def join_distances(
    mtato_rows: Sequence[Dict[str, Any]], ska_index: Dict[Tuple[str, str], common.SkaPair]
) -> Dict[str, Any]:
    """Join M-Tato's pairwise distances with SKA2's, pair by pair."""
    joined: List[Dict[str, Any]] = []
    mtato_index = {common.pair_key(r["a"], r["b"]): r for r in mtato_rows}
    for key, row in mtato_index.items():
        ska = ska_index.get(key)
        snps = row.get("snps")
        per_site = row.get("distance")
        comparable = row.get("comparable_bp")
        if snps is None and per_site is not None and comparable:
            snps = per_site * comparable  # substitutions/site -> SNP count
        joined.append(
            {
                "genome_a": key[0],
                "genome_b": key[1],
                "cluster": row.get("cluster"),
                "mtato_snps": snps,
                "mtato_snps_derived": row.get("snps") is None and snps is not None,
                "mtato_distance": per_site,
                "mtato_comparable_bp": comparable,
                "mtato_focus": row.get("focus"),
                "ska2_snps": ska.distance if ska else None,
                "ska2_match_count": ska.match_count if ska else None,
                "in_both": ska is not None,
            }
        )
    only_ska = [
        {"genome_a": k[0], "genome_b": k[1], "ska2_snps": v.distance}
        for k, v in ska_index.items()
        if k not in mtato_index
    ]
    for entry in joined:
        a, b = entry["mtato_snps"], entry["ska2_snps"]
        entry["abs_diff"] = (
            abs(float(a) - float(b)) if a is not None and b is not None else None
        )
    both = [e for e in joined if e["in_both"] and e["abs_diff"] is not None]
    diffs = [e["abs_diff"] for e in both]
    return {
        "pairs": sorted(joined, key=lambda e: (e["genome_a"], e["genome_b"])),
        "pairs_only_in_ska2": sorted(only_ska, key=lambda e: (e["genome_a"], e["genome_b"])),
        "summary": {
            "mtato_pairs": len(mtato_index),
            "ska2_pairs": len(ska_index),
            "pairs_in_both": len(both),
            "pairs_only_in_mtato": len([e for e in joined if not e["in_both"]]),
            "pairs_only_in_ska2": len(only_ska),
            "exact_matches": sum(1 for d in diffs if d == 0),
            "within_1_snp": sum(1 for d in diffs if d <= 1),
            "within_2_snp": sum(1 for d in diffs if d <= 2),
            "mean_abs_diff_snps": common.mean(diffs),
            "median_abs_diff_snps": common.median(diffs),
            "max_abs_diff_snps": max(diffs) if diffs else None,
            "p95_abs_diff_snps": common.quantile(diffs, 0.95),
            "pearson": common.pearson(
                [e["mtato_snps"] for e in both], [e["ska2_snps"] for e in both]
            ),
            "spearman": common.spearman(
                [e["mtato_snps"] for e in both], [e["ska2_snps"] for e in both]
            ),
        },
    }


def compare_trees(run_tree: Path, reference_tree: Path, collapse: float) -> Dict[str, Any]:
    """Robinson-Foulds distance between M-Tato's tree and a reference tree."""
    out: Dict[str, Any] = {
        "run_tree": str(run_tree),
        "reference_tree": str(reference_tree),
        "collapse_threshold": collapse,
    }
    for label, path in (("run", run_tree), ("reference", reference_tree)):
        if not path.exists():
            out["status"] = "skipped"
            out["reason"] = f"{label} tree not found: {path}"
            return out
    try:
        a = common.read_newick_file(run_tree)
        b = common.read_newick_file(reference_tree)
    except common.PrereqError as exc:
        out["status"] = "failed"
        out["reason"] = str(exc)
        return out
    result = common.rf_distance(a, b, collapse=collapse)
    out.update(result.as_dict())
    out["status"] = "ok" if result.ok else "failed"
    if not result.ok and not out.get("reason"):
        out["reason"] = "; ".join(result.problems)
    return out


# ---------------------------------------------------------------------------
# Dataset plumbing
# ---------------------------------------------------------------------------


def dataset_truth_paths(dataset: Dict[str, Any]) -> Dict[str, Optional[str]]:
    truth = dataset.get("truth") or {}
    dist = truth.get("distances") or {}
    return {
        "ska2_dist": dist.get("path") if dist.get("format") == "ska2" else None,
        "tree": truth.get("tree"),
        "dir": truth.get("dir"),
        "expected_values": truth.get("expected_values"),
    }


def load_dataset_by_name(work: Path, name: str) -> Dict[str, Any]:
    datasets = common.find_datasets(work, [name])
    if not datasets:
        raise common.PrereqError(f"no dataset named {name!r} under {work}")
    return datasets[0]


# ---------------------------------------------------------------------------
# Sub-command implementations
# ---------------------------------------------------------------------------


def cmd_status(args: argparse.Namespace) -> int:
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    print("Outside baselines (PROJECT.md section 11) - availability")
    print(f"platform: {platform.system()} {platform.release()} ({platform.machine()})")
    print()
    for tool in (TOOL_MASH, TOOL_SKA2):
        info = tool_status(tool)
        mark = "FOUND" if info["found"] else "MISSING"
        print(f"{tool:<6} {mark:<8} {info['path'] or '-'}")
        if info.get("version"):
            print(f"       version: {info['version']}")
        if info.get("reason"):
            print(f"       note:    {info['reason']}")
    print()
    if not work.is_dir():
        print(f"no benchmark work directory at {work}; run scripts/real_datasets.py first")
        return 0
    print(f"pre-computed reference data under {work}:")
    found_any = False
    for path in sorted(work.glob("*/dataset.json")):
        obj, problem = common.load_json_optional(path)
        if problem or not isinstance(obj, dict):
            print(f"  {path.parent.name}: unreadable dataset.json ({problem})")
            continue
        truth = dataset_truth_paths(obj)
        items = []
        if truth["ska2_dist"]:
            exists = Path(truth["ska2_dist"]).exists()
            items.append(f"ska2 distances: {truth['ska2_dist']} [{'ok' if exists else 'MISSING'}]")
            found_any |= exists
        if truth["tree"]:
            exists = Path(truth["tree"]).exists()
            items.append(f"reference tree: {truth['tree']} [{'ok' if exists else 'MISSING'}]")
            found_any |= exists
        if truth["expected_values"]:
            exists = Path(truth["expected_values"]).exists()
            items.append(f"expected values: {truth['expected_values']} [{'ok' if exists else 'MISSING'}]")
        print(f"  {obj.get('name')}:")
        for item in items or ["  (no reference data declared)"]:
            print(f"      {item}")
    baseline_dir = work / "baselines"
    if baseline_dir.is_dir():
        for path in sorted(baseline_dir.glob("*/*.tsv")):
            print(f"      cached baseline output: {path}")
            found_any = True
    if not found_any:
        print("\nnothing usable found: with no mash, no ska2 and no pre-computed "
              "table, the outside baselines cannot be produced")
    return 0


def _out_dir(args: argparse.Namespace, dataset_name: str) -> Path:
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    return common.ensure_dir(work / "baselines" / dataset_name)


def cmd_ska2(args: argparse.Namespace) -> int:
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    dataset = load_dataset_by_name(work, args.dataset)
    truth = dataset_truth_paths(dataset)
    out_dir = _out_dir(args, args.dataset)
    source: Optional[Path] = Path(args.ska_dist) if args.ska_dist else None
    provenance = "command line --ska-dist"
    if source is None and truth["ska2_dist"]:
        source = Path(truth["ska2_dist"])
        provenance = f"dataset.json truth.distances ({dataset.get('name')})"
    ran = False
    if source is None or not source.exists():
        if args.no_run:
            raise common.PrereqError(
                "no pre-computed SKA2 distance table was found "
                f"(looked at {truth['ska2_dist']!r}) and --no-run was given"
            )
        genomes = common.dataset_asm_files(dataset)
        source, notes, meta = run_ska2(genomes, out_dir, k=args.k, threads=args.threads or os.cpu_count() or 1)
        ran = source is not None
        for note in notes:
            print(f"NOTE: {note}")
        if source is None:
            print(
                "SKA2 unavailable and no pre-computed table: nothing to do.\n"
                "  supply one with --ska-dist <path> (a 'ska distance' TSV)"
            )
            common.dump_json(out_dir / "ska2.json", {"status": "unavailable", "notes": notes, "tool": meta.get("tool")})
            return 0
    else:
        meta = {"ran": False, "table": str(source)}
    ska = common.parse_ska2_distance_tsv(source)
    normalised = normalise_ska(ska, out_dir)
    stats = ska.stats()
    print(f"SKA2 distances: {source}")
    print(f"  provenance:   {provenance}" + (" (ran just now)" if ran else " (pre-computed)"))
    print(f"  rows:         {stats['rows']} ({stats['rows_with_distance']} with a distance)")
    print(f"  samples:      {stats['sample_count']}")
    print(f"  distance min/median/mean/max: "
          f"{common.fmt(stats['min'])} / {common.fmt(stats['median'])} / "
          f"{common.fmt(stats['mean'])} / {common.fmt(stats['max'])}")
    print(f"  normalised:   {normalised}")
    for problem in ska.problems:
        print(f"  NOTE: {problem}")
    doc = {
        "status": "ok",
        "tool": "ska2",
        "ran_here": ran,
        "source": str(source),
        "provenance": provenance,
        "normalised": str(normalised),
        "columns_detected": ska.columns,
        "stats": stats,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if args.tree:
        newick, notes = nj_tree_from_rows(
            [{"a": p.sample1, "b": p.sample2, "distance": p.distance} for p in ska.pairs]
        )
        if newick:
            tree_path = out_dir / "ska2.nj.tree.nwk"
            tree_path.write_text(newick + "\n", encoding="utf-8")
            doc["nj_tree"] = str(tree_path)
            print(f"  NJ tree:      {tree_path}")
        for note in notes:
            print(f"  NOTE: {note}")
    common.dump_json(out_dir / "ska2.json", doc)
    return 0


def cmd_mash(args: argparse.Namespace) -> int:
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    dataset = load_dataset_by_name(work, args.dataset)
    out_dir = _out_dir(args, args.dataset)
    precomputed: Optional[Path] = Path(args.mash_dist) if args.mash_dist else None
    if precomputed is None:
        for candidate in (out_dir / "mash.dist.tab", out_dir / "mash.dist.tsv"):
            if candidate.exists():
                precomputed = candidate
                break
    rows: List[Dict[str, Any]] = []
    notes: List[str] = []
    meta: Dict[str, Any] = {}
    if precomputed and precomputed.exists():
        rows, problems = parse_mash_table(precomputed)
        notes.extend(problems)
        meta = {"ran": False, "table": str(precomputed), "rows": len(rows)}
        print(f"Mash distances (pre-computed): {precomputed}")
    elif args.no_run:
        raise common.PrereqError(
            "no pre-computed mash distance table was found and --no-run was given; "
            f"supply one with --mash-dist <path> (looked in {out_dir})"
        )
    else:
        genomes = common.dataset_asm_files(dataset)
        rows, run_notes, meta = run_mash(
            genomes, out_dir, k=args.k, sketch_size=args.sketch_size,
            threads=args.threads or os.cpu_count() or 1,
        )
        notes.extend(run_notes)
        for note in run_notes:
            print(f"NOTE: {note}")
    if not rows:
        print(
            "Mash unavailable and no pre-computed table: the fast blurry reference "
            "cannot be produced.\n"
            "  this is reported, not hidden; install mash or pass --mash-dist <path>"
        )
        common.dump_json(
            out_dir / "mash.json",
            {"status": "unavailable", "notes": notes, "tool": meta.get("tool")},
        )
        return 0
    table_path = out_dir / "mash.dist.tsv"
    common.write_tsv(
        table_path,
        ["sample1", "sample2", "distance", "pvalue", "matching_hashes"],
        [[r["a"], r["b"], r["distance"], r.get("pvalue"), r.get("matching")] for r in rows],
    )
    distances = [r["distance"] for r in rows if r.get("distance") is not None]
    newick, nj_notes = nj_tree_from_rows(rows)
    notes.extend(nj_notes)
    tree_path = None
    if newick:
        tree_path = out_dir / "mash.nj.tree.nwk"
        tree_path.write_text(newick + "\n", encoding="utf-8")
    print(f"  rows: {len(rows)}  distance min/median/max: "
          f"{common.fmt(min(distances) if distances else None)} / "
          f"{common.fmt(common.median(distances))} / "
          f"{common.fmt(max(distances) if distances else None)}")
    print(f"  table: {table_path}")
    if tree_path:
        print(f"  NJ tree (neighbour-joining built by this script): {tree_path}")
    for note in nj_notes:
        print(f"  NOTE: {note}")
    common.dump_json(
        out_dir / "mash.json",
        {
            "status": "ok",
            "tool": "mash",
            "ran_here": bool(meta.get("ran")),
            "meta": meta,
            "table": str(table_path),
            "nj_tree": str(tree_path) if tree_path else None,
            "rows": len(rows),
            "distance_min": min(distances) if distances else None,
            "distance_median": common.median(distances),
            "distance_max": max(distances) if distances else None,
            "notes": notes,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    )
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    dataset = load_dataset_by_name(work, args.dataset)
    truth = dataset_truth_paths(dataset)
    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        raise common.PrereqError(
            f"--run-dir is not a directory: {run_dir}\n"
            "       point it at an mtato output folder (the one holding "
            "cluster_<id>.dist.tsv), e.g. bench-work/matrix/<dataset>/<config>"
        )
    out_dir = _out_dir(args, args.dataset)

    ska_source = Path(args.ska_dist) if args.ska_dist else (
        Path(truth["ska2_dist"]) if truth["ska2_dist"] else None
    )
    if ska_source is None or not ska_source.exists():
        raise common.PrereqError(
            f"no SKA2 distance table to compare against (looked at {ska_source}).\n"
            "       pass --ska-dist <path>, or run scripts/real_datasets.py so the "
            "dataset folder carries its own truth/"
        )
    ska = common.parse_ska2_distance_tsv(ska_source)
    ska_index = ska.as_dict()

    mtato_rows, problems = common.load_cluster_distances(run_dir)
    for problem in problems:
        print(f"NOTE: {problem}")
    if not mtato_rows:
        raise common.PrereqError(
            f"no pairwise distances could be read from {run_dir}\n"
            "       expected cluster_<id>.dist.tsv files with two genome columns and "
            f"a SNP or distance column; problems were: {problems}"
        )

    joined = join_distances(mtato_rows, ska_index)
    summary = dict(joined["summary"])
    tsv_path = out_dir / f"ska2_vs_mtato__{run_dir.name}.tsv"
    common.write_tsv(
        tsv_path,
        [
            "genome_a", "genome_b", "cluster", "mtato_snps", "ska2_snps", "abs_diff",
            "mtato_distance", "mtato_comparable_bp", "mtato_focus", "in_both",
            "mtato_snps_derived_from_per_site",
        ],
        [
            [e["genome_a"], e["genome_b"], e["cluster"], e["mtato_snps"], e["ska2_snps"],
             e["abs_diff"], e["mtato_distance"], e["mtato_comparable_bp"], e["mtato_focus"],
             e["in_both"], e["mtato_snps_derived"]]
            for e in joined["pairs"]
        ],
    )
    if joined["pairs_only_in_ska2"]:
        common.write_tsv(
            out_dir / f"ska2_pairs_missing_from_mtato__{run_dir.name}.tsv",
            ["genome_a", "genome_b", "ska2_snps"],
            [[e["genome_a"], e["genome_b"], e["ska2_snps"]] for e in joined["pairs_only_in_ska2"]],
        )

    print(f"M-Tato distances: {run_dir} ({len(mtato_rows)} pairs from "
          f"{len({r['cluster'] for r in mtato_rows})} cluster file(s))")
    print(f"SKA2 distances:   {ska_source} ({len(ska_index)} pairs)")
    print(f"  pairs in both:        {summary['pairs_in_both']}")
    print(f"  only in M-Tato:       {summary['pairs_only_in_mtato']}")
    print(f"  only in SKA2:         {summary['pairs_only_in_ska2']}")
    if summary["pairs_in_both"]:
        print(f"  exact matches:        {summary['exact_matches']}")
        print(f"  within 1 SNP:         {summary['within_1_snp']}")
        print(f"  within 2 SNPs:        {summary['within_2_snp']}")
        print(f"  mean |diff|:          {common.fmt(summary['mean_abs_diff_snps'])} SNPs")
        print(f"  median |diff|:        {common.fmt(summary['median_abs_diff_snps'])} SNPs")
        print(f"  max |diff|:           {common.fmt(summary['max_abs_diff_snps'])} SNPs")
        print(f"  p95 |diff|:           {common.fmt(summary['p95_abs_diff_snps'])} SNPs")
        print(f"  Pearson r:            {common.fmt(summary['pearson'], '.5f')}")
        print(f"  Spearman rho:         {common.fmt(summary['spearman'], '.5f')}")
    print(f"  table: {tsv_path}")

    doc: Dict[str, Any] = {
        "status": "ok" if summary["pairs_in_both"] else "no_overlap",
        "dataset": args.dataset,
        "run_dir": str(run_dir),
        "run_config": run_dir.name,
        "ska2_source": str(ska_source),
        "summary": summary,
        "problems": problems + list(ska.problems),
        "tsv": str(tsv_path),
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if not summary["pairs_in_both"]:
        doc["reason"] = (
            "no pair of genome names is present in both tables; M-Tato uses the "
            "assembly file name as the genome name, SKA2 uses the sample name it "
            "was given - they must agree for this join to mean anything"
        )
        print(f"NOTE: {doc['reason']}")

    # Tree comparison against the dataset's reference tree, when there is one.
    if args.truth_tree or truth["tree"]:
        reference = Path(args.truth_tree) if args.truth_tree else Path(truth["tree"])
        run_tree = Path(args.run_tree) if args.run_tree else run_dir / "tree.nwk"
        comparable_bp = common.comparable_bp_for(run_dir)
        threshold, explanation = common.guess_collapse_threshold(
            None, collapse_snps=args.collapse_snps, comparable_bp=comparable_bp, mode="auto"
        )
        tree_cmp = compare_trees(run_tree, reference, threshold)
        tree_cmp["collapse_explanation"] = explanation
        doc["tree_vs_reference"] = tree_cmp
        print(f"  tree vs reference ({reference.name}): "
              f"{tree_cmp.get('status')} rf={common.fmt(tree_cmp.get('rf'))} "
              f"rf_norm={common.fmt(tree_cmp.get('rf_normalized'), '.4f')} "
              f"shared_leaves={common.fmt(tree_cmp.get('leaves_shared'))}")
        if tree_cmp.get("reason"):
            print(f"    NOTE: {tree_cmp['reason']}")
        for problem in tree_cmp.get("problems", []) or []:
            print(f"    NOTE: {problem}")

    json_path = out_dir / f"compare__{run_dir.name}.json"
    common.dump_json(json_path, doc)
    print(f"  json:  {json_path}")
    return 0 if summary["pairs_in_both"] else 1


def cmd_all(args: argparse.Namespace) -> int:
    codes: Dict[str, int] = {}
    for name, func in (("ska2", cmd_ska2), ("mash", cmd_mash)):
        print(f"--- {name} ---")
        sub = argparse.Namespace(**{**vars(args), "no_run": args.no_run})
        try:
            codes[name] = func(sub)
        except common.PrereqError as exc:
            print(f"NOTE: {name} skipped: {exc}")
            codes[name] = 0
        print()
    if args.run_dir:
        print("--- compare ---")
        codes["compare"] = cmd_compare(args)
    else:
        print("--- compare --- skipped: no --run-dir given")
    return max(codes.values()) if codes else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="baselines.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--work-dir", type=Path, default=None, help="bench-work directory")
    sub = p.add_subparsers(dest="command", required=True)

    common_args = argparse.ArgumentParser(add_help=False)
    common_args.add_argument("--k", type=int, default=31, help="k-mer size for the tool (default 31)")
    common_args.add_argument("--threads", type=int, default=None, help="threads for the tool")
    common_args.add_argument(
        "--no-run",
        action="store_true",
        help="never invoke the external tool; only consume pre-computed results",
    )
    common_args.add_argument("--work-dir", type=Path, default=None, help=argparse.SUPPRESS)
    common_args.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

    dataset_arg = argparse.ArgumentParser(add_help=False)
    dataset_arg.add_argument(
        "--dataset", required=True, help="dataset name under bench-work/ (see --list of real_datasets.py)"
    )

    s = sub.add_parser("status", parents=[common_args], help="tool and data availability")
    s.set_defaults(func=cmd_status)

    s = sub.add_parser(
        "ska2", parents=[common_args, dataset_arg], help="SKA2 distances (exact reference)"
    )
    s.add_argument("--ska-dist", type=Path, default=None, help="pre-computed 'ska distance' TSV")
    s.add_argument("--tree", action="store_true", help="also build an NJ tree from the distances")
    s.set_defaults(func=cmd_ska2)

    s = sub.add_parser(
        "mash", parents=[common_args, dataset_arg], help="Mash + NJ (fast blurry reference)"
    )
    s.add_argument("--mash-dist", type=Path, default=None, help="pre-computed mash distance table")
    s.add_argument("--sketch-size", type=int, default=1000, help="mash sketch size (default 1000)")
    s.set_defaults(func=cmd_mash)

    s = sub.add_parser(
        "compare", parents=[common_args, dataset_arg], help="join SKA2 with an M-Tato run"
    )
    s.add_argument("--run-dir", required=True, type=Path, help="mtato output folder")
    s.add_argument("--ska-dist", type=Path, default=None, help="pre-computed 'ska distance' TSV")
    s.add_argument("--truth-tree", type=Path, default=None, help="reference tree (default: dataset truth)")
    s.add_argument("--run-tree", type=Path, default=None, help="M-Tato tree (default: <run-dir>/tree.nwk)")
    s.add_argument("--collapse-snps", type=float, default=0.5, help="RF collapse threshold in SNPs")
    s.set_defaults(func=cmd_compare)

    s = sub.add_parser("all", parents=[common_args, dataset_arg], help="ska2 + mash + compare")
    s.add_argument("--ska-dist", type=Path, default=None)
    s.add_argument("--mash-dist", type=Path, default=None)
    s.add_argument("--run-dir", type=Path, default=None, help="mtato output folder for the comparison")
    s.add_argument("--truth-tree", type=Path, default=None)
    s.add_argument("--run-tree", type=Path, default=None)
    s.add_argument("--collapse-snps", type=float, default=0.5)
    s.add_argument("--tree", action="store_true")
    s.add_argument("--sketch-size", type=int, default=1000)
    s.set_defaults(func=cmd_all)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.work_dir is None:
        args.work_dir = common.default_work_dir()
    for attr in ("ska_dist", "mash_dist", "truth_tree", "run_tree", "run_dir"):
        if getattr(args, attr, None) is None and hasattr(args, attr):
            setattr(args, attr, None)
    if not hasattr(args, "no_run"):
        args.no_run = False
    if not hasattr(args, "k"):
        args.k = 31
    if not hasattr(args, "threads"):
        args.threads = None
    return args.func(args)


if __name__ == "__main__":
    common.guard(main)()
