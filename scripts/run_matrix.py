#!/usr/bin/env python3
"""Run the full configuration x dataset matrix of PROJECT.md section 11.

For every dataset and every one of the six configurations -- the baseline
(everything on) plus one run per experiment switch (``--no-block-skip``,
``--flat``, ``--no-piece-cache``, ``--no-profile-grouping``,
``--single-cluster``) -- this script

1. invokes ``mtato run`` into its own output folder,
2. invokes ``mtato bench`` against the dataset's truth,
3. records wall time, peak memory, counters, timings and the bench result,
4. appends the cell to ``<matrix-dir>/summary.json``.

Two rules shape the implementation:

* **Resumable.** A configuration whose output folder already contains
  ``run.json`` is not re-run; it is reported as ``cached`` and its numbers are
  still read into the summary. Same for ``bench.json``. Use ``--force`` /
  ``--force-bench`` to redo them.
* **Never silently swallow a failure.** A configuration that fails is written
  into the summary with ``status: "failed"``, its return code, the tail of its
  stderr and its command line, and it appears in the printed table. It is never
  omitted: PROJECT.md section 2 lists "a benchmark that only reports the cases
  where the idea wins" as an anti-goal.

Runs are **sequential on purpose**: wall time and peak memory of one
configuration must not be contaminated by another running next to it.

Usage::

    scripts/run_matrix.py --dry-run                      # show the command lines
    scripts/run_matrix.py --datasets ncbi_pd_cluster     # one dataset, all configs
    scripts/run_matrix.py                                # everything in bench-work/
    scripts/run_matrix.py --configs baseline flat --force
"""

from __future__ import annotations

import argparse
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

SUMMARY_VERSION = 1
SUMMARY_NAME = "summary.json"

#: How long one ``mtato run`` may take before it is killed. Section 10.4 targets
#: 5 minutes to exact on 1000 genomes, so the default is generous.
DEFAULT_TIMEOUT_S = 6 * 3600.0


# ---------------------------------------------------------------------------
# Command construction
# ---------------------------------------------------------------------------


def mtato_version(binary: Path) -> Tuple[Optional[str], Optional[str]]:
    """``mtato --version`` output, or ``(None, reason)``."""
    try:
        res = common.run_child([binary, "--version"], capture=True, timeout=60, poll_memory=False)
    except OSError as exc:
        return None, f"could not execute {binary}: {exc}"
    out, err = common.captured(res)
    text = (out or err or "").strip()
    if res.returncode != 0:
        return None, f"`{binary} --version` exited {res.returncode}: {text[:200]}"
    return (text.splitlines()[0] if text else None), None


def run_command(
    binary: Path,
    inputs: Sequence[Path],
    outdir: Path,
    config: common.Config,
    args: argparse.Namespace,
) -> List[str]:
    """The exact ``mtato run`` command line for one cell."""
    cmd: List[str] = [str(binary), "run", *[str(p) for p in inputs], "-o", str(outdir)]
    if args.threads:
        cmd += ["--threads", str(args.threads)]
    if args.k:
        cmd += ["--k", str(args.k)]
    if args.densities:
        cmd += ["--densities", args.densities]
    if args.cluster_dist is not None:
        cmd += ["--cluster-dist", str(args.cluster_dist)]
    cmd += list(config.flags)
    cmd += list(args.extra_run_arg or [])
    return cmd


def truth_arguments(dataset: Dict[str, Any]) -> Tuple[List[str], Optional[str]]:
    """``mtato bench`` truth arguments for a dataset, or ``([], reason)``.

    ``mtato bench --truth`` requires a **simulator truth folder** (``tree_*.nwk``,
    ``pairwise_snps.tsv``, ``clusters.tsv``). A real dataset's reference data -- a
    SKA2 distance table and an IQ-TREE tree -- is not that layout, and handing it
    over would produce a bench failure that looks like an M-Tato bug. So for real
    datasets this returns no arguments and a reason naming what does the
    comparison instead: ``baselines.py compare`` (SKA2 join + RF against the
    reference tree) and ``accuracy_over_time.py --rf-mode internal``.
    """
    truth = dataset.get("truth") or {}
    if not isinstance(truth, dict):
        return [], "dataset.json has no usable 'truth' object"
    tdir = truth.get("dir")
    if tdir and Path(tdir).is_dir():
        if common.is_simulator_truth_dir(Path(tdir)):
            return ["--truth", str(tdir)], None
        return [], (
            f"truth folder {tdir} is not a `mtato simulate` truth folder (it has no "
            "tree_*.nwk / pairwise_snps.tsv), which `mtato bench --truth` requires; "
            "for this real dataset use `scripts/baselines.py compare` (SKA2 join and "
            "RF against the reference tree) and "
            "`scripts/accuracy_over_time.py --rf-mode internal`"
        )
    have_reference = bool(truth.get("tree") or truth.get("distances"))
    if have_reference:
        return [], (
            "this dataset has reference data (a tree and/or a distance table) but no "
            "`mtato simulate` truth folder, so `mtato bench` cannot score it; "
            "use `scripts/baselines.py compare` and "
            "`scripts/accuracy_over_time.py --rf-mode internal` instead"
        )
    return [], "dataset declares no truth of any kind"


def bench_command(
    binary: Path,
    outdir: Path,
    dataset: Dict[str, Any],
    config: common.Config,
    truth_args: Sequence[str],
    args: argparse.Namespace,
) -> List[str]:
    """The exact ``mtato bench`` command line for one cell.

    Matches the implemented interface::

        mtato bench --output <dir> --truth <truth dir> [-j <json>]
                    [--collapse-below <float>] [--quiet]
    """
    cmd: List[str] = [str(binary), "bench", "--output", str(outdir)]
    cmd += list(truth_args)
    cmd += ["--json", str(outdir / "bench.json")]
    if args.collapse_snps is not None:
        cmd += ["--collapse-below", str(args.collapse_snps)]
    cmd += ["--quiet"]
    cmd += list(args.extra_bench_arg or [])
    return cmd


def dataset_inputs(dataset: Dict[str, Any]) -> List[Path]:
    """What to hand ``mtato run`` for this dataset."""
    explicit = dataset.get("inputs")
    if isinstance(explicit, list) and explicit:
        paths = [Path(p) for p in explicit]
        missing = [str(p) for p in paths if not p.exists()]
        if missing:
            raise common.PrereqError(
                f"dataset {dataset.get('name')!r}: declared inputs do not exist: {missing}"
            )
        return paths
    return [Path(dataset["asm_dir"])]


# ---------------------------------------------------------------------------
# Cell execution
# ---------------------------------------------------------------------------


def summarise_run(outdir: Path, config: common.Config) -> Dict[str, Any]:
    """Read ``run.json`` from a finished cell into the summary."""
    out: Dict[str, Any] = {
        "run_json_present": False,
        "problems": [],
        "switch_verified": None,
        "timings": {},
        "counters": {},
    }
    run_json_path = outdir / "run.json"
    if not run_json_path.exists():
        out["problems"].append(
            f"no run.json in {outdir}: the run did not finish writing its summary"
        )
        return out
    rj = common.load_run_json(run_json_path, required=False)
    if rj is None:
        out["problems"].append(f"run.json in {outdir} is unreadable or not a JSON object")
        return out
    out["run_json_present"] = True
    out["problems"].extend(rj.unique_problems())
    out["version"] = rj.version
    out["cancelled"] = rj.cancelled
    out["genomes"] = rj.genomes
    out["clusters"] = rj.clusters
    out["total_length"] = rj.total_length
    out["focus_counts"] = rj.focus_counts
    out["peak_memory_bytes"] = rj.peak_memory_bytes
    out["timings"] = {key: rj.timing_ms(key) for key in common.TIMING_KEYS}
    out["timings_s"] = {
        "first_tree_s": rj.first_tree_s,
        "exact_tree_s": rj.exact_tree_s,
        "total_s": rj.total_s,
    }
    counters: Dict[str, Any] = {}
    for key in common.COUNTER_KEYS:
        counters[key] = rj.counter(key)
    counters["block_identity"] = rj.block_identity()
    out["counters"] = counters
    out["notes"] = rj.raw.get("notes") if isinstance(rj.raw.get("notes"), list) else []
    out["params_echo"] = rj.raw.get("params") if isinstance(rj.raw.get("params"), dict) else None

    # Did the switch actually take effect? The verdict is meaningless if a flag
    # was accepted and then ignored, so this is checked rather than assumed.
    if config.param:
        seen = rj.param(config.param)
        if seen is None:
            out["switch_verified"] = None
            out["problems"].append(
                f"run.json params has no {config.param!r}, so the {config.name} switch "
                "could not be verified to have taken effect"
            )
        else:
            out["switch_verified"] = bool(seen) is config.param_expected
            if not out["switch_verified"]:
                out["problems"].append(
                    f"SWITCH NOT APPLIED: --{config.name.replace('_','-')} should set "
                    f"params.{config.param}={config.param_expected} but run.json says {seen!r}"
                )
    else:
        # baseline: every optimisation must be on
        expected_on = {
            "block_skip": True,
            "hierarchical": True,
            "piece_cache": True,
            "profile_grouping": True,
            "single_cluster": False,
        }
        seen = {k: rj.param(k) for k in expected_on}
        if any(v is None for v in seen.values()):
            out["switch_verified"] = None
        else:
            out["switch_verified"] = all(bool(seen[k]) is v for k, v in expected_on.items())
            if not out["switch_verified"]:
                out["problems"].append(
                    f"BASELINE IS NOT THE FULL IDEA: params say {seen}, expected {expected_on}"
                )
    if rj.cancelled:
        out["problems"].append("run.json reports cancelled=true: the run did not complete")
    return out


def memory_summary(run_info: Dict[str, Any], child: Optional[common.ChildResult]) -> Dict[str, Any]:
    """Peak memory, with its source and whether it is only a lower bound."""
    from_json = run_info.get("peak_memory_bytes")
    if from_json:
        return {
            "bytes": int(from_json),
            "source": "run.json peak_memory_bytes",
            "is_floor": False,
        }
    polled = child.peak_rss_bytes_polled if child else None
    if polled:
        return {
            "bytes": int(polled),
            "source": "VmHWM sampled from /proc while the run was in flight",
            "is_floor": True,
        }
    return {"bytes": None, "source": None, "is_floor": True}


def summarise_bench(path: Path) -> Dict[str, Any]:
    bench = common.load_bench_json(path, required=False)
    if bench is None:
        return {
            "status": "failed",
            "path": str(path),
            "reason": f"{path} is missing or unreadable",
            "problems": [],
        }
    return {
        "status": "ok",
        "path": str(path),
        "reason": None,
        # rf is the SUM of the per-cluster RF distances and rf_normalized their
        # MEAN: the engine scores each cluster tree against its own truth tree, so
        # there is no single full-tree RF to report.
        "rf": bench.rf,
        "rf_is_sum_over_clusters": True,
        "rf_normalized": bench.rf_normalized,
        "rf_normalized_is_mean_over_clusters": True,
        "trees_scored": bench.n_trees_scored,
        "distance_mae_snps": bench.distance_mae,
        "distance_max_error_snps": bench.distance_max_error,
        "pairs_in_truth": bench.pairs_in_truth,
        "pairs_compared": bench.pairs_compared,
        "pairs_exact": bench.pairs_exact,
        "pairs_missing_from_output": bench.pairs_missing_from_output,
        "fraction_pairs_exact": bench.fraction_exact,
        "worst_pair": bench.worst_pair,
        "cluster_f1": bench.cluster_f1,
        "clusters_recovered_exactly": bench.clusters_recovered,
        "clusters_true": bench.clusters_true,
        "problems": bench.problems + bench.bench_problems,
        "raw_keys": sorted(bench.raw.keys()),
    }


def materialise_simulated(
    dataset: Dict[str, Any], binary: Path, args: argparse.Namespace, matrix_dir: Path
) -> Optional[str]:
    """Create a simulated dataset's assemblies and truth, if it declares how.

    ``dataset.json`` may carry ``{"simulate": {"command": ["simulate", ...]}}``;
    the command is run once, with ``-o <dataset dir>`` appended when it does not
    already contain one. Returns a problem string, or ``None`` on success.
    """
    spec = dataset.get("simulate")
    if not isinstance(spec, dict):
        return None
    command = spec.get("command")
    if not isinstance(command, list) or not command:
        return f"dataset {dataset.get('name')!r}: simulate.command is not a list"
    dest = Path(dataset.get("_dir") or (common.default_work_dir() / str(dataset.get("name"))))
    cmd = [str(binary), *[str(c) for c in command]]
    if not any(str(c) in ("-o", "--out", "--output") for c in cmd):
        cmd += ["-o", str(dest)]
    if args.dry_run:
        print("  [simulate] " + " ".join(cmd))
        return None
    log = dest / "simulate.log"
    dest.mkdir(parents=True, exist_ok=True)
    res = common.run_child(
        cmd,
        stdout_path=log,
        stderr_path=dest / "simulate.err",
        timeout=args.timeout,
        poll_memory=False,
    )
    if not res.ok:
        return (
            f"`{' '.join(cmd)}` failed with code {res.returncode}; "
            f"see {log} and {dest / 'simulate.err'}: {res.stderr_tail(400)}"
        )
    return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_matrix.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--datasets", nargs="+", default=None, help="dataset names (default: all)")
    p.add_argument(
        "--configs",
        nargs="+",
        default=None,
        choices=[c.name for c in common.CONFIGS],
        help="configurations to run (default: all six; the baseline is always included)",
    )
    p.add_argument("--work-dir", type=Path, default=None, help="where dataset.json folders live")
    p.add_argument(
        "--matrix-dir",
        type=Path,
        default=None,
        help="where run outputs and summary.json go (default: <work-dir>/matrix)",
    )
    p.add_argument("--mtato", default=None, help="path to the mtato binary")
    p.add_argument("--threads", type=int, default=None, help="passed through as --threads")
    p.add_argument("--k", type=int, default=None, help="passed through as --k")
    p.add_argument("--densities", default=None, help="passed through as --densities")
    p.add_argument("--cluster-dist", type=float, default=None, help="passed through as --cluster-dist")
    p.add_argument(
        "--extra-run-arg",
        action="append",
        default=[],
        help="extra argument for every `mtato run` (repeatable)",
    )
    p.add_argument(
        "--extra-bench-arg",
        action="append",
        default=[],
        help="extra argument for every `mtato bench` (repeatable)",
    )
    p.add_argument(
        "--collapse-snps",
        type=float,
        default=0.5,
        help="branch length below which the truth tree is collapsed for RF "
        "(section 11: 'half a SNP'); passed to `mtato bench` (default 0.5)",
    )
    p.add_argument(
        "--comparable-bp",
        type=int,
        default=None,
        help="comparable site count used to convert SNPs to substitutions/site",
    )
    p.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT_S, help="per-run timeout, seconds")
    p.add_argument("--force", action="store_true", help="re-run configurations that have a run.json")
    p.add_argument("--force-bench", action="store_true", help="re-run bench even if bench.json exists")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="print the exact command lines and change nothing",
    )
    p.add_argument("--json", type=Path, default=None, help="also write the summary to this path")
    p.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    matrix_dir = Path(args.matrix_dir) if args.matrix_dir else work / "matrix"

    datasets = common.find_datasets(work, args.datasets)
    configs = common.select_configs(args.configs)

    binary: Optional[Path] = None
    version: Optional[str] = None
    version_problem: Optional[str] = None
    try:
        binary = common.resolve_binary(args.mtato, required=not args.dry_run)
    except common.PrereqError as exc:
        if not args.dry_run:
            raise
        print(f"NOTE: {exc}", file=sys.stderr)
    if binary is not None:
        version, version_problem = mtato_version(binary)
        if version_problem:
            print(f"NOTE: {version_problem}", file=sys.stderr)
    binary_display = str(binary) if binary else "<mtato: NOT FOUND>"

    print(f"mtato:      {binary_display}" + (f" ({version})" if version else ""))
    print(f"work dir:   {work}")
    print(f"matrix dir: {matrix_dir}")
    print(f"datasets:   {', '.join(d['name'] for d in datasets)}")
    print(f"configs:    {', '.join(c.name for c in configs)}")
    if args.dry_run:
        print("MODE:       dry run - nothing will be executed or written\n")
    else:
        print()

    if not args.dry_run:
        common.ensure_dir(matrix_dir)

    summary_path = matrix_dir / SUMMARY_NAME

    # Merge with an existing summary so that re-running one dataset or one
    # configuration never erases the others. Losing an earlier *failed* cell
    # would be exactly the anti-goal of PROJECT.md section 2, so preserved cells
    # keep their status, error and command line.
    planned = {(str(d.get("name")), c.name) for d in datasets for c in configs}
    preserved_cells: List[Dict[str, Any]] = []
    preserved_datasets: List[Dict[str, Any]] = []
    if not args.dry_run and summary_path.exists():
        old, problem = common.load_json_optional(summary_path)
        if isinstance(old, dict):
            for cell in old.get("cells") or []:
                key = (str(cell.get("dataset")), cell.get("config"))
                if key not in planned:
                    cell = dict(cell)
                    cell["preserved_from_previous_run"] = True
                    preserved_cells.append(cell)
            this_names = {str(d.get("name")) for d in datasets}
            for entry in old.get("datasets") or []:
                if str(entry.get("name")) not in this_names:
                    preserved_datasets.append(entry)
            if preserved_cells:
                print(
                    f"preserving {len(preserved_cells)} cell(s) from the previous summary "
                    "that this invocation does not cover\n"
                )
        elif problem:
            print(f"NOTE: ignoring unreadable {summary_path}: {problem}", file=sys.stderr)

    cells: List[Dict[str, Any]] = list(preserved_cells)
    summary: Dict[str, Any] = {
        "schema_version": SUMMARY_VERSION,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "dry_run": bool(args.dry_run),
        "mtato": {
            "path": str(binary) if binary else None,
            "version": version,
            "version_problem": version_problem,
        },
        "host": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "cpu_count": os.cpu_count(),
        },
        "work_dir": str(work),
        "matrix_dir": str(matrix_dir),
        "request": {
            "threads": args.threads,
            "k": args.k,
            "densities": args.densities,
            "cluster_dist": args.cluster_dist,
            "extra_run_arg": list(args.extra_run_arg),
            "extra_bench_arg": list(args.extra_bench_arg),
            "collapse_snps": args.collapse_snps,
            "comparable_bp": args.comparable_bp,
            "timeout_s": args.timeout,
            "force": bool(args.force),
            "force_bench": bool(args.force_bench),
        },
        "configs": [
            {"name": c.name, "flags": list(c.flags), "removes": c.removes} for c in configs
        ],
        "datasets": [
            {
                "name": d.get("name"),
                "kind": d.get("kind"),
                "organism": d.get("organism"),
                "genome_count": d.get("genome_count"),
                "total_bytes": d.get("total_bytes"),
                "asm_dir": d.get("asm_dir"),
                "truth": d.get("truth"),
                "caveats": d.get("caveats", []),
                "problems": d.get("problems", []),
            }
            for d in datasets
        ]
        + preserved_datasets,
        "cells": cells,
    }

    failures: List[str] = []
    started_all = time.monotonic()

    for dataset in datasets:
        dname = str(dataset.get("name"))
        sim_problem = materialise_simulated(dataset, binary or Path("mtato"), args, matrix_dir)
        if sim_problem:
            failures.append(f"{dname}: {sim_problem}")
            print(f"[{dname}] simulation step FAILED: {sim_problem}", file=sys.stderr)
        try:
            inputs = dataset_inputs(dataset)
        except common.PrereqError as exc:
            failures.append(f"{dname}: {exc}")
            cells.append(
                {
                    "dataset": dname,
                    "config": None,
                    "status": "failed",
                    "stage": "inputs",
                    "error": str(exc),
                }
            )
            _flush(summary, summary_path, args, cells, failures, started_all)
            continue

        bargs, breason = truth_arguments(dataset)
        print(f"=== {dname} ({dataset.get('genome_count')} genomes) ===")
        if breason:
            print(f"    bench: SKIPPED - {breason}")
        else:
            print(f"    bench truth args: {' '.join(bargs)}")

        for config in configs:
            outdir = matrix_dir / dname / config.name
            cell: Dict[str, Any] = {
                "dataset": dname,
                "config": config.name,
                "flags": list(config.flags),
                "out_dir": str(outdir),
                "status": "pending",
                "reused": False,
                "problems": [],
            }
            cells.append(cell)
            run_cmd = run_command(binary or Path("mtato"), inputs, outdir, config, args)
            cell["run_command"] = run_cmd
            bench_cmd = (
                bench_command(binary or Path("mtato"), outdir, dataset, config, bargs, args)
                if bargs
                else None
            )
            cell["bench_command"] = bench_cmd

            have_run_json = (outdir / "run.json").exists()
            have_bench_json = (outdir / "bench.json").exists()

            if args.dry_run:
                cell["status"] = "dry-run"
                cell["reused"] = have_run_json and not args.force
                tag = "cached" if cell["reused"] else "would run"
                print(f"  [{config.name}] {tag}")
                print(f"      {' '.join(run_cmd)}")
                if bench_cmd:
                    btag = "cached" if (have_bench_json and not args.force_bench) else "would run"
                    print(f"      {btag}: {' '.join(bench_cmd)}")
                else:
                    print(f"      bench skipped: {breason}")
                continue

            common.ensure_dir(outdir)
            child: Optional[common.ChildResult] = None

            # ---- mtato run -------------------------------------------------
            if have_run_json and not args.force:
                cell["status"] = "cached"
                cell["reused"] = True
                print(f"  [{config.name}] reusing existing run.json in {outdir}")
            else:
                if outdir.exists() and args.force:
                    for stale in ("run.json", "bench.json", "events.jsonl", "tree.nwk"):
                        p = outdir / stale
                        if p.exists():
                            p.unlink()
                print(f"  [{config.name}] running: {' '.join(run_cmd)}")
                child = common.run_child(
                    run_cmd,
                    stdout_path=outdir / "mtato.run.log",
                    stderr_path=outdir / "mtato.run.err",
                    timeout=args.timeout,
                )
                cell["returncode"] = child.returncode
                cell["wall_s"] = round(child.wall_s, 3)
                cell["timed_out"] = child.timed_out
                if child.timed_out:
                    cell["status"] = "failed"
                    cell["error"] = f"timed out after {args.timeout:.0f}s and was killed"
                elif child.returncode != 0:
                    cell["status"] = "failed"
                    cell["error"] = f"`mtato run` exited with code {child.returncode}"
                    cell["stderr_tail"] = child.stderr_tail()
                else:
                    cell["status"] = "ok"
                print(
                    f"      -> {cell['status']} in {common.human_seconds(cell.get('wall_s'))}"
                    + (f": {cell.get('error')}" if cell.get("error") else "")
                )
                if cell["status"] == "failed" and not (outdir / "run.json").exists():
                    failures.append(f"{dname}/{config.name}: {cell.get('error')}")
                    cell["memory"] = memory_summary({}, child)
                    _flush(summary, summary_path, args, cells, failures, started_all)
                    continue

            # ---- read run.json --------------------------------------------
            run_info = summarise_run(outdir, config)
            cell["run"] = run_info
            cell["problems"].extend(run_info.get("problems", []))
            if not run_info.get("run_json_present"):
                if cell["status"] == "ok":
                    cell["status"] = "failed"
                    cell["error"] = "`mtato run` exited 0 but wrote no run.json"
                failures.append(f"{dname}/{config.name}: no usable run.json")
            cell["memory"] = memory_summary(run_info, child)
            if cell["memory"]["source"] is None:
                cell["problems"].append(
                    "peak memory unknown: run.json has no peak_memory_bytes and /proc "
                    "sampling was unavailable"
                )
            timings_s = run_info.get("timings_s") or {}
            cell["first_tree_s"] = timings_s.get("first_tree_s")
            cell["exact_tree_s"] = timings_s.get("exact_tree_s")
            cell["total_s"] = timings_s.get("total_s")
            cell["peak_memory_bytes"] = cell["memory"]["bytes"]
            cell["peak_memory_is_floor"] = bool(
                cell["memory"]["is_floor"] and cell["memory"]["bytes"] is not None
            )

            # ---- mtato bench ----------------------------------------------
            if bench_cmd is None:
                cell["bench"] = {"status": "skipped", "reason": breason, "path": None}
            elif (outdir / "bench.json").exists() and not args.force_bench:
                cell["bench"] = summarise_bench(outdir / "bench.json")
                cell["bench"]["reused"] = True
                print(f"      bench: reusing {outdir / 'bench.json'}")
            else:
                print(f"      bench: {' '.join(bench_cmd)}")
                bchild = common.run_child(
                    bench_cmd,
                    stdout_path=outdir / "mtato.bench.log",
                    stderr_path=outdir / "mtato.bench.err",
                    timeout=args.timeout,
                    poll_memory=False,
                )
                cell["bench_command_returncode"] = bchild.returncode
                cell["bench_wall_s"] = round(bchild.wall_s, 3)
                if (outdir / "bench.json").exists():
                    cell["bench"] = summarise_bench(outdir / "bench.json")
                    if bchild.returncode != 0:
                        cell["bench"]["problems"].append(
                            f"`mtato bench` exited {bchild.returncode} but wrote bench.json"
                        )
                else:
                    cell["bench"] = {
                        "status": "failed",
                        "path": str(outdir / "bench.json"),
                        "reason": (
                            f"`mtato bench` exited {bchild.returncode} and wrote no bench.json: "
                            + (bchild.stderr_tail(400) or "no stderr")
                        ),
                        "problems": [],
                    }
                    failures.append(f"{dname}/{config.name}: bench failed")
                print(f"      bench: {cell['bench']['status']}")

            _flush(summary, summary_path, args, cells, failures, started_all)
        print()

    _flush(summary, summary_path, args, cells, failures, started_all, final=True)
    if args.json:
        common.dump_json(args.json, summary)
        print(f"summary also written to {args.json}")

    _print_table(summary)
    reported = summary.get("cells", [])
    ok = sum(1 for c in reported if c.get("status") in ("ok", "cached"))
    failed = sum(1 for c in reported if c.get("status") == "failed")
    bench_failed = sum(
        1 for c in reported if (c.get("bench") or {}).get("status") == "failed"
    )
    bench_skipped = sum(
        1 for c in reported if (c.get("bench") or {}).get("status") == "skipped"
    )
    print(
        f"\ncells: {len(reported)} total | runs: {ok} usable, {failed} FAILED, "
        f"{len(reported) - ok - failed} other | "
        f"bench: {bench_failed} FAILED, {bench_skipped} skipped (no truth)"
    )
    if bench_failed:
        print(
            f"       a failed bench does not invalidate the run: timings, memory and "
            f"counters for those {bench_failed} cell(s) are still usable, but their "
            "accuracy is unknown and is reported as unknown"
        )
    if summary.get("preserved_cells"):
        print(
            f"({summary['preserved_cells']} of those cells were preserved from an "
            "earlier invocation covering other datasets/configurations)"
        )
    all_failures = summary.get("failures", [])
    if all_failures:
        print("\nFAILURES (every one of these is in summary.json, none is omitted):")
        for f in all_failures:
            print(f"  - {f}")
        return 1
    return 0


def _flush(
    summary: Dict[str, Any],
    path: Path,
    args: argparse.Namespace,
    cells: List[Dict[str, Any]],
    failures: List[str],
    started_all: float,
    final: bool = False,
) -> None:
    """Write the summary after every cell so an interrupted run is still usable.

    The failure list is *derived from every cell*, including the ones preserved
    from an earlier invocation, so a failure can never fall out of the summary by
    re-running a subset of the matrix.
    """
    if args.dry_run:
        return
    dataset_order = [str(d.get("name")) for d in summary.get("datasets", [])]
    config_order = [c.name for c in common.CONFIGS]

    def sort_key(cell: Dict[str, Any]) -> Tuple[int, int, str, str]:
        d = str(cell.get("dataset"))
        c = str(cell.get("config"))
        return (
            dataset_order.index(d) if d in dataset_order else len(dataset_order),
            config_order.index(c) if c in config_order else len(config_order),
            d,
            c,
        )

    ordered = sorted(cells, key=sort_key)
    derived: List[str] = []
    for cell in ordered:
        label = f"{cell.get('dataset')}/{cell.get('config') or '-'}"
        if cell.get("status") == "failed":
            derived.append(f"{label}: {cell.get('error') or 'failed'}")
        bench = cell.get("bench") or {}
        if bench.get("status") == "failed":
            derived.append(f"{label}: bench failed")
    all_failures = list(dict.fromkeys([f for f in derived if f] + [f for f in failures if f]))

    totals: Dict[str, int] = {}
    for cell in ordered:
        status = str(cell.get("status"))
        totals[status] = totals.get(status, 0) + 1
    summary["cells"] = ordered
    summary["totals"] = {"cells": len(ordered), **totals}
    summary["failures"] = all_failures
    summary["preserved_cells"] = sum(1 for c in ordered if c.get("preserved_from_previous_run"))
    summary["elapsed_s"] = round(time.monotonic() - started_all, 1)
    summary["complete"] = final
    summary["updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    common.dump_json(path, summary)


def _print_table(summary: Dict[str, Any]) -> None:
    """A one-line-per-cell table. Failed cells are shown, never dropped."""
    cells = summary.get("cells", [])
    if not cells:
        return
    print(
        f"{'dataset':<18} {'config':<21} {'status':<9} {'wall':>9} "
        f"{'1st tree':>9} {'exact':>9} {'peak mem':>10} {'RF':>6}  notes"
    )
    print("-" * 112)
    for cell in cells:
        run = cell.get("run") or {}
        mem = cell.get("memory") or {}
        bench = cell.get("bench") or {}
        notes: List[str] = []
        if cell.get("reused"):
            notes.append("cached")
        if run.get("switch_verified") is False:
            notes.append("SWITCH NOT APPLIED")
        if mem.get("is_floor") and mem.get("bytes"):
            notes.append("mem is a floor")
        if bench.get("status") == "failed":
            notes.append("bench failed")
        elif bench.get("status") == "skipped":
            notes.append("bench skipped (no truth)")
        if cell.get("error"):
            notes.append(str(cell["error"])[:48])
        print(
            f"{str(cell.get('dataset'))[:18]:<18} {str(cell.get('config'))[:21]:<21} "
            f"{str(cell.get('status'))[:9]:<9} "
            f"{common.human_seconds(cell.get('wall_s')):>9} "
            f"{common.human_seconds(cell.get('first_tree_s')):>9} "
            f"{common.human_seconds(cell.get('exact_tree_s')):>9} "
            f"{common.human_bytes(mem.get('bytes')):>10} "
            f"{common.fmt(bench.get('rf'), '.4g'):>6}  {'; '.join(notes)[:60]}"
        )


if __name__ == "__main__":
    common.guard(main)()
