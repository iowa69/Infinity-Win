#!/usr/bin/env python3
"""Write ``REPORT.md`` at the repository root from the benchmark matrix results.

PROJECT.md section 11 is a deliverable, not an appendix: the report has to say,
with numbers, whether each part of the idea pays. Three rules shape this script.

1. **Every result is reported, including the ones where the idea does not pay.**
   Failed configurations, cells with no truth, switches that were accepted by the
   CLI but did not take effect, missing timings, memory that is only a lower
   bound: all of it appears in the report, in a dedicated section, rather than
   being dropped. Section 2 lists "a benchmark that only reports the cases where
   the idea wins" as an anti-goal.
2. **The verdict follows the stated rule and nothing else.** A part of the idea
   pays if turning it off makes time-to-exact or peak memory at least 20% worse on
   at least two of the three datasets. Where a metric is missing the verdict says
   "not measurable" instead of guessing.
3. **An unknown number is never printed as zero.** ``n/a`` means the tool did not
   report it, and the reason is listed next to it.

The key figure -- tree accuracy versus elapsed time, one curve per configuration
-- is read out of the ``events.jsonl`` stream by ``accuracy_over_time.py``, which
this script imports and drives. The numbers behind the figure are always written
as CSV; the PNG is only drawn when matplotlib is importable.

Usage::

    scripts/make_report.py                       # from bench-work/matrix/summary.json
    scripts/make_report.py --matrix-dir DIR --out REPORT.md
    scripts/make_report.py --no-figures
"""

from __future__ import annotations

import argparse
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

import accuracy_over_time as aot  # noqa: E402
import common  # noqa: E402

REPORT_VERSION = 1


# ---------------------------------------------------------------------------
# Gathering measurements
# ---------------------------------------------------------------------------


def first_tree_from_events(run_dir: Path) -> Tuple[Optional[float], Optional[str]]:
    """Time to the first tree, from the event stream, in seconds.

    A fallback for when ``run.json`` does not carry ``timings.first_tree_ms`` yet:
    the first ``backbone_tree`` or ``cluster_focus`` event *is* the first tree the
    user can see.
    """
    path = run_dir / "events.jsonl"
    if not path.exists():
        return None, None
    try:
        log = common.read_event_log(path)
    except common.PrereqError:
        return None, None
    best: Optional[float] = None
    for ev in log.events:
        if common.event_name(ev) not in ("backbone_tree", "cluster_focus"):
            continue
        t_ms = common.event_t_ms(ev)
        if t_ms is None:
            continue
        best = t_ms if best is None else min(best, t_ms)
        break
    if best is None:
        return None, None
    return best / 1000.0, "first tree event in events.jsonl"


def measurements_for_cell(cell: Dict[str, Any]) -> Dict[str, Any]:
    """The numbers the report and the verdict rule use for one cell.

    Times and memory are the **median over the usable repeats**; the full
    ``{n, median, min, max}`` is kept under ``<key>_stats`` for the report's
    ranges. A summary from before ``--repeats`` gives n = 1.
    """
    run = cell.get("run") or {}
    timings_s = run.get("timings_s") or {}
    stats = common.cell_repeat_stats(cell)
    out: Dict[str, Any] = {
        "dataset": cell.get("dataset"),
        "config": cell.get("config"),
        "status": cell.get("status"),
        "reused": bool(cell.get("reused")),
        "wall_s": cell.get("wall_s"),
        "first_tree_s": timings_s.get("first_tree_s"),
        "first_tree_source": "run.json timings.first_tree_ms"
        if timings_s.get("first_tree_s") is not None
        else None,
        "exact_tree_s": timings_s.get("exact_tree_s"),
        "exact_tree_source": "run.json timings.exact_tree_ms"
        if timings_s.get("exact_tree_s") is not None
        else None,
        "peak_memory_bytes": cell.get("peak_memory_bytes"),
        "peak_memory_is_floor": bool(cell.get("peak_memory_is_floor")),
        "memory_source": (cell.get("memory") or {}).get("source"),
        "switch_verified": run.get("switch_verified"),
        "counters": run.get("counters") or {},
        "focus_counts": run.get("focus_counts"),
        "genomes": run.get("genomes"),
        "clusters": run.get("clusters"),
        "cancelled": run.get("cancelled"),
        "bench": cell.get("bench") or {},
        "problems": list(cell.get("problems") or []),
        "error": cell.get("error"),
        "out_dir": cell.get("out_dir"),
        "repeats_total": len(cell.get("repeats") or []) or (1 if cell.get("status") else 0),
    }
    for key in common.REPEAT_METRICS:
        out[f"{key}_stats"] = stats[key]
        if stats[key]["n"]:
            out[key] = stats[key]["median"]
    out["repeats_usable"] = stats["exact_tree_s"]["n"] or stats["wall_s"]["n"]
    if stats["peak_memory_bytes"].get("any_floor"):
        out["peak_memory_is_floor"] = True
    # Fall back to the event stream, then to the process wall time. Both are
    # labelled, because a proxy is not the measurement it stands in for.
    run_dir = Path(cell["out_dir"]) if cell.get("out_dir") else None
    if out["first_tree_s"] is None and run_dir and run_dir.is_dir():
        seconds, source = first_tree_from_events(run_dir)
        if seconds is not None:
            out["first_tree_s"] = seconds
            out["first_tree_source"] = source
    if out["exact_tree_s"] is None and out.get("wall_s") is not None and cell.get("status") in (
        "ok",
        "cached",
    ):
        out["exact_tree_s"] = out["wall_s"]
        out["exact_tree_s_stats"] = dict(out.get("wall_s_stats") or {})
        out["exact_tree_source"] = (
            "process wall time (proxy: run.json carried no timings.exact_tree_ms)"
        )
        out["problems"].append(
            "time-to-exact is the whole `infinity run` wall time, used because "
            "run.json has no timings.exact_tree_ms; it includes output writing"
        )
    bench = out["bench"]
    out["rf"] = bench.get("rf")
    out["rf_normalized"] = bench.get("rf_normalized")
    out["trees_scored"] = bench.get("trees_scored")
    out["distance_mae_snps"] = bench.get("distance_mae_snps")
    out["distance_max_error_snps"] = bench.get("distance_max_error_snps")
    out["pairs_in_truth"] = bench.get("pairs_in_truth")
    out["pairs_compared"] = bench.get("pairs_compared")
    out["pairs_exact"] = bench.get("pairs_exact")
    out["pairs_missing_from_output"] = bench.get("pairs_missing_from_output")
    out["fraction_pairs_exact"] = bench.get("fraction_pairs_exact")
    out["worst_pair"] = bench.get("worst_pair")
    out["cluster_f1"] = bench.get("cluster_f1")
    out["clusters_recovered_exactly"] = bench.get("clusters_recovered_exactly")
    out["clusters_true"] = bench.get("clusters_true")
    out["bench_status"] = bench.get("status")
    out["bench_reason"] = bench.get("reason")
    return out


def resolve_out_dir(cell: Dict[str, Any], matrix_dir: Optional[Path]) -> Optional[str]:
    """The cell's output folder as a path that exists from here, if it can be found.

    ``run_matrix.py`` records ``out_dir`` as it was given, which is relative to
    the folder it was started in. When that does not resolve from the current
    folder, the canonical ``<matrix>/<dataset>/<config>`` layout is tried.
    """
    recorded = cell.get("out_dir")
    if recorded and Path(recorded).is_dir():
        return str(recorded)
    if matrix_dir is not None and cell.get("dataset") and cell.get("config"):
        candidate = Path(matrix_dir) / str(cell["dataset"]) / str(cell["config"])
        if candidate.is_dir():
            return str(candidate)
    return recorded


def collect(
    summary: Dict[str, Any], matrix_dir: Optional[Path] = None
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """``{dataset: {config: measurements}}``."""
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for cell in summary.get("cells", []):
        dataset = cell.get("dataset")
        config = cell.get("config")
        if not dataset or not config:
            continue
        resolved = resolve_out_dir(cell, matrix_dir)
        if resolved != cell.get("out_dir"):
            cell = {**cell, "out_dir": resolved}
        out.setdefault(str(dataset), {})[str(config)] = measurements_for_cell(cell)
    return out


# ---------------------------------------------------------------------------
# Acceptance checks against published expected values
# ---------------------------------------------------------------------------


def _overall_status(checks: Sequence[Any]) -> str:
    """PASS / FAIL / PARTIAL / not measured, from a list of ExpectedCheck."""
    if not checks:
        return "not measured"
    n_pass = sum(1 for c in checks if c.passed is True)
    n_fail = sum(1 for c in checks if c.passed is False)
    n_na = sum(1 for c in checks if c.passed is None)
    if n_pass + n_fail == 0:
        return "not measured"
    if n_fail:
        return "FAIL"
    if n_na:
        return "PARTIAL (some pairs not measured)"
    return "PASS"


def run_expected_checks(
    work: Path,
    summary: Dict[str, Any],
    by_dataset: Dict[str, Dict[str, Dict[str, Any]]],
    expected_dir: Optional[Path] = None,
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Compare Infinity-Win's pairwise SNP counts with published expected values."""
    results: List[Dict[str, Any]] = []
    problems: List[str] = []
    expected_dir = Path(expected_dir) if expected_dir else common.repo_root() / "scripts" / "expected"
    if not expected_dir.is_dir():
        problems.append(f"no expected-values folder at {expected_dir}")
        return results, problems
    expected_files = sorted(expected_dir.glob("*.json"))
    if not expected_files:
        problems.append(f"no expected-values JSON files in {expected_dir}")
    for path in expected_files:
        expected, problem = common.load_json_optional(path)
        if problem or not isinstance(expected, dict):
            problems.append(f"{path.name}: {problem or 'not a JSON object'}")
            continue
        dataset_name = expected.get("dataset")
        entry: Dict[str, Any] = {
            "file": str(path),
            "dataset": dataset_name,
            "title": expected.get("title"),
            "tolerance": expected.get("abs_tolerance_snps"),
            "tolerance_rationale": expected.get("tolerance_rationale"),
            "checks_by_config": {},
        }
        cells = by_dataset.get(dataset_name, {})
        if not cells:
            entry["status"] = "not run"
            entry["problems"] = [
                f"dataset {dataset_name!r} has no cells in the matrix summary"
            ]
            results.append(entry)
            continue
        for config, meas in sorted(cells.items()):
            out_dir = meas.get("out_dir")
            observed: Dict[Tuple[str, str], Optional[float]] = {}
            observed_raw: Dict[Tuple[str, str], Optional[float]] = {}
            observed_filtered: Dict[Tuple[str, str], Optional[float]] = {}
            cell_problems: List[str] = []
            filter_info: Dict[str, Any] = {}
            dist_basis = common.distance_basis([])
            if out_dir and Path(out_dir).is_dir():
                rows, row_problems = common.load_cluster_distances(Path(out_dir))
                cell_problems.extend(row_problems)
                dist_basis = common.distance_basis(rows)
                observed = {
                    common.pair_key(r["a"], r["b"]): r.get("snps") for r in rows
                }
                # Once the engine writes raw_snps, `snps` is its filtered
                # core-genome count; both are checked against the same values.
                observed_raw = {
                    common.pair_key(r["a"], r["b"]): r.get("raw_snps")
                    for r in rows
                    if r.get("raw_snps") is not None
                }
                if not observed:
                    cell_problems.append(
                        "no pairwise distances could be read from this output folder"
                    )
                # Secondary basis: the same comparison after the filter that every
                # published reference number was produced with. Declared in the
                # expected-values file, not hard-coded here.
                basis = (expected.get("comparison_basis") or {}).get("secondary") or {}
                if basis.get("min_pivot_contig_bp") or basis.get("exclude_flags"):
                    recomputed = common.recompute_pairwise_snps(
                        Path(out_dir),
                        min_pivot_contig_bp=basis.get("min_pivot_contig_bp"),
                        exclude_flags=tuple(basis.get("exclude_flags") or ()),
                    )
                    observed_filtered = dict(recomputed.pairs)
                    cell_problems.extend(recomputed.problems)
                    filter_info = {
                        "label": basis.get("label"),
                        "min_pivot_contig_bp": basis.get("min_pivot_contig_bp"),
                        "exclude_flags": list(basis.get("exclude_flags") or []),
                        "variants_total": recomputed.variants_total,
                        "variants_used": recomputed.variants_used,
                        "variants_excluded": recomputed.variants_excluded,
                        "caveat": basis.get("caveat"),
                    }
            else:
                cell_problems.append(f"output folder missing: {out_dir}")
            checks = common.check_expected_pairs(expected, observed)
            checks_raw = (
                common.check_expected_pairs(expected, observed_raw) if observed_raw else []
            )
            checks_filtered = (
                common.check_expected_pairs(expected, observed_filtered)
                if observed_filtered
                else []
            )
            n_pass = sum(1 for c in checks if c.passed is True)
            n_fail = sum(1 for c in checks if c.passed is False)
            n_na = sum(1 for c in checks if c.passed is None)
            entry["checks_by_config"][config] = {
                "checks": [c.as_dict() for c in checks],
                "passed": n_pass,
                "failed": n_fail,
                "not_measured": n_na,
                "status": (
                    "not measured" if n_pass + n_fail == 0
                    else ("PASS" if n_fail == 0 and n_na == 0 else
                          ("FAIL" if n_fail else "PARTIAL (some pairs not measured)"))
                ),
                "snps_label": dist_basis["snps_label"],
                "checks_raw": [c.as_dict() for c in checks_raw],
                "status_raw": _overall_status(checks_raw),
                "checks_filtered": [c.as_dict() for c in checks_filtered],
                "status_filtered": _overall_status(checks_filtered),
                "filter": filter_info,
                "problems": cell_problems,
            }
        results.append(entry)
    return results, problems


# ---------------------------------------------------------------------------
# Accuracy-over-time figure
# ---------------------------------------------------------------------------


def build_accuracy_curves(
    by_dataset: Dict[str, Dict[str, Dict[str, Any]]],
    work: Path,
    figures_dir: Path,
    *,
    collapse_snps: float,
    max_points: int,
    enabled: bool = True,
    summary: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Drive ``accuracy_over_time.py`` over every usable cell.

    Truth is looked up in ``<work>/<dataset>/dataset.json`` and then in the
    matrix summary, which is the only record of where a simulated dataset kept
    outside the work folder has its truth.
    """
    doc: Dict[str, Any] = {
        "figures_dir": str(figures_dir),
        "datasets": {},
        "problems": [],
        "enabled": enabled,
    }
    if not enabled:
        doc["problems"].append("figure generation disabled with --no-figures")
        return doc
    common.ensure_dir(figures_dir)
    all_series: Dict[str, List[Dict[str, Any]]] = {}
    labels = aot.labels_from_summary(summary)
    empty_reasons: Dict[str, str] = {}
    for dataset, cells in sorted(by_dataset.items()):
        truth_tree, truth_dir, truth_origin = common.dataset_truth(work, dataset, summary)
        per_dataset: Dict[str, Any] = {
            "runs": {},
            "truth_tree": str(truth_tree) if truth_tree else None,
            "truth_dir": str(truth_dir) if truth_dir else None,
            "truth_origin": truth_origin,
            "simulator_truth": bool(truth_dir and common.is_simulator_truth_dir(truth_dir)),
            "label": labels.get(dataset, dataset),
        }
        # Why this dataset may have no accuracy curve at all, stated once and
        # precisely: no reference declared, and/or too few taxa for RF to exist.
        genomes = next((m.get("genomes") for m in cells.values() if m.get("genomes")), None)
        truth_note = None
        if truth_tree is None and truth_dir is None:
            truth_note = (
                "this dataset declares no reference tree and no simulator truth "
                "folder, so tree accuracy cannot be scored"
            )
            if genomes is not None and genomes < 4:
                truth_note += (
                    f"; with {genomes} genomes the Robinson-Foulds distance is "
                    "undefined in any case (it needs at least four taxa), so this "
                    "dataset contributes distance accuracy and the acceptance "
                    "check only"
                )
        per_dataset["truth_note"] = truth_note
        series: Dict[str, List[Dict[str, Any]]] = {}
        for config, meas in sorted(cells.items()):
            out_dir = meas.get("out_dir")
            if not out_dir or not Path(out_dir).is_dir():
                per_dataset["runs"][config] = {"status": "skipped", "reason": "no output folder"}
                continue
            try:
                res = aot.analyse_run_dir(
                    Path(out_dir),
                    dataset=dataset,
                    config=config,
                    out_dir=figures_dir,
                    truth_tree=truth_tree,
                    truth_dir=truth_dir,
                    rf_from=None,
                    rf_mode="auto",
                    collapse_snps=collapse_snps,
                    collapse_mode="auto",
                    collapse_length=None,
                    comparable_bp=None,
                    max_points=max_points,
                    infinity=None,
                )
            except common.PrereqError as exc:
                res = {"status": "failed", "problems": [str(exc)], "points": 0, "rows": []}
            rows = res.pop("rows", [])
            res.pop("dataset", None)
            res.pop("config", None)
            if rows and meas.get("status") == "failed":
                # The event stream of a run that was killed or crashed still
                # scores, but the curve stops where the run stopped.
                res["status"] = "partial (run failed)"
                res.setdefault("problems", []).insert(
                    0, f"the run failed ({meas.get('error') or 'see section 9'}); the curve "
                    "covers only the events written before it stopped"
                )
            if rows:
                ordered = sorted(
                    (r for r in rows if r.get("t_ms") is not None),
                    key=lambda r: float(r["t_ms"]),
                )
                if ordered:
                    # "first" is the first point whose RF is actually defined:
                    # while only part of the tree exists RF can be undefined
                    # (fewer than four shared taxa), and that must not be shown
                    # as an accuracy of zero.
                    defined = [r for r in ordered if r.get("rf_normalized") is not None]
                    res["points_without_rf"] = len(ordered) - len(defined)
                    if defined:
                        res["first"] = {
                            "t_s": float(defined[0]["t_ms"]) / 1000.0,
                            "rf": defined[0].get("rf"),
                            "rf_normalized": defined[0].get("rf_normalized"),
                        }
                        res["last"] = {
                            "t_s": float(defined[-1]["t_ms"]) / 1000.0,
                            "rf": defined[-1].get("rf"),
                            "rf_normalized": defined[-1].get("rf_normalized"),
                            "rf_between_clusters": defined[-1].get("rf_between_clusters"),
                            "rf_within_clusters": defined[-1].get("rf_within_clusters"),
                        }
                    res["t_span_s"] = [
                        float(ordered[0]["t_ms"]) / 1000.0,
                        float(ordered[-1]["t_ms"]) / 1000.0,
                    ]
            if truth_note and not res.get("points"):
                # Put the dataset-level reason first: it is the one that matters.
                res["problems"] = [truth_note] + [
                    p for p in res.get("problems", []) if p != truth_note
                ]
            per_dataset["runs"][config] = res
            if rows:
                series[config] = rows
                all_series[f"{dataset}/{config}"] = rows
        if not series:
            first_problem = next(
                (r.get("problems")[0] for r in per_dataset["runs"].values() if r.get("problems")),
                "no configuration produced accuracy points",
            )
            empty_reasons[dataset] = truth_note or first_problem
        if series:
            png = figures_dir / f"accuracy_over_time__{dataset}.png"
            written, reason = aot.plot_curves(
                series, png, title=f"Tree accuracy over time - {dataset}",
                labels=labels, default_dataset=dataset,
            )
            per_dataset["figure"] = str(png) if written else None
            per_dataset["figure_error"] = None if written else reason
            if not written:
                doc["problems"].append(f"{dataset}: no figure ({reason}); CSVs were still written")
        else:
            per_dataset["figure"] = None
            per_dataset["figure_error"] = "no configuration produced accuracy points"
        per_dataset["csv"] = str(figures_dir / "accuracy_over_time.csv")
        doc["datasets"][dataset] = per_dataset
    if all_series:
        combined = figures_dir / "accuracy_over_time.csv"
        rows = [r for series in all_series.values() for r in series]
        aot.write_rows(combined, rows)
        doc["combined_csv"] = str(combined)
        doc["combined_rows"] = len(rows)
        png = figures_dir / "accuracy_over_time__all.png"
        written, reason = aot.plot_curves(
            all_series, png, title="Tree accuracy over time - every dataset",
            labels=labels, empty=empty_reasons,
        )
        doc["combined_figure"] = str(png) if written else None
        if not written:
            doc["problems"].append(f"combined figure not written: {reason}")
    return doc


def _residual(last: Dict[str, Any]) -> str:
    """``"RF (between / within clusters)"`` for the last point of a curve."""
    if last.get("rf") in (None, ""):
        return "n/a"
    between, within = last.get("rf_between_clusters"), last.get("rf_within_clusters")
    if between in (None, "") or within in (None, ""):
        return common.fmt(last.get("rf"), ".4g")
    return f"{common.fmt(last.get('rf'), '.4g')} ({between} / {within})"


def curve_endpoints(doc: Dict[str, Any]) -> List[List[Any]]:
    """A text stand-in for the figure: first and last defined point of every curve."""
    rows: List[List[Any]] = []
    for dataset, info in sorted((doc.get("datasets") or {}).items()):
        for config, run in sorted((info.get("runs") or {}).items()):
            points = run.get("points") or 0
            first = run.get("first") or {}
            last = run.get("last") or {}
            # A curve with no points must say why, not just look blank.
            why = ""
            if not points:
                problems = run.get("problems") or []
                why = problems[0] if problems else run.get("reason") or "no reason recorded"
            rows.append(
                [
                    dataset,
                    config,
                    run.get("status"),
                    points,
                    run.get("points_without_rf"),
                    run.get("rf_source") or "n/a",
                    common.fmt((run.get("collapse") or {}).get("threshold"), ".3g"),
                    common.fmt(first.get("t_s"), ".4g"),
                    common.fmt(first.get("rf_normalized"), ".4f"),
                    common.fmt(last.get("t_s"), ".4g"),
                    common.fmt(last.get("rf_normalized"), ".4f"),
                    _residual(last),
                    why[:200],
                ]
            )
    return rows


# ---------------------------------------------------------------------------
# Baselines
# ---------------------------------------------------------------------------


def collect_baselines(work: Path) -> Dict[str, Any]:
    """Read whatever ``baselines.py`` left behind."""
    base = work / "baselines"
    doc: Dict[str, Any] = {"dir": str(base), "datasets": {}, "problems": []}
    if not base.is_dir():
        doc["problems"].append(
            f"no baselines folder at {base}; run scripts/baselines.py all --dataset <name> "
            "--run-dir <matrix output folder>"
        )
        return doc
    for dataset_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        entry: Dict[str, Any] = {}
        for name in ("ska2.json", "mash.json"):
            path = dataset_dir / name
            if path.exists():
                obj, problem = common.load_json_optional(path)
                entry[name.split(".")[0]] = obj if obj is not None else {"error": problem}
        for path in sorted(dataset_dir.glob("compare__*.json")):
            obj, problem = common.load_json_optional(path)
            entry.setdefault("compare", {})[path.stem] = obj if obj is not None else {"error": problem}
        doc["datasets"][dataset_dir.name] = entry
    if not doc["datasets"]:
        doc["problems"].append(f"{base} holds no dataset folders")
    return doc


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------


def rel(path: Optional[str], root: Path) -> str:
    if not path:
        return "n/a"
    p = Path(path)
    try:
        return str(p.relative_to(root))
    except ValueError:
        return str(p)


def variants_fingerprint(outdir: Path) -> Dict[str, Any]:
    """A comparable fingerprint of every variant row in an output folder.

    All ``cluster_<id>.variants.tsv`` rows are pooled and sorted, *ignoring which
    cluster file each row came from*: cluster numbering and the number of clusters
    can legitimately shift between configurations (``--single-cluster`` collapses
    everything into one), so keying by file name would report a difference where
    the variant set is the same. Pooling gives the invariant that actually
    matters - the same set of variants overall.

    Used to tell two kinds of disagreement apart: a switch that changes the
    *variants* contradicts DECISIONS.md D16, while one that leaves them identical
    and only moves the distances is the coverage-bookkeeping difference D16
    records.
    """
    outdir = Path(outdir)
    rows: List[str] = []
    files = 0
    problems: List[str] = []
    if not outdir.is_dir():
        return {"rows": [], "files": 0, "problems": [f"no such folder: {outdir}"]}
    for path in sorted(outdir.glob("cluster_*.variants.tsv")):
        files += 1
        try:
            table = common.read_table(path)
        except common.PrereqError as exc:
            problems.append(str(exc))
            continue
        for row in table.rows:
            rows.append("|".join(row.get(col, "") for col in table.header))
    rows.sort()
    return {"rows": rows, "files": files, "problems": problems}


def diff_variants(base: Dict[str, Any], other: Dict[str, Any]) -> Dict[str, Any]:
    """Multiset difference between two variant fingerprints."""
    from collections import Counter

    base_counts = Counter(base.get("rows") or [])
    other_counts = Counter(other.get("rows") or [])
    only_base = sum((base_counts - other_counts).values())
    only_other = sum((other_counts - base_counts).values())
    return {
        "identical": base_counts == other_counts,
        "rows_baseline": sum(base_counts.values()),
        "rows_config": sum(other_counts.values()),
        "rows_only_in_baseline": only_base,
        "rows_only_in_config": only_other,
        "files_baseline": base.get("files"),
        "files_config": other.get("files"),
        "problems": list(base.get("problems") or []) + list(other.get("problems") or []),
    }


#: Switches DECISIONS.md D16 states are exactly result-neutral. A distance
#: difference from any of these is an undocumented defect.
MUST_BE_EXACTLY_NEUTRAL = ("flat", "no_piece_cache", "no_profile_grouping")
#: DECISIONS.md D16: `--no-block-skip` keeps the variants identical but changes
#: the coverage bookkeeping, so its distances may legitimately differ. That is a
#: finding for section 11, not a defect - but it must be reported, never hidden.
DOCUMENTED_NON_NEUTRAL = ("no_block_skip",)
#: `--single-cluster` changes the pivot, so its pair set and coordinates differ
#: by construction.
EXPECTED_TO_DIFFER = ("single_cluster",)


def check_switch_invariance(
    by_dataset: Dict[str, Dict[str, Dict[str, Any]]]
) -> Dict[str, Any]:
    """Test PROJECT.md section 10.2 property 6 at the benchmark level.

    "Flat mode, hierarchical mode, and every cache switch give identical
    variants, distances, and trees." Every switch is supposed to change only how
    much work is done, never the answer. Each configuration's pairwise SNP counts
    and variant lists are compared with the baseline's, and the result is
    classified as one of:

    ``identical``            reproduced the baseline exactly
    ``DEFECT``               differs, and DECISIONS.md does not say it may
    ``documented deviation`` differs, as DECISIONS.md D16 records for
                             ``--no-block-skip`` (variants identical, coverage not)
    ``undocumented``         ``--no-block-skip`` changed the *variants*, which D16
                             says must stay identical - a new problem
    ``expected``             ``--single-cluster``, which changes the pivot
    """
    out: Dict[str, Any] = {}
    for dataset, cells in sorted(by_dataset.items()):
        baseline = cells.get("baseline") or {}
        baseline_dir = baseline.get("out_dir")
        entry: Dict[str, Any] = {"configs": {}, "problems": []}
        if not baseline_dir or not Path(baseline_dir).is_dir():
            entry["problems"].append(
                f"no baseline output folder to compare against ({baseline_dir})"
            )
            out[dataset] = entry
            continue
        base_rows, base_problems = common.load_cluster_distances(Path(baseline_dir))
        entry["problems"].extend(f"baseline: {p}" for p in base_problems)
        base_index = {
            common.pair_key(r["a"], r["b"]): r.get("snps") for r in base_rows
        }
        base_raw = {
            common.pair_key(r["a"], r["b"]): r.get("raw_snps")
            for r in base_rows
            if r.get("raw_snps") is not None
        }
        base_variants = variants_fingerprint(Path(baseline_dir))
        for config, meas in sorted(cells.items()):
            if config == "baseline":
                continue
            other_dir = meas.get("out_dir")
            record: Dict[str, Any] = {"status": "not run"}
            if not other_dir or not Path(other_dir).is_dir():
                record["problems"] = [f"no output folder ({other_dir})"]
                entry["configs"][config] = record
                continue
            rows, problems = common.load_cluster_distances(Path(other_dir))
            record["problems"] = list(problems)
            index = {common.pair_key(r["a"], r["b"]): r.get("snps") for r in rows}
            raw_index = {
                common.pair_key(r["a"], r["b"]): r.get("raw_snps")
                for r in rows
                if r.get("raw_snps") is not None
            }
            shared = sorted(set(base_index) & set(index))
            differing = []
            max_diff = 0
            for key in shared:
                a, b = base_index[key], index[key]
                if a is None or b is None:
                    continue
                if a != b:
                    differing.append((key, a, b))
                    max_diff = max(max_diff, abs(int(a) - int(b)))
            # Property 6 covers every distance the engine writes, so once it
            # reports raw counts next to the filtered ones, both must match.
            raw_differing = [
                (key, base_raw[key], raw_index[key])
                for key in sorted(set(base_raw) & set(raw_index))
                if base_raw[key] != raw_index[key]
            ]
            other_variants = variants_fingerprint(Path(other_dir))
            vdiff = diff_variants(base_variants, other_variants)
            variants_identical = bool(vdiff["identical"])
            record["problems"].extend(vdiff["problems"])
            record.update(
                {
                    "pairs_compared": len(shared),
                    "pairs_only_in_baseline": len(set(base_index) - set(index)),
                    "pairs_only_in_config": len(set(index) - set(base_index)),
                    "pairs_differing": len(differing),
                    "max_abs_diff_snps": max_diff if differing else None,
                    "raw_pairs_differing": len(raw_differing) if (base_raw or raw_index) else None,
                    "variants_identical": variants_identical,
                    "variant_rows_baseline": vdiff["rows_baseline"],
                    "variant_rows_config": vdiff["rows_config"],
                    "variant_rows_only_in_baseline": vdiff["rows_only_in_baseline"],
                    "variant_rows_only_in_config": vdiff["rows_only_in_config"],
                    "examples": [
                        {"pair": list(k), "baseline": a, "config": b}
                        for k, a, b in differing[:5]
                    ],
                }
            )
            if raw_differing:
                record["examples"].extend(
                    {"pair": list(k), "baseline": a, "config": b, "basis": "raw"}
                    for k, a, b in raw_differing[: max(0, 5 - len(record["examples"]))]
                )
                differing = differing + raw_differing
            if not shared:
                record["status"] = "no comparable pairs"
            elif not differing and variants_identical:
                record["status"] = "identical"
            elif config in EXPECTED_TO_DIFFER:
                record["status"] = "expected"
            elif config in DOCUMENTED_NON_NEUTRAL:
                # D16 records a distance/coverage difference with the variant list
                # unchanged. If the variants moved too, that is beyond D16.
                record["status"] = (
                    "documented deviation" if variants_identical else "UNDOCUMENTED"
                )
            else:
                # Property 6 covers variants, distances and trees, so a variant-only
                # difference is just as much a failure as a distance difference.
                record["status"] = "DEFECT"
                if not differing:
                    record["problems"].append(
                        "pairwise distances match the baseline but the variant list "
                        "does not; property 6 requires both to be identical"
                    )
            entry["configs"][config] = record
        out[dataset] = entry
    return out


def render_invariance(doc: Dict[str, Any]) -> str:
    out: List[str] = []
    if not doc:
        return "_No datasets to compare._"
    marks = {
        "identical": "identical",
        "DEFECT": "**DEFECT**",
        "documented deviation": "documented deviation (D16)",
        "UNDOCUMENTED": "**UNDOCUMENTED**",
        "expected": "expected (different pivot)",
        "no comparable pairs": "no comparable pairs",
        "not run": "not run",
    }
    for dataset, entry in sorted(doc.items()):
        out.append(f"### {dataset}")
        out.append("")
        rows: List[List[Any]] = []
        for config, rec in sorted((entry.get("configs") or {}).items()):
            status = rec.get("status")
            variant_note = {True: "yes", False: "**no**", None: "n/a"}[
                rec.get("variants_identical")
            ]
            if rec.get("variants_identical") is False:
                variant_note += (
                    f" ({common.fmt(rec.get('variant_rows_only_in_baseline'))} rows only "
                    f"in baseline, {common.fmt(rec.get('variant_rows_only_in_config'))} "
                    f"only in `{config}`; {common.fmt(rec.get('variant_rows_baseline'))} "
                    f"vs {common.fmt(rec.get('variant_rows_config'))} total)"
                )
            rows.append(
                [
                    config,
                    marks.get(status, status),
                    common.fmt(rec.get("pairs_compared")),
                    common.fmt(rec.get("pairs_differing")),
                    common.fmt(rec.get("max_abs_diff_snps"), ".4g"),
                    common.fmt(rec.get("raw_pairs_differing")),
                    variant_note,
                    common.fmt(rec.get("pairs_only_in_baseline")),
                    common.fmt(rec.get("pairs_only_in_config")),
                ]
            )
        out.append(
            common.markdown_table(
                ["config", "distances vs baseline", "pairs compared",
                 "pairs differing", "max |diff| SNPs", "raw pairs differing",
                 "variant lists identical", "pairs only in baseline", "pairs only in config"],
                rows,
            )
        )
        out.append("")
        for config, rec in sorted((entry.get("configs") or {}).items()):
            if rec.get("examples"):
                out.append(f"Examples of differing pairs under `{config}`:")
                out.append("")
                for ex in rec["examples"]:
                    basis = " (raw_snps)" if ex.get("basis") == "raw" else ""
                    out.append(
                        f"- {ex['pair'][0]} vs {ex['pair'][1]}{basis}: baseline {ex['baseline']} "
                        f"SNPs, `{config}` {ex['config']} SNPs"
                    )
                out.append("")
        for problem in entry.get("problems") or []:
            out.append(f"> {problem}")
    return "\n".join(out)


def verdict_mark(passed: Optional[bool]) -> str:
    """Bold PASS/FAIL, or an unbolded 'not measured' - never a silent pass."""
    if passed is True:
        return "**PASS**"
    if passed is False:
        return "**FAIL**"
    return "not measured"


def fmt_stat(m: Dict[str, Any], key: str, *, memory: bool = False) -> str:
    """``median [min-max]`` over the repeats, or the single value when n = 1."""
    human = common.human_bytes if memory else common.human_seconds
    stats = m.get(f"{key}_stats") or {}
    value = m.get(key)
    if value is None:
        return "n/a"
    if stats.get("n", 0) > 1 and stats.get("min") is not None:
        return f"{human(value)} [{human(stats['min'])}-{human(stats['max'])}]"
    return human(value)


def render_accuracy_table(cells: Dict[str, Dict[str, Any]]) -> str:
    """Accuracy against truth, per configuration.

    The engine scores each cluster tree against its own truth tree, so `RF` is
    the **sum** over clusters and `RF (norm)` their **mean**; both are labelled.
    A run with no scorable truth shows `n/a` and the reason, never a zero.
    """
    header = [
        "config", "bench", "clusters scored", "RF (sum over clusters)",
        "RF (mean normalised)", "pairs in truth", "pairs compared",
        "pairs exact", "fraction exact", "max |err| SNPs", "mean |err| SNPs",
        "cluster F1", "clusters recovered",
    ]
    rows: List[List[Any]] = []
    order = [c.name for c in common.CONFIGS if c.name in cells]
    order += [c for c in sorted(cells) if c not in order]
    for config in order:
        m = cells[config]
        status = m.get("bench_status") or "n/a"
        if status == "skipped":
            status = "skipped (no simulator truth)"
        recovered = None
        if m.get("clusters_recovered_exactly") is not None and m.get("clusters_true") is not None:
            recovered = f"{m['clusters_recovered_exactly']}/{m['clusters_true']}"
        rows.append(
            [
                config,
                status,
                common.fmt(m.get("trees_scored")),
                common.fmt(m.get("rf"), ".4g"),
                common.fmt(m.get("rf_normalized"), ".4f"),
                common.fmt(m.get("pairs_in_truth")),
                common.fmt(m.get("pairs_compared")),
                common.fmt(m.get("pairs_exact")),
                common.fmt(m.get("fraction_pairs_exact"), ".4f"),
                common.fmt(m.get("distance_max_error_snps"), ".4g"),
                common.fmt(m.get("distance_mae_snps"), ".4g"),
                common.fmt(m.get("cluster_f1"), ".4f"),
                recovered if recovered else "n/a",
            ]
        )
    return common.markdown_table(header, rows)


def render_dataset_table(cells: Dict[str, Dict[str, Any]]) -> str:
    header = [
        "config", "status", "runs", "wall", "time to first tree", "time to exact tree",
        "peak memory", "block identity", "pieces compared", "cache hits",
        "alignments", "bases aligned", "breakpoints", "switch applied",
    ]
    rows: List[List[Any]] = []
    order = [c.name for c in common.CONFIGS if c.name in cells]
    order += [c for c in sorted(cells) if c not in order]
    for config in order:
        m = cells[config]
        counters = m.get("counters") or {}
        matched = counters.get("blocks_matched") or []
        matched_total = sum(v for v in matched if isinstance(v, int)) if matched else None
        looked = counters.get("blocks_looked_up")
        identity = (
            matched_total / looked
            if isinstance(looked, int) and looked and matched_total is not None
            else None
        )
        status = m.get("status")
        if status == "failed":
            status = f"FAILED ({m.get('error') or 'see failures'})"
        runs = f"{common.fmt(m.get('repeats_usable'))}/{common.fmt(m.get('repeats_total'))}"
        rows.append(
            [
                config,
                status,
                runs,
                fmt_stat(m, "wall_s"),
                fmt_stat(m, "first_tree_s"),
                fmt_stat(m, "exact_tree_s"),
                fmt_stat(m, "peak_memory_bytes", memory=True)
                + (" (floor)" if m.get("peak_memory_is_floor") else ""),
                common.fmt(identity, ".3%") if identity is not None else "n/a",
                common.fmt(counters.get("pieces_compared")),
                common.fmt(counters.get("cache_hits")),
                common.fmt(counters.get("alignments_run")),
                common.fmt(counters.get("bases_aligned")),
                common.fmt(counters.get("breakpoints")),
                {True: "yes", False: "**NO**", None: "n/a"}[m.get("switch_verified")],
            ]
        )
    return common.markdown_table(header, rows)


def render_verdict(verdicts: Sequence[common.Verdict]) -> Tuple[str, str]:
    """The verdict summary line per switch, and the detailed per-metric tables."""
    summary_rows: List[List[Any]] = []
    detail: List[str] = []
    VERDICT_LABELS = {
        "pays": "**pays**",
        "does_not_pay": "does not pay",
        "undetermined": "undetermined",
        "not_measurable": "not measurable",
    }
    for v in verdicts:
        label = VERDICT_LABELS.get(v.state, v.state)
        summary_rows.append(
            [
                f"`{v.switch}`",
                v.removes or "-",
                label,
                f"{len(v.datasets_worse)} of {v.datasets_total}",
                f"{len(v.datasets_measurable)} of {v.datasets_total}",
                v.reason,
            ]
        )
        detail.append(f"#### `{v.switch}` - {v.removes}")
        detail.append("")
        rows = []
        for m in v.metrics:
            human = common.human_bytes if m.metric.startswith("peak memory") else common.human_seconds

            def side(value: Optional[float], rng: Any, n: Optional[int]) -> str:
                if value is None:
                    return "n/a"
                text = human(value)
                if n and n > 1 and rng and None not in rng:
                    text += f" [{human(rng[0])}-{human(rng[1])}]"
                return text + (f" (n={n})" if n else "")

            if m.state == "too_short":
                worse = "too short to judge"
            elif m.state == "excluded":
                worse = "excluded"
            else:
                worse = {True: "**yes**", False: "no", None: "n/a"}[m.worse]
            overlap = {True: "yes", False: "no", None: "n/a"}[m.ranges_overlap]
            rows.append(
                [
                    m.dataset,
                    m.metric,
                    side(m.baseline, m.baseline_range, m.n_baseline),
                    side(m.switched_off, m.switched_off_range, m.n_switched_off),
                    common.fmt(m.ratio, ".3f") + "x" if m.ratio is not None else "n/a",
                    common.fmt(common.pct_change(m.switched_off, m.baseline), "+.1f") + "%"
                    if m.ratio is not None
                    else "n/a",
                    overlap,
                    worse,
                    m.note if m.state not in ("too_short", "excluded") else "",
                ]
            )
        detail.append(
            common.markdown_table(
                ["dataset", "metric", "baseline median [min-max]",
                 "switch off median [min-max]", "ratio of medians", "change",
                 "ranges overlap", f">= {int(round((v.threshold - 1) * 100))}% worse?"],
                [r[:8] for r in rows],
            )
        )
        notes = [r[8] for r in rows if r[8]]
        if notes:
            detail.append("")
            detail.append("Notes: " + "; ".join(sorted(set(notes))))
        detail.append("")
        detail.append(f"**Verdict: {label.replace('**', '')}.** {v.reason}")
        detail.append("")
    summary = common.markdown_table(
        ["switch", "what turning it off removes", "verdict", "datasets >=20% worse",
         "datasets judged", "why"],
        summary_rows,
    )
    return summary, "\n".join(detail)


def render_expected(entries: Sequence[Dict[str, Any]]) -> str:
    out: List[str] = []
    if not entries:
        out.append("_No expected-values files were found, so no acceptance check was run._")
        return "\n".join(out)
    for entry in entries:
        out.append(f"### {entry.get('dataset')} - {entry.get('title') or entry['file']}")
        out.append("")
        out.append(f"Source file: `{rel(entry['file'], common.repo_root())}`")
        if entry.get("tolerance") is not None:
            out.append(
                f"Tolerance: +/- {entry['tolerance']} SNP. {entry.get('tolerance_rationale') or ''}"
            )
        if entry.get("status") == "not run":
            out.append("")
            out.append(
                "**Not run:** " + "; ".join(entry.get("problems") or [])
            )
            out.append("")
            continue
        by_config = entry.get("checks_by_config") or {}
        if not by_config:
            out.append("")
            out.append("_No configuration produced a check._")
            out.append("")
            continue
        # The bases: what the engine reported (its filtered core-genome count and,
        # once it writes one, its raw count) and the scripts' own recomputation
        # after the filter every published reference number was produced with.
        any_raw = any((info.get("checks_raw") or []) for info in by_config.values())
        any_filtered = any((info.get("checks_filtered") or []) for info in by_config.values())
        if any_raw:
            out.append("")
            out.append(
                "The engine reports two bases in `cluster_<id>.dist.tsv`: **raw** "
                "(`raw_snps`, every SNP between the two assemblies) and **core** "
                "(`snps`, its filtered core-genome count). Both are checked against "
                "the same published values; the published numbers are themselves "
                "filtered, so the core count is the like-for-like comparison and the "
                "raw count shows what the filter removed."
            )
        if any_filtered:
            filt = next(
                (info.get("filter") or {} for info in by_config.values() if info.get("filter")),
                {},
            )
            out.append("")
            out.append(
                f"**Script-filtered** is "
                f"{filt.get('label') or 'a filtered recomputation'}"
                + (
                    f" (pivot contigs >= {filt.get('min_pivot_contig_bp')} bp)"
                    if filt.get("min_pivot_contig_bp")
                    else ""
                )
                + ", recomputed from `cluster_<id>.variants.tsv` by the benchmark "
                "scripts, independently of the engine. "
                + (
                    "It is kept as a cross-check of the engine's own filter."
                    if any_raw
                    else "It exists because this engine build reports only an "
                    "unfiltered count, while every published reference number here "
                    "*is* filtered - so the unfiltered comparison is not apples to "
                    "apples."
                )
            )
            if filt.get("caveat"):
                out.append("")
                out.append(f"> {filt['caveat']}")
        rows = []
        for config, info in sorted(by_config.items()):
            checks = info.get("checks") or []

            def index(entries: Sequence[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
                return {
                    (tuple(c.get("pair") or []) and "-".join(c["pair"])) or c.get("name"): c
                    for c in entries
                }

            filtered = index(info.get("checks_filtered") or [])
            raw = index(info.get("checks_raw") or [])
            for check in checks:
                pair = "-".join(check.get("pair") or []) or check.get("name")
                row = [
                    config,
                    pair,
                    common.fmt(check.get("expected"), ".4g"),
                    common.fmt(check.get("observed"), ".4g"),
                    common.fmt(
                        None
                        if check.get("observed") is None or check.get("expected") is None
                        else float(check["observed"]) - float(check["expected"]),
                        "+.4g",
                    ),
                    verdict_mark(check.get("passed")),
                ]
                if any_raw:
                    rcheck = raw.get(pair) or {}
                    row += [common.fmt(rcheck.get("observed"), ".4g"), verdict_mark(rcheck.get("passed"))]
                if any_filtered:
                    fcheck = filtered.get(pair) or {}
                    row += [common.fmt(fcheck.get("observed"), ".4g"), verdict_mark(fcheck.get("passed"))]
                row += [f"+/- {common.fmt(check.get('tolerance'), '.4g')}", check.get("source") or ""]
                rows.append(row)
            filt_info = info.get("filter") or {}
            note = "; ".join(info.get("problems") or [])[:160]
            if any_filtered and filt_info:
                note = (
                    f"variants {common.fmt(filt_info.get('variants_used'))} used / "
                    f"{common.fmt(filt_info.get('variants_excluded'))} excluded by the "
                    f"script filter. {note}"
                ).strip()
            overall = [f"`{config}` overall", "", "", "", "", f"**{info.get('status')}**"]
            if any_raw:
                overall += ["", f"**{info.get('status_raw')}**"]
            if any_filtered:
                overall += ["", f"**{info.get('status_filtered')}**"]
            overall += ["", note]
            rows.append(overall)
        out.append("")
        snps_header = (
            "engine core (`snps`)" if any_raw else "as reported (`snps`)"
        )
        head = ["config", "pair", "expected", snps_header, "diff", "verdict"]
        if any_raw:
            head += ["engine raw (`raw_snps`)", "verdict"]
        if any_filtered:
            head += ["script-filtered", "verdict"]
        head += ["tolerance", "provenance / note"]
        out.append(common.markdown_table(head, rows))
        out.append("")
    return "\n".join(out)


def render_baselines(doc: Dict[str, Any], root: Path) -> str:
    out: List[str] = []
    datasets = doc.get("datasets") or {}
    if not datasets:
        out.append("_No baseline results found._")
        for problem in doc.get("problems", []):
            out.append(f"> {problem}")
        return "\n".join(out)
    for name, entry in sorted(datasets.items()):
        out.append(f"### {name}")
        out.append("")
        rows: List[List[Any]] = []
        ska = entry.get("ska2")
        if isinstance(ska, dict):
            stats = ska.get("stats") or {}
            rows.append(
                [
                    "SKA2 (exact, alignment-free)",
                    ska.get("status"),
                    "ran here" if ska.get("ran_here") else "pre-computed",
                    f"{common.fmt(stats.get('rows'))} rows / {common.fmt(stats.get('sample_count'))} samples",
                    f"min {common.fmt(stats.get('min'))}, median {common.fmt(stats.get('median'))}, "
                    f"max {common.fmt(stats.get('max'))}",
                    rel(ska.get("source"), root),
                ]
            )
        mash = entry.get("mash")
        if isinstance(mash, dict):
            rows.append(
                [
                    "Mash + neighbour-joining (fast, blurry)",
                    mash.get("status"),
                    "ran here" if mash.get("ran_here") else "pre-computed",
                    f"{common.fmt(mash.get('rows'))} rows",
                    f"median {common.fmt(mash.get('distance_median'), '.4g')}",
                    "; ".join(mash.get("notes") or [])[:160] or rel(mash.get("table"), root),
                ]
            )
        for key, cmp_doc in sorted((entry.get("compare") or {}).items()):
            if not isinstance(cmp_doc, dict):
                continue
            summary = cmp_doc.get("summary") or {}
            rows.append(
                [
                    f"Infinity-Win vs SKA2 ({cmp_doc.get('run_config') or key})",
                    cmp_doc.get("status"),
                    "-",
                    f"{common.fmt(summary.get('pairs_in_both'))} pairs in both "
                    f"(only Infinity-Win: {common.fmt(summary.get('pairs_only_in_infinity'))}, "
                    f"only SKA2: {common.fmt(summary.get('pairs_only_in_ska2'))})",
                    f"mean |diff| {common.fmt(summary.get('mean_abs_diff_snps'), '.4g')} SNPs, "
                    f"max {common.fmt(summary.get('max_abs_diff_snps'), '.4g')}, "
                    f"r={common.fmt(summary.get('pearson'), '.4f')}",
                    rel(cmp_doc.get("tsv"), root),
                ]
            )
            tree_cmp = cmp_doc.get("tree_vs_reference") or {}
            if tree_cmp:
                rows.append(
                    [
                        f"  tree vs reference ({cmp_doc.get('run_config') or key})",
                        tree_cmp.get("status"),
                        "-",
                        f"RF {common.fmt(tree_cmp.get('rf'), '.4g')} "
                        f"(norm {common.fmt(tree_cmp.get('rf_normalized'), '.4f')}), "
                        f"{common.fmt(tree_cmp.get('leaves_shared'))} shared leaves",
                        f"collapse {common.fmt(tree_cmp.get('collapse_threshold'), '.3g')}",
                        tree_cmp.get("reason") or rel(tree_cmp.get("reference_tree"), root),
                    ]
                )
        out.append(
            common.markdown_table(
                ["baseline / comparison", "status", "how obtained", "size", "headline numbers", "file / note"],
                rows,
            )
        )
        out.append("")
    for problem in doc.get("problems", []):
        out.append(f"> {problem}")
    return "\n".join(out)


def _mean_genome_mb(cells_for: Dict[str, Dict[str, Any]], summary_cells: Sequence[Dict[str, Any]], name: str) -> Optional[float]:
    for cell in summary_cells:
        run = cell.get("run") or {}
        if str(cell.get("dataset")) == name and run.get("total_length") and run.get("genomes"):
            return float(run["total_length"]) / float(run["genomes"]) / 1e6
    return None


def performance_targets_text(
    datasets_meta: Sequence[Dict[str, Any]],
    by_dataset: Dict[str, Dict[str, Dict[str, Any]]],
    summary_cells: Sequence[Dict[str, Any]] = (),
) -> str:
    """Say what the datasets are - simulated or real, how big - next to the targets."""
    simulated: List[str] = []
    real: List[str] = []
    sim_sizes: List[float] = []
    for d in datasets_meta:
        name = str(d.get("name"))
        if name not in by_dataset:
            continue
        mb = _mean_genome_mb(by_dataset.get(name, {}), summary_cells, name)
        count = d.get("genome_count")
        size = (f", ~{mb:.1f} Mb each" if mb >= 0.1 else f", ~{mb * 1000:.0f} kb each") if mb else ""
        desc = f"`{name}` ({common.fmt(count)} genomes{size})"
        if d.get("kind") == "simulated":
            simulated.append(desc)
            if mb:
                sim_sizes.append(mb)
        else:
            real.append(desc)
    parts = [
        "Targets are reported against, not gated. The reference machine is 8 cores "
        "with an SSD and the reference data set is 1,000 simulated **5 Mb** genomes "
        "in 20 clusters."
    ]
    if simulated:
        size_note = ""
        if sim_sizes and max(sim_sizes) < 4.5:
            size_note = (
                f" Their genomes are ~{min(sim_sizes):.1f}"
                + (f"-{max(sim_sizes):.1f}" if max(sim_sizes) - min(sim_sizes) >= 0.1 else "")
                + " Mb, not 5 Mb, so per-genome work is several times lower than in "
                "the target data set."
            )
        parts.append(
            "The simulated datasets here - " + ", ".join(simulated) + " - are "
            "**simulated** genomes from `infinity simulate`, scaling the genome *count* "
            "as section 11 asks." + size_note
        )
    if real:
        parts.append(
            "The real datasets - " + ", ".join(real) + " - are real assemblies but far "
            "fewer than 1,000 genomes."
        )
    parts.append(
        "None of them is the reference data set, so the numbers below are shown "
        "alongside the targets for orientation only; they neither meet nor miss them."
    )
    return " ".join(parts)


def render_report(
    summary: Dict[str, Any],
    by_dataset: Dict[str, Dict[str, Dict[str, Any]]],
    verdicts: Sequence[common.Verdict],
    expected_entries: Sequence[Dict[str, Any]],
    invariance_doc: Dict[str, Any],
    accuracy_doc: Dict[str, Any],
    baselines_doc: Dict[str, Any],
    datasets_meta: Sequence[Dict[str, Any]],
    args: argparse.Namespace,
) -> str:
    root = common.repo_root()
    out: List[str] = []
    cells = summary.get("cells", [])
    totals = summary.get("totals", {})
    n_ok = sum(1 for c in cells if c.get("status") in ("ok", "cached"))
    n_failed = sum(1 for c in cells if c.get("status") == "failed")
    infinity = summary.get("infinity") or {}

    out.append("# Infinity-Win benchmark report")
    out.append("")
    out.append(
        "Does each part of the idea pay? This report is generated by "
        "`scripts/make_report.py` from the configuration x dataset matrix of "
        "PROJECT.md section 11. **Every** result is reported, including the "
        "configurations that failed and the switches that do not pay: section 2 "
        'lists "a benchmark that only reports the cases where the idea wins" as '
        "an anti-goal."
    )
    out.append("")
    out.append(f"- Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z').strip()}")
    out.append(f"- Infinity-Win binary: `{infinity.get('path') or 'not resolved'}`"
               + (f" ({infinity.get('version')})" if infinity.get("version") else ""))
    out.append(f"- Host: {platform.platform()}, Python {platform.python_version()}, "
               f"{os.cpu_count()} CPU(s)")
    out.append(f"- Matrix summary: `{rel(summary.get('_path'), root)}`")
    out.append(f"- Cells: {len(cells)} total - {n_ok} usable, **{n_failed} failed**, "
               f"{len(cells) - n_ok - n_failed} other")
    out.append(f"- Configurations: {', '.join('`' + c['name'] + '`' for c in summary.get('configs', []))}")
    out.append("")
    if summary.get("dry_run"):
        out.append("> **This summary came from a `--dry-run`: nothing was executed.**")
        out.append("")

    # ---- verdict up front ------------------------------------------------
    out.append("## 1. Verdict: does each part of the idea pay?")
    out.append("")
    out.append(
        f"Rule (PROJECT.md section 11): a part of the idea pays if turning it off "
        f"makes **time-to-exact** or **peak memory** at least "
        f"{int(round((common.VERDICT_THRESHOLD - 1) * 100))}% worse on at least "
        f"{common.VERDICT_MIN_DATASETS} of the three datasets."
    )
    out.append("")
    repeats_requested = (summary.get("request") or {}).get("repeats")
    run_counts = sorted(
        {m.get("repeats_usable") for cells_for in by_dataset.values() for m in cells_for.values()
         if m.get("repeats_usable")}
    )
    counts_text = (
        f"{run_counts[0]}" if len(run_counts) == 1
        else (f"{run_counts[0]}-{run_counts[-1]}" if run_counts else "no")
    )
    excluded = next((v.datasets_excluded for v in verdicts if v.datasets_excluded), {})
    min_s = next((v.min_baseline_s for v in verdicts), 0.0)
    out.append("How the rule is applied to noisy measurements:")
    out.append("")
    out.append(
        f"- **Medians of repeated runs.** Every configuration was run "
        f"{repeats_requested or counts_text} time(s) (usable runs per cell: {counts_text}); "
        "the rule compares the **median** of the switched-off runs with the median of "
        "the baseline runs, and every table shows min-max next to the median. "
        + ("With a single run per cell there is no spread to show, so a verdict near "
           "the 20% line is weak evidence." if run_counts == [1] else "")
    )
    if min_s:
        out.append(
            f"- **Too short to judge.** Where the baseline median time-to-exact is under "
            f"{min_s:g} s, the time metric is shown but marked *too short to judge* and "
            "counts neither way: at that scale a 20% difference is process start-up and "
            "page-cache noise. Peak memory on the same dataset is still judged."
        )
    for name, why in sorted(excluded.items()):
        out.append(
            f"- **`{name}` is not part of the performance verdict** ({why}). Its "
            "timings are still reported in section 3, and it is the acceptance "
            "check of section 5."
        )
    out.append("")
    summary_table, verdict_detail = render_verdict(verdicts)
    out.append(summary_table)
    out.append("")
    n_pays = sum(1 for v in verdicts if v.state == "pays")
    n_not = sum(1 for v in verdicts if v.state == "does_not_pay")
    n_undetermined = sum(1 for v in verdicts if v.state == "undetermined")
    n_unmeasurable = sum(1 for v in verdicts if v.state == "not_measurable")
    out.append(
        f"**{n_pays} of {len(verdicts)} parts of the idea pay** by this rule; "
        f"{n_not} do not; {n_undetermined} are undetermined; {n_unmeasurable} could "
        "not be measured at all."
    )
    if n_undetermined or n_unmeasurable:
        out.append("")
        out.append(
            '> Neither "undetermined" nor "not measurable" is a pass. "Not measurable" '
            "means the metrics the rule needs (time-to-exact and peak memory) were "
            "missing for every dataset; \"undetermined\" means fewer than "
            f"{common.VERDICT_MIN_DATASETS} datasets had usable numbers, so the "
            '"at least two of the three datasets" part of the rule cannot be decided '
            "either way. Both happen while the CLI is still being built."
        )
    out.append("")
    out.append("## 2. Datasets")
    out.append("")
    rows = []
    for d in datasets_meta:
        truth = d.get("truth") or {}
        truth_bits = []
        if truth.get("dir"):
            truth_bits.append("truth folder")
        if truth.get("tree"):
            truth_bits.append("reference tree")
        if truth.get("distances"):
            truth_bits.append(f"{(truth.get('distances') or {}).get('format')} distances")
        if truth.get("expected_values"):
            truth_bits.append("published expected values")
        cells_for = by_dataset.get(str(d.get("name")), {})
        rows.append(
            [
                f"`{d.get('name')}`",
                d.get("kind"),
                d.get("organism"),
                common.fmt(d.get("genome_count")),
                common.human_bytes(d.get("total_bytes")),
                ", ".join(truth_bits) or "**none**",
                f"{sum(1 for c in cells_for.values() if c.get('status') in ('ok', 'cached'))}/"
                f"{len(cells_for)} configs usable",
            ]
        )
    out.append(
        common.markdown_table(
            ["dataset", "kind", "organism", "genomes", "size", "truth available", "matrix"],
            rows,
        )
    )
    out.append("")
    for d in datasets_meta:
        caveats = list(d.get("caveats") or [])
        problems = list(d.get("problems") or [])
        if caveats or problems:
            out.append(f"- `{d.get('name')}`")
            for c in caveats:
                out.append(f"  - **Caveat:** {c}")
            for p in problems:
                out.append(f"  - Note: {p}")
    out.append("")

    # ---- per dataset tables ---------------------------------------------
    out.append("## 3. Results per dataset")
    out.append("")
    out.append(
        "`wall` is the measured process time of `infinity run`. `time to first tree` "
        "and `time to exact tree` come from `run.json` timings when present, and "
        "otherwise from the first tree event in `events.jsonl` / from the process "
        "wall time as a labelled proxy. `peak memory` is `run.json "
        "peak_memory_bytes` when present; `(floor)` means it was sampled from "
        "`/proc` while the run was in flight and is therefore a lower bound. "
        "`n/a` means the tool did not report the value - it is never zero."
    )
    out.append("")
    for d in datasets_meta:
        name = str(d.get("name"))
        cells_for = by_dataset.get(name, {})
        out.append(f"### {name}")
        out.append("")
        if not cells_for:
            out.append("_No cells for this dataset in the matrix summary._")
            out.append("")
            continue
        out.append(render_dataset_table(cells_for))
        out.append("")
        out.append("**Accuracy against truth** (from `infinity bench`; `skipped` means the "
                   "dataset has no simulator truth folder, so the engine could not "
                   "score it and `scripts/baselines.py compare` is the comparison to "
                   "read instead):")
        out.append("")
        out.append(render_accuracy_table(cells_for))
        out.append("")
        worst = [
            (config, m.get("worst_pair"))
            for config, m in sorted(cells_for.items())
            if m.get("worst_pair")
        ]
        if worst:
            out.append("Worst-scoring pair per configuration (true SNPs vs reported):")
            out.append("")
            for config, wp in worst:
                out.append(f"- `{config}`: {wp[0]} vs {wp[1]} - true {wp[2]}, reported {wp[3]}")
            out.append("")
        focus_rows = []
        for config, m in sorted(cells_for.items()):
            fc = m.get("focus_counts") or []
            focus_rows.append(
                [
                    config,
                    common.fmt(m.get("genomes")),
                    common.fmt(m.get("clusters")),
                    common.fmt(fc[0] if len(fc) > 0 else None),
                    common.fmt(fc[1] if len(fc) > 1 else None),
                    common.fmt(fc[2] if len(fc) > 2 else None),
                    common.fmt(m.get("bench_status")),
                    (m.get("bench_reason") or "")[:120],
                ]
            )
        out.append(
            common.markdown_table(
                ["config", "genomes", "clusters", "focus 0", "focus 1", "focus 2",
                 "bench", "why bench did not score"],
                focus_rows,
            )
        )
        out.append("")

    # ---- verdict detail --------------------------------------------------
    out.append("## 4. Verdict detail, per switch and per dataset")
    out.append("")
    out.append(verdict_detail)

    # ---- acceptance check ------------------------------------------------
    out.append("## 5. Acceptance check against published values")
    out.append("")
    out.append(
        "These numbers were produced by analyses that predate Infinity-Win (SKA2 on "
        "reads; reference mapping after filtering) and are an acceptance check, "
        "not a tuning target. Provenance for every pair is in "
        "`scripts/expected/`."
    )
    out.append("")
    out.append(render_expected(expected_entries))
    out.append("")

    # ---- switch invariance ----------------------------------------------
    out.append("## 6. Switch invariance: do the switches change the answer?")
    out.append("")
    out.append(
        "PROJECT.md section 10.2 property 6 requires that flat mode, hierarchical "
        "mode and every cache switch give **identical** variants, distances and "
        "trees: each switch removes work, not accuracy. This compares every "
        "configuration's pairwise SNP counts *and* its variant lists with the "
        "baseline's. A switch that is faster but wrong has not been shown to pay, "
        "so a difference is classified rather than quietly averaged away:"
    )
    out.append("")
    out.append(
        "| classification | meaning |\n|---|---|\n"
        "| identical | reproduced the baseline exactly |\n"
        f"| **DEFECT** | differs, and nothing in DECISIONS.md says it may "
        f"({', '.join('`' + m + '`' for m in MUST_BE_EXACTLY_NEUTRAL)}) |\n"
        "| documented deviation (D16) | `no_block_skip` changed distances while its "
        "variant list stayed identical - exactly what DECISIONS.md D16 records |\n"
        "| **UNDOCUMENTED** | `no_block_skip` changed the *variants* too, which D16 "
        "says must stay identical |\n"
        "| expected (different pivot) | `single_cluster` compares every genome to one "
        "pivot, so its pair set and coordinates differ by construction |\n"
    )
    counts: Dict[str, int] = {}
    for entry in invariance_doc.values():
        for rec in (entry.get("configs") or {}).values():
            status = str(rec.get("status"))
            counts[status] = counts.get(status, 0) + 1
    n_defect = counts.get("DEFECT", 0) + counts.get("UNDOCUMENTED", 0)
    n_documented = counts.get("documented deviation", 0)
    n_identical = counts.get("identical", 0)
    if n_defect:
        out.append(
            f"> **{n_defect} configuration/dataset combination(s) differ from the "
            "baseline with no documented reason.** Section 10.2 property 6 is not "
            "satisfied, so the timing differences in section 3 cannot be read as "
            "pure optimisation gains until this is fixed."
        )
        out.append("")
    if n_documented:
        out.append(
            f"> `no_block_skip` differs on {n_documented} dataset(s) with its variant "
            "list unchanged. That is the coverage-bookkeeping difference recorded in "
            "DECISIONS.md D16, not a defect - but per D16 it is reported here "
            "alongside its speed difference, and it means `--no-block-skip` is **not** "
            "a pure performance switch: block skipping costs a measurable amount of "
            "coverage on fragmented assemblies."
        )
        out.append("")
    if not n_defect and not n_documented:
        out.append(
            f"> Every optimisation switch reproduced the baseline exactly "
            f"({n_identical} combination(s)): section 10.2 property 6 holds."
        )
        out.append("")
    out.append(render_invariance(invariance_doc))
    out.append("")

    # ---- accuracy over time ---------------------------------------------
    out.append("## 7. The key figure: tree accuracy versus elapsed time")
    out.append("")
    out.append(
        "Accuracy at a given instant is the normalised Robinson-Foulds distance "
        "between the **whole tree the program was showing** and the truth tree, "
        "over the **full leaf set** - every genome of the dataset that the truth "
        "knows - after collapsing branches shorter than half a SNP. The curves are "
        "built by replaying `events.jsonl`: after every `backbone_tree` / "
        "`cluster_focus` event the displayed tree is re-assembled as PROJECT.md "
        "section 6.9 describes (the backbone, with each pivot leaf replaced by its "
        "cluster), where a cluster that has not emitted a tree yet is an "
        "unresolved star of its members and every other cluster is its latest "
        "emitted tree. So each curve starts near 1 (a star resolves nothing) and "
        "falls as clusters become exact."
    )
    out.append("")
    out.append(
        "- **Normalisation:** RF / (splits in the displayed tree + splits in the "
        "truth): 1 = no split right, 0 = identical topology. Splits are "
        "unrooted bipartitions of the full leaf set."
    )
    out.append(
        "- **Simulated truth:** `infinity simulate` writes one truth tree per "
        "cluster and no full tree. Every cluster root is derived from the "
        "reference independently, so the full truth is those trees joined at one "
        "star: each cluster is a clade, and no order of clusters is claimed."
    )
    out.append(
        "- **Why the simulated curves do not end at exactly 0:** the full tree also "
        "claims where each cluster hangs off the rest. Two things the per-cluster "
        "RF of `infinity bench` ignores are counted here: splits *between* "
        "clusters that the backbone resolves although the truth is a star, and "
        "splits *within* a cluster that move because the cluster tree is rooted "
        "at its midpoint (section 6.9) rather than at the true cluster root. The "
        "`final RF` column below splits the residual into those two parts "
        "(between / within)."
    )
    out.append(
        "- **Real data:** scored against the dataset's reference tree, restricted "
        "to the genomes that are both in the run and in that tree; extra "
        "reference/outgroup leaves are left out."
    )
    out.append("")
    figures: List[str] = []
    if accuracy_doc.get("combined_figure"):
        figures.append(accuracy_doc["combined_figure"])
    for dataset, info in sorted((accuracy_doc.get("datasets") or {}).items()):
        if info.get("figure"):
            figures.append(info["figure"])
    if figures:
        for f in figures:
            out.append(f"![accuracy over time]({rel(f, root)})")
            out.append("")
        out.append(
            f"_Figures live under `{rel(accuracy_doc.get('figures_dir'), root)}/`, "
            "which is gitignored benchmark scratch space; re-run "
            "`scripts/make_report.py` to regenerate them._"
        )
    else:
        out.append(
            "**No figure could be drawn.** "
            + "; ".join(accuracy_doc.get("problems") or [])
        )
    out.append("")
    csv_path = accuracy_doc.get("combined_csv")
    if csv_path:
        out.append(
            f"Every number behind the curves is in "
            f"`{rel(csv_path, root)}` ({common.fmt(accuracy_doc.get('combined_rows'))} rows), "
            "so nothing is lost when plotting is unavailable."
        )
        out.append("")
    endpoint_rows = curve_endpoints(accuracy_doc)
    if endpoint_rows:
        out.append(
            "Curve endpoints (the same data as the figure, as text). `first` is the "
            "first instant whose RF is defined; RF is undefined (never zero) when "
            "fewer than four genomes can be scored - those points are counted in "
            "`no RF`. `collapse threshold` is the truth tree's, in its own units."
        )
        out.append("")
        out.append(
            common.markdown_table(
                ["dataset", "config", "status", "points", "no RF", "RF source",
                 "collapse threshold", "first t (s)", "first RF (norm)",
                 "last t (s)", "last RF (norm)", "final RF (between / within)",
                 "why empty (if it is)"],
                endpoint_rows,
            )
        )
        out.append("")
    for problem in accuracy_doc.get("problems") or []:
        out.append(f"> {problem}")
    out.append("")

    # ---- baselines -------------------------------------------------------
    out.append("## 8. Outside baselines")
    out.append("")
    out.append(
        "Linux only, not part of the program (PROJECT.md section 11). Neither "
        "Mash nor SKA2 is installed on the machine this was generated on, so both "
        "degrade to pre-computed results; that is recorded per row below rather "
        "than hidden."
    )
    out.append("")
    out.append(render_baselines(baselines_doc, root))
    out.append("")

    # ---- failures --------------------------------------------------------
    out.append("## 9. Failures, gaps and data-quality notes")
    out.append("")
    out.append(
        "Everything that did not work, or that worked only partially, is listed "
        "here. Nothing in this section is omitted from the tables above."
    )
    out.append("")
    failure_rows = []
    for cell in cells:
        problems = list(cell.get("problems") or [])
        if cell.get("status") == "failed":
            detail = f"FAILED: {cell.get('error') or 'unknown error'}"
            # The child's own stderr is usually the actual cause, so it belongs in
            # the report and not only in the per-cell log file.
            tail = (cell.get("stderr_tail") or "").strip()
            if tail:
                detail += f" | stderr: {tail.splitlines()[-1][:200]}"
            problems.insert(0, detail)
        bench = cell.get("bench") or {}
        if bench.get("status") == "failed":
            problems.append(f"bench failed: {bench.get('reason')}")
        elif bench.get("status") == "skipped":
            problems.append(f"bench skipped: {bench.get('reason')}")
        for p in bench.get("problems") or []:
            problems.append(f"bench: {p}")
        run = cell.get("run") or {}
        if run.get("switch_verified") is False:
            problems.append("the switch did not take effect according to run.json params")
        if problems:
            failure_rows.append(
                [
                    cell.get("dataset"),
                    cell.get("config") or "-",
                    cell.get("status"),
                    "; ".join(dict.fromkeys(problems))[:600],
                ]
            )
    if failure_rows:
        out.append(common.markdown_table(["dataset", "config", "status", "what went wrong / what is missing"], failure_rows))
    else:
        out.append("_No failures and no missing values were recorded._")
    out.append("")
    if summary.get("failures"):
        out.append("Failure list recorded by `run_matrix.py`:")
        out.append("")
        for f in summary["failures"]:
            out.append(f"- {f}")
        out.append("")

    # ---- performance targets --------------------------------------------
    out.append("## 10. Performance targets (PROJECT.md section 10.4)")
    out.append("")
    out.append(performance_targets_text(datasets_meta, by_dataset, cells))
    out.append("")
    target_rows = []
    for name, cells_for in sorted(by_dataset.items()):
        base = cells_for.get("baseline") or {}
        target_rows.append(
            [
                name,
                common.fmt(base.get("genomes")),
                fmt_stat(base, "first_tree_s"),
                "<= 60 s",
                fmt_stat(base, "exact_tree_s"),
                "<= 5 min",
                fmt_stat(base, "peak_memory_bytes", memory=True),
                "<= 4 GB",
            ]
        )
    out.append(
        common.markdown_table(
            ["dataset (baseline)", "genomes", "first tree", "target", "exact tree",
             "target", "peak memory", "target"],
            target_rows,
        )
    )
    out.append("")

    # ---- reproduction ---------------------------------------------------
    out.append("## 11. How to reproduce")
    out.append("")
    out.append("```bash")
    out.append("export PATH=\"$HOME/.cargo/bin:$PATH\"")
    out.append("cargo build --release                     # binaries land in /tmp/infinity-target/release")
    out.append("")
    out.append("scripts/real_datasets.py --datasets all   # copy + verify the real datasets")
    out.append("scripts/run_matrix.py --dry-run           # show every command line first")
    out.append("scripts/run_matrix.py                     # 6 configurations x every dataset")
    out.append("scripts/baselines.py all --dataset ncbi_pd_cluster \\")
    out.append("    --run-dir bench-work/matrix/ncbi_pd_cluster/baseline")
    out.append("scripts/make_report.py                    # writes REPORT.md")
    out.append("```")
    out.append("")
    out.append(
        f"The matrix command lines actually used are stored per cell in "
        f"`{rel(summary.get('_path'), root)}` under `cells[].run_command` and "
        "`cells[].bench_command`, so any single cell can be re-run by hand."
    )
    out.append("")
    out.append("---")
    out.append("")
    out.append(
        f"_Generated by `scripts/make_report.py` (report schema v{REPORT_VERSION}) on "
        f"{time.strftime('%Y-%m-%d %H:%M:%S %Z').strip()}. Numbers marked `n/a` were "
        "not reported by the tool and were never substituted with zero._"
    )
    out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="make_report.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--work-dir", type=Path, default=None, help="bench-work directory")
    p.add_argument("--matrix-dir", type=Path, default=None, help="matrix directory (default <work>/matrix)")
    p.add_argument("--summary", type=Path, default=None, help="summary.json (default <matrix>/summary.json)")
    p.add_argument("--out", type=Path, default=None, help="report path (default <repo>/REPORT.md)")
    p.add_argument("--figures-dir", type=Path, default=None, help="where figures/CSVs go")
    p.add_argument("--baselines-dir", type=Path, default=None, help="baseline outputs (default <work>/baselines)")
    p.add_argument(
        "--expected-dir",
        type=Path,
        default=None,
        help="folder of expected-values JSON files for the acceptance check "
        "(default <repo>/scripts/expected)",
    )
    p.add_argument("--no-figures", action="store_true", help="do not regenerate accuracy curves")
    p.add_argument("--collapse-snps", type=float, default=0.5, help="RF collapse threshold in SNPs")
    p.add_argument(
        "--min-judge-s",
        type=float,
        default=common.VERDICT_MIN_BASELINE_S,
        help="time-to-exact is only judged when the baseline median is at least this "
        f"many seconds (default {common.VERDICT_MIN_BASELINE_S:g}; 0 judges everything)",
    )
    p.add_argument("--max-points", type=int, default=400, help="cap scored points per configuration")
    p.add_argument("--stdout", action="store_true", help="also print the report to stdout")
    p.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    matrix_dir = Path(args.matrix_dir) if args.matrix_dir else work / "matrix"
    summary_path = Path(args.summary) if args.summary else matrix_dir / "summary.json"
    out_path = Path(args.out) if args.out else common.repo_root() / "REPORT.md"
    figures_dir = Path(args.figures_dir) if args.figures_dir else matrix_dir / "figures"

    if not summary_path.exists():
        raise common.PrereqError(
            f"no matrix summary at {summary_path}\n"
            "       run scripts/run_matrix.py first (it writes summary.json after every "
            "cell, so a partially finished matrix still produces a report)"
        )
    summary = common.load_json(summary_path)
    if not isinstance(summary, dict):
        raise common.PrereqError(f"{summary_path} is not a JSON object")
    summary["_path"] = str(summary_path)
    if not summary.get("cells"):
        print(
            f"WARNING: {summary_path} has no cells; the report will document an "
            "empty/failed matrix rather than pretend there are results",
            file=sys.stderr,
        )

    by_dataset = collect(summary, summary_path.parent)

    # dataset metadata, straight from the dataset.json files
    datasets_meta: List[Dict[str, Any]] = []
    for entry in summary.get("datasets", []):
        name = entry.get("name")
        path = work / str(name) / "dataset.json"
        obj, _ = common.load_json_optional(path)
        merged = dict(entry)
        if isinstance(obj, dict):
            for key in ("kind", "organism", "genome_count", "total_bytes", "truth",
                        "caveats", "problems", "licence", "description"):
                if obj.get(key) is not None:
                    merged[key] = obj[key]
        datasets_meta.append(merged)
    for name in sorted(by_dataset):
        if not any(d.get("name") == name for d in datasets_meta):
            datasets_meta.append({"name": name, "kind": None, "genome_count": None})

    # verdict per switch. The rule counts "at least two of the three datasets",
    # so the denominator is every dataset the matrix knows about, including the
    # ones that produced no usable cell at all.
    verdicts: List[common.Verdict] = []
    excluded = common.performance_exclusions(datasets_meta)
    dataset_names = sorted({str(d.get("name")) for d in datasets_meta if d.get("name")} | set(by_dataset))
    for config in common.CONFIGS:
        if config.is_baseline:
            continue
        baseline_by = {
            d: by_dataset[d].get("baseline", {}) for d in dataset_names if "baseline" in by_dataset.get(d, {})
        }
        off_by = {
            d: by_dataset[d][config.name]
            for d in dataset_names
            if config.name in by_dataset.get(d, {})
        }
        verdicts.append(
            common.evaluate_verdict(
                config.name,
                config.removes or "",
                baseline_by,
                off_by,
                dataset_names=dataset_names,
                min_baseline_s=args.min_judge_s,
                excluded=excluded,
            )
        )

    expected_entries, expected_problems = run_expected_checks(
        work, summary, by_dataset, expected_dir=args.expected_dir
    )
    accuracy_doc = build_accuracy_curves(
        by_dataset,
        work,
        figures_dir,
        collapse_snps=args.collapse_snps,
        max_points=args.max_points,
        enabled=not args.no_figures,
        summary=summary,
    )
    for problem in expected_problems:
        accuracy_doc.setdefault("problems", []).append(problem)
    baselines_doc = collect_baselines(work)
    invariance_doc = check_switch_invariance(by_dataset)

    text = render_report(
        summary,
        by_dataset,
        verdicts,
        expected_entries,
        invariance_doc,
        accuracy_doc,
        baselines_doc,
        datasets_meta,
        args,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")
    print(f"report written to {out_path} ({len(text.splitlines())} lines)")
    if args.stdout:
        print()
        print(text)
    n_failed = sum(1 for c in summary.get("cells", []) if c.get("status") == "failed")
    if n_failed:
        print(
            f"NOTE: {n_failed} matrix cell(s) FAILED and are reported as failed in "
            "section 8 of the report; they are not omitted",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    common.guard(main)()
