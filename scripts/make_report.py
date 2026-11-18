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
    """The numbers the report and the verdict rule use for one cell."""
    run = cell.get("run") or {}
    timings_s = run.get("timings_s") or {}
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
    }
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
        out["exact_tree_source"] = (
            "process wall time (proxy: run.json carried no timings.exact_tree_ms)"
        )
        out["problems"].append(
            "time-to-exact is the whole `mtato run` wall time, used because "
            "run.json has no timings.exact_tree_ms; it includes output writing"
        )
    bench = out["bench"]
    out["rf"] = bench.get("rf")
    out["rf_normalized"] = bench.get("rf_normalized")
    out["distance_mae_snps"] = bench.get("distance_mae_snps")
    out["distance_max_error_snps"] = bench.get("distance_max_error_snps")
    out["pairs_compared"] = bench.get("pairs_compared")
    out["pairs_exact"] = bench.get("pairs_exact")
    out["bench_status"] = bench.get("status")
    out["bench_reason"] = bench.get("reason")
    return out


def collect(summary: Dict[str, Any]) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """``{dataset: {config: measurements}}``."""
    out: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for cell in summary.get("cells", []):
        dataset = cell.get("dataset")
        config = cell.get("config")
        if not dataset or not config:
            continue
        out.setdefault(str(dataset), {})[str(config)] = measurements_for_cell(cell)
    return out


# ---------------------------------------------------------------------------
# Acceptance checks against published expected values
# ---------------------------------------------------------------------------


def run_expected_checks(
    work: Path, summary: Dict[str, Any], by_dataset: Dict[str, Dict[str, Dict[str, Any]]]
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Compare M-Tato's pairwise SNP counts with published expected values."""
    results: List[Dict[str, Any]] = []
    problems: List[str] = []
    expected_dir = common.repo_root() / "scripts" / "expected"
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
            cell_problems: List[str] = []
            if out_dir and Path(out_dir).is_dir():
                rows, row_problems = common.load_cluster_distances(Path(out_dir))
                cell_problems.extend(row_problems)
                observed = {
                    common.pair_key(r["a"], r["b"]): r.get("snps") for r in rows
                }
                if not observed:
                    cell_problems.append(
                        "no pairwise distances could be read from this output folder"
                    )
            else:
                cell_problems.append(f"output folder missing: {out_dir}")
            checks = common.check_expected_pairs(expected, observed)
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
) -> Dict[str, Any]:
    """Drive ``accuracy_over_time.py`` over every usable cell."""
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
    for dataset, cells in sorted(by_dataset.items()):
        dataset_dir = work / dataset
        truth_tree: Optional[Path] = None
        truth_dir: Optional[Path] = None
        if dataset_dir.is_dir():
            obj, _ = common.load_json_optional(dataset_dir / "dataset.json")
            if isinstance(obj, dict):
                truth = obj.get("truth") or {}
                tree = truth.get("tree")
                if tree and Path(tree).exists():
                    truth_tree = Path(tree)
                tdir = truth.get("dir")
                if tdir and Path(tdir).is_dir():
                    truth_dir = Path(tdir)
        per_dataset: Dict[str, Any] = {
            "runs": {},
            "truth_tree": str(truth_tree) if truth_tree else None,
            "truth_dir": str(truth_dir) if truth_dir else None,
            "simulator_truth": bool(truth_dir and common.is_simulator_truth_dir(truth_dir)),
        }
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
                    mtato=None,
                )
            except common.PrereqError as exc:
                res = {"status": "failed", "problems": [str(exc)], "points": 0, "rows": []}
            rows = res.pop("rows", [])
            res.pop("dataset", None)
            res.pop("config", None)
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
                        }
                    res["t_span_s"] = [
                        float(ordered[0]["t_ms"]) / 1000.0,
                        float(ordered[-1]["t_ms"]) / 1000.0,
                    ]
            per_dataset["runs"][config] = res
            if rows:
                series[config] = rows
                all_series[f"{dataset}/{config}"] = rows
        if series:
            png = figures_dir / f"accuracy_over_time__{dataset}.png"
            written, reason = aot.plot_curves(
                series, png, title=f"Tree accuracy over time - {dataset}"
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
        written, reason = aot.plot_curves(all_series, png, title="Tree accuracy over time - all datasets")
        doc["combined_figure"] = str(png) if written else None
        if not written:
            doc["problems"].append(f"combined figure not written: {reason}")
    return doc


def curve_endpoints(doc: Dict[str, Any]) -> List[List[Any]]:
    """A text stand-in for the figure: first and last defined point of every curve."""
    rows: List[List[Any]] = []
    for dataset, info in sorted((doc.get("datasets") or {}).items()):
        for config, run in sorted((info.get("runs") or {}).items()):
            points = run.get("points") or 0
            first = run.get("first") or {}
            last = run.get("last") or {}
            rows.append(
                [
                    dataset,
                    config,
                    run.get("status"),
                    points,
                    run.get("points_without_rf"),
                    run.get("rf_source"),
                    common.fmt((run.get("collapse") or {}).get("threshold"), ".3g"),
                    common.fmt(first.get("t_s"), ".4g"),
                    common.fmt(first.get("rf_normalized"), ".4f"),
                    common.fmt(last.get("t_s"), ".4g"),
                    common.fmt(last.get("rf_normalized"), ".4f"),
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


def verdict_mark(passed: Optional[bool]) -> str:
    """Bold PASS/FAIL, or an unbolded 'not measured' - never a silent pass."""
    if passed is True:
        return "**PASS**"
    if passed is False:
        return "**FAIL**"
    return "not measured"


def render_dataset_table(cells: Dict[str, Dict[str, Any]]) -> str:
    header = [
        "config", "status", "wall", "time to first tree", "time to exact tree",
        "peak memory", "block identity", "pieces compared", "cache hits",
        "alignments", "bases aligned", "breakpoints", "RF", "RF (norm)",
        "SNP MAE", "switch applied",
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
        rows.append(
            [
                config,
                status,
                common.human_seconds(m.get("wall_s")),
                common.human_seconds(m.get("first_tree_s")),
                common.human_seconds(m.get("exact_tree_s")),
                common.human_bytes(m.get("peak_memory_bytes"))
                + (" (floor)" if m.get("peak_memory_is_floor") else ""),
                common.fmt(identity, ".3%") if identity is not None else "n/a",
                common.fmt(counters.get("pieces_compared")),
                common.fmt(counters.get("cache_hits")),
                common.fmt(counters.get("alignments_run")),
                common.fmt(counters.get("bases_aligned")),
                common.fmt(counters.get("breakpoints")),
                common.fmt(m.get("rf"), ".4g"),
                common.fmt(m.get("rf_normalized"), ".4f"),
                common.fmt(m.get("distance_mae_snps"), ".4g"),
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
            if m.metric.startswith("peak memory"):
                base = common.human_bytes(m.baseline)
                off = common.human_bytes(m.switched_off)
            else:
                base = common.human_seconds(m.baseline)
                off = common.human_seconds(m.switched_off)
            worse = {True: "**yes**", False: "no", None: "n/a"}[m.worse]
            rows.append(
                [
                    m.dataset,
                    m.metric,
                    base,
                    off,
                    common.fmt(m.ratio, ".3f") + "x" if m.ratio is not None else "n/a",
                    common.fmt(common.pct_change(m.switched_off, m.baseline), "+.1f") + "%"
                    if m.ratio is not None
                    else "n/a",
                    worse,
                    m.note,
                ]
            )
        detail.append(
            common.markdown_table(
                ["dataset", "metric", "baseline", "switch off", "ratio", "change",
                 f">= {int(round((v.threshold - 1) * 100))}% worse?"],
                [r[:7] for r in rows],
            )
        )
        notes = [r[7] for r in rows if r[7]]
        if notes:
            detail.append("")
            detail.append("Notes: " + "; ".join(sorted(set(notes))))
        detail.append("")
        detail.append(f"**Verdict: {label.replace('**', '')}.** {v.reason}")
        detail.append("")
    summary = common.markdown_table(
        ["switch", "what turning it off removes", "verdict", "datasets >=20% worse",
         "datasets measurable", "why"],
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
        rows = []
        for config, info in sorted(by_config.items()):
            checks = info.get("checks") or []
            for check in checks:
                pair = "-".join(check.get("pair") or []) or check.get("name")
                rows.append(
                    [
                        config,
                        pair,
                        common.fmt(check.get("expected"), ".4g"),
                        common.fmt(check.get("observed"), ".4g"),
                        (
                            common.fmt(
                                None
                                if check.get("observed") is None
                                or check.get("expected") is None
                                else float(check["observed"]) - float(check["expected"]),
                                "+.4g",
                            )
                        ),
                        f"+/- {common.fmt(check.get('tolerance'), '.4g')}",
                        verdict_mark(check.get("passed")),
                        check.get("source") or "",
                    ]
                )
            status = info.get("status")
            rows.append(
                [f"`{config}` overall", "", "", "", "", "", f"**{status}**",
                 "; ".join(info.get("problems") or [])[:160]]
            )
        out.append("")
        out.append(
            common.markdown_table(
                ["config", "pair", "expected", "observed", "diff", "tolerance", "verdict", "provenance / note"],
                rows,
            )
        )
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
                    f"M-Tato vs SKA2 ({cmp_doc.get('run_config') or key})",
                    cmp_doc.get("status"),
                    "-",
                    f"{common.fmt(summary.get('pairs_in_both'))} pairs in both "
                    f"(only M-Tato: {common.fmt(summary.get('pairs_only_in_mtato'))}, "
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


def render_report(
    summary: Dict[str, Any],
    by_dataset: Dict[str, Dict[str, Dict[str, Any]]],
    verdicts: Sequence[common.Verdict],
    expected_entries: Sequence[Dict[str, Any]],
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
    mtato = summary.get("mtato") or {}

    out.append("# M-Tato benchmark report")
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
    out.append(f"- M-Tato binary: `{mtato.get('path') or 'not resolved'}`"
               + (f" ({mtato.get('version')})" if mtato.get("version") else ""))
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
        "`wall` is the measured process time of `mtato run`. `time to first tree` "
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
        "These numbers were produced by analyses that predate M-Tato (SKA2 on "
        "reads; reference mapping after filtering) and are an acceptance check, "
        "not a tuning target. Provenance for every pair is in "
        "`scripts/expected/`."
    )
    out.append("")
    out.append(render_expected(expected_entries))
    out.append("")

    # ---- accuracy over time ---------------------------------------------
    out.append("## 6. The key figure: tree accuracy versus elapsed time")
    out.append("")
    out.append(
        "Accuracy at a given instant is the Robinson-Foulds distance between the "
        "tree M-Tato had emitted by then and the truth tree, after collapsing "
        "branches shorter than half a SNP. The curves are built by replaying "
        "`events.jsonl` and re-assembling the full tree after every "
        "`backbone_tree` / `cluster_focus` event (PROJECT.md section 6.9)."
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
            "first instant whose RF is defined; while only part of the tree exists "
            "there may be fewer than four shared taxa and RF is then undefined "
            "rather than zero - those points are counted in `no RF`."
        )
        out.append("")
        out.append(
            common.markdown_table(
                ["dataset", "config", "status", "points", "no RF", "RF source",
                 "collapse threshold", "first t (s)", "first RF (norm)",
                 "last t (s)", "last RF (norm)"],
                endpoint_rows,
            )
        )
        out.append("")
    for problem in accuracy_doc.get("problems") or []:
        out.append(f"> {problem}")
    out.append("")

    # ---- baselines -------------------------------------------------------
    out.append("## 7. Outside baselines")
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
    out.append("## 8. Failures, gaps and data-quality notes")
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
    out.append("## 9. Performance targets (PROJECT.md section 10.4)")
    out.append("")
    out.append(
        "Targets are reported against, not gated. The reference machine is 8 "
        "cores with an SSD and the reference data set is 1,000 simulated 5 Mb "
        "genomes in 20 clusters; the datasets in this report are the real ones, "
        "so the numbers are not directly comparable to the targets and are shown "
        "alongside them for orientation only."
    )
    out.append("")
    target_rows = []
    for name, cells_for in sorted(by_dataset.items()):
        base = cells_for.get("baseline") or {}
        target_rows.append(
            [
                name,
                common.fmt(base.get("genomes")),
                common.human_seconds(base.get("first_tree_s")),
                "<= 60 s",
                common.human_seconds(base.get("exact_tree_s")),
                "<= 5 min",
                common.human_bytes(base.get("peak_memory_bytes")),
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
    out.append("## 10. How to reproduce")
    out.append("")
    out.append("```bash")
    out.append("export PATH=\"$HOME/.cargo/bin:$PATH\"")
    out.append("cargo build --release                     # binaries land in /tmp/mtato-target/release")
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
    p.add_argument("--no-figures", action="store_true", help="do not regenerate accuracy curves")
    p.add_argument("--collapse-snps", type=float, default=0.5, help="RF collapse threshold in SNPs")
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

    by_dataset = collect(summary)

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
            )
        )

    expected_entries, expected_problems = run_expected_checks(work, summary, by_dataset)
    accuracy_doc = build_accuracy_curves(
        by_dataset,
        work,
        figures_dir,
        collapse_snps=args.collapse_snps,
        max_points=args.max_points,
        enabled=not args.no_figures,
    )
    for problem in expected_problems:
        accuracy_doc.setdefault("problems", []).append(problem)
    baselines_doc = collect_baselines(work)

    text = render_report(
        summary,
        by_dataset,
        verdicts,
        expected_entries,
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
