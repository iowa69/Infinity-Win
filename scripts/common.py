#!/usr/bin/env python3
"""Shared helpers for the Infinity-Win benchmark and reporting scripts.

Design rules this module exists to enforce:

* **Standard library only.** Every script must run on a bare Python 3.8+
  install. ``matplotlib``/``numpy`` are optional and are imported lazily by the
  plotting code only; when they are missing the numbers are still written as
  CSV/TSV so nothing is lost.
* **Parse defensively.** The Rust CLI is being built concurrently with these
  scripts (PROJECT.md section 13). A column, a key or a whole file may be
  missing while the tool is still incomplete. Readers here never raise on a
  missing field: they return ``None`` and record a human-readable *problem*
  string that the report prints. Section 0 rule 8 ("never present an
  approximate number as exact") means an unknown value stays unknown.
* **Fail loudly on missing prerequisites.** Anything that stops a script from
  doing its job at all raises :class:`PrereqError`, which every entry point
  turns into one actionable line on stderr instead of a traceback.

Nothing in here imports the Rust crates or shells out to ``cargo``.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
)

# ---------------------------------------------------------------------------
# Errors and process plumbing
# ---------------------------------------------------------------------------


class PrereqError(Exception):
    """A prerequisite is missing. The message must say how to fix it."""


def die(message: str, code: int = 2) -> None:
    """Print one actionable line and exit. Never raises a traceback."""
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(code)


def guard(main: Callable[[Optional[Sequence[str]]], int]) -> Callable[[], None]:
    """Wrap a ``main(argv)`` so missing prerequisites print, not traceback.

    Unexpected exceptions are summarised to one line unless ``--debug`` was
    passed or ``INFINITY_SCRIPTS_DEBUG=1`` is set; the full traceback is then
    printed, because a real bug is easier to find with one.
    """

    def runner() -> None:
        argv = sys.argv[1:]
        debug = bool(os.environ.get("INFINITY_SCRIPTS_DEBUG")) or "--debug" in argv
        try:
            code = main(argv)
        except PrereqError as exc:
            die(str(exc))
        except BrokenPipeError:
            raise SystemExit(0)
        except KeyboardInterrupt:
            print("interrupted", file=sys.stderr)
            raise SystemExit(130)
        except SystemExit:
            raise
        except Exception as exc:  # noqa: BLE001 - deliberate catch-all
            if debug:
                raise
            print(
                f"ERROR: unexpected {type(exc).__name__}: {exc}\n"
                "       re-run with --debug (or INFINITY_SCRIPTS_DEBUG=1) for a traceback",
                file=sys.stderr,
            )
            raise SystemExit(1)
        raise SystemExit(code or 0)

    return runner


@dataclass
class ChildResult:
    """Outcome of one subprocess we benchmarked."""

    cmd: List[str]
    returncode: Optional[int]
    wall_s: float
    timed_out: bool = False
    #: Peak RSS sampled from ``/proc/<pid>/status`` (``VmHWM``), in bytes.
    #: ``VmHWM`` is a high-water mark, so any sample is a true *lower bound* on
    #: the real peak. ``None`` when the platform does not expose it.
    peak_rss_bytes_polled: Optional[int] = None
    memory_is_floor: bool = True
    stdout_path: Optional[Path] = None
    stderr_path: Optional[Path] = None

    @property
    def ok(self) -> bool:
        return self.returncode == 0 and not self.timed_out

    def stderr_tail(self, limit: int = 1200) -> str:
        if not self.stderr_path or not self.stderr_path.exists():
            return ""
        try:
            text = self.stderr_path.read_text(errors="replace")
        except OSError:
            return ""
        text = text.strip()
        return text[-limit:] if len(text) > limit else text


def _poll_vmhwm(pid: int, stop: threading.Event, out: List[Optional[int]]) -> None:
    """Sample ``VmHWM`` of ``pid`` until ``stop`` is set. Linux only."""
    path = f"/proc/{pid}/status"
    best: Optional[int] = None
    while not stop.is_set():
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith("VmHWM:"):
                        kb = int(line.split()[1])
                        best = kb * 1024 if best is None else max(best, kb * 1024)
                        break
        except (OSError, ValueError, IndexError):
            pass  # process gone, or a partial read: keep the best we have
        stop.wait(0.05)
    out[0] = best


def run_child(
    cmd: Sequence[str],
    *,
    stdout_path: Optional[Path] = None,
    stderr_path: Optional[Path] = None,
    timeout: Optional[float] = None,
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    poll_memory: bool = True,
    capture: bool = False,
) -> ChildResult:
    """Run ``cmd`` and time it.

    ``capture=True`` keeps the output in memory instead of a file (used for
    short helper calls such as ``infinity --version``). Long benchmark runs should
    pass file paths so the log survives a crash.
    """
    cmd = [str(c) for c in cmd]
    if capture:
        stdout_path = stderr_path = None
    started = time.monotonic()
    out_fh = err_fh = None
    proc: Optional[subprocess.Popen] = None
    peak: List[Optional[int]] = [None]
    stop = threading.Event()
    poller: Optional[threading.Thread] = None
    captured_out = captured_err = ""
    try:
        if capture:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=str(cwd) if cwd else None,
                env=env,
                text=True,
            )
        else:
            if stdout_path is not None:
                stdout_path.parent.mkdir(parents=True, exist_ok=True)
                out_fh = open(stdout_path, "wb")
            if stderr_path is not None:
                stderr_path.parent.mkdir(parents=True, exist_ok=True)
                err_fh = open(stderr_path, "wb")
            proc = subprocess.Popen(
                cmd,
                stdout=out_fh if out_fh else subprocess.DEVNULL,
                stderr=err_fh if err_fh else (out_fh or subprocess.DEVNULL),
                cwd=str(cwd) if cwd else None,
                env=env,
            )
        if poll_memory and proc.pid and Path(f"/proc/{proc.pid}/status").exists():
            poller = threading.Thread(
                target=_poll_vmhwm, args=(proc.pid, stop, peak), daemon=True
            )
            poller.start()
        timed_out = False
        try:
            if capture:
                captured_out, captured_err = proc.communicate(timeout=timeout)
            else:
                proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            proc.kill()
            if capture:
                captured_out, captured_err = proc.communicate()
            else:
                proc.wait()
        rc = proc.returncode
    finally:
        stop.set()
        if poller is not None:
            poller.join(timeout=1.0)
        for fh in (out_fh, err_fh):
            if fh is not None:
                try:
                    fh.close()
                except OSError:
                    pass
    elapsed = time.monotonic() - started
    result = ChildResult(
        cmd=cmd,
        returncode=rc,
        wall_s=elapsed,
        timed_out=timed_out,
        peak_rss_bytes_polled=peak[0],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    if capture:
        # Stash captured text where callers can reach it without a file.
        result.stdout_path = None
        result.stderr_path = None
        _CAPTURED[id(result)] = (captured_out or "", captured_err or "")
    return result


_CAPTURED: Dict[int, Tuple[str, str]] = {}


def captured(result: ChildResult) -> Tuple[str, str]:
    """Return ``(stdout, stderr)`` for a ``run_child(..., capture=True)``."""
    return _CAPTURED.get(id(result), ("", ""))


# ---------------------------------------------------------------------------
# Paths and binary resolution
# ---------------------------------------------------------------------------


def repo_root(start: Optional[Path] = None) -> Path:
    """The Infinity-Win repository root (the folder holding ``Cargo.toml``)."""
    here = (start or Path(__file__)).resolve()
    for parent in [here, *here.parents]:
        if (parent / "Cargo.toml").exists() and (parent / "PROJECT.md").exists():
            return parent
    return Path(__file__).resolve().parent.parent


def default_work_dir() -> Path:
    """``bench-work/`` in the repository root.

    The root partition of the build machine is 99% full, so every generated
    artefact lives here (and ``bench-work/`` is gitignored).
    """
    env = os.environ.get("INFINITY_BENCH_WORK")
    return Path(env) if env else repo_root() / "bench-work"


def target_dirs() -> List[Path]:
    """Candidate cargo target directories, most likely first."""
    root = repo_root()
    out: List[Path] = []
    for cand in (
        os.environ.get("CARGO_TARGET_DIR"),
        "/tmp/infinity-target",  # set in .cargo/config.toml on this machine
        str(root / "target"),
    ):
        if not cand:
            continue
        p = Path(cand)
        if p not in out:
            out.append(p)
    return out


def infinity_candidates(binary: str = "infinity") -> List[Path]:
    """Every place the ``infinity`` binary might be, in preference order."""
    out: List[Path] = []
    env_bin = os.environ.get("INFINITY_BIN")
    if env_bin:
        out.append(Path(env_bin))
    for td in target_dirs():
        out.append(td / "release" / binary)
    for td in target_dirs():
        out.append(td / "debug" / binary)
    return out


def resolve_binary(
    explicit: Optional[str] = None, binary: str = "infinity", required: bool = True
) -> Optional[Path]:
    """Resolve the CLI binary.

    Order: ``--infinity`` argument, ``$INFINITY_BIN``, ``$CARGO_TARGET_DIR/release``,
    ``/tmp/infinity-target/release`` (this machine's target dir), ``target/release``,
    ``target/debug``, then ``binary`` on ``PATH``.
    """
    tried: List[str] = []
    if explicit:
        p = Path(explicit).expanduser()
        tried.append(str(p))
        if p.is_file() and os.access(p, os.X_OK):
            return p
        if required:
            raise PrereqError(
                f"{binary} binary not found or not executable at --infinity value {p}\n"
                f"       build it with: cargo build --release   (target dir: {target_dirs()[0]})"
            )
        return None
    for cand in infinity_candidates(binary):
        tried.append(str(cand))
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand
    on_path = shutil.which(binary)
    tried.append(f"{binary} on PATH")
    if on_path:
        return Path(on_path)
    if required:
        raise PrereqError(
            f"no {binary} binary found. Tried:\n         "
            + "\n         ".join(tried)
            + f"\n       build it with: cargo build --release\n"
            f"       or point at an existing one with --infinity <path>"
        )
    return None


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def human_bytes(n: Optional[float]) -> str:
    if n is None:
        return "n/a"
    n = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024 or unit == "TiB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TiB"


def human_seconds(s: Optional[float]) -> str:
    if s is None:
        return "n/a"
    s = float(s)
    if s < 1:
        return f"{s * 1000:.0f} ms"
    if s < 60:
        return f"{s:.2f} s"
    m, sec = divmod(s, 60)
    if m < 60:
        return f"{int(m)}m{sec:04.1f}s"
    h, m = divmod(int(m), 60)
    return f"{h}h{m:02d}m{sec:02.0f}s"


def fmt(value: Any, spec: str = ".3g") -> str:
    """Format a number, or ``n/a`` when it is unknown (never ``0``)."""
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int,)):
        return f"{value:d}"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "n/a"
        return format(value, spec)
    return str(value)


def ratio(off: Optional[float], base: Optional[float]) -> Optional[float]:
    """``off / base``, or ``None`` when either side is unknown or zero.

    A ratio is the input to the verdict rule of PROJECT.md section 11, so an
    unknown side must stay unknown rather than become 1.0.
    """
    if off is None or base is None:
        return None
    try:
        base_f = float(base)
        off_f = float(off)
    except (TypeError, ValueError):
        return None
    if base_f <= 0:
        return None
    return off_f / base_f


def pct_change(off: Optional[float], base: Optional[float]) -> Optional[float]:
    r = ratio(off, base)
    return None if r is None else (r - 1.0) * 100.0


def median(values: Sequence[float]) -> Optional[float]:
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return statistics.median(vals) if vals else None


def mean(values: Sequence[float]) -> Optional[float]:
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return statistics.fmean(vals) if vals else None


def quantile(values: Sequence[float], q: float) -> Optional[float]:
    vals = sorted(v for v in values if v is not None and not math.isnan(v))
    if not vals:
        return None
    if len(vals) == 1:
        return vals[0]
    pos = (len(vals) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[int(pos)]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def pearson(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    pairs = [
        (float(a), float(b))
        for a, b in zip(xs, ys)
        if a is not None and b is not None and not math.isnan(a) and not math.isnan(b)
    ]
    if len(pairs) < 2:
        return None
    xa = [p[0] for p in pairs]
    ya = [p[1] for p in pairs]
    mx, my = statistics.fmean(xa), statistics.fmean(ya)
    num = sum((x - mx) * (y - my) for x, y in pairs)
    dx = math.sqrt(sum((x - mx) ** 2 for x in xa))
    dy = math.sqrt(sum((y - my) ** 2 for y in ya))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def spearman(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    def ranks(vals: Sequence[float]) -> List[float]:
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        out = [0.0] * len(vals)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for t in range(i, j + 1):
                out[order[t]] = avg
            i = j + 1
        return out

    pairs = [(a, b) for a, b in zip(xs, ys) if a is not None and b is not None]
    if len(pairs) < 2:
        return None
    return pearson(ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs]))


# ---------------------------------------------------------------------------
# JSON and table IO
# ---------------------------------------------------------------------------


def load_json(path: Path) -> Any:
    if not path.exists():
        raise PrereqError(f"missing required JSON file: {path}")
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except json.JSONDecodeError as exc:
        raise PrereqError(f"{path} is not valid JSON ({exc})") from exc
    except OSError as exc:
        raise PrereqError(f"cannot read {path}: {exc}") from exc


def load_json_optional(path: Path) -> Tuple[Any, Optional[str]]:
    """``load_json`` that never raises: returns ``(obj_or_None, problem)``."""
    if not path.exists():
        return None, f"not found: {path}"
    try:
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh), None
    except (json.JSONDecodeError, OSError) as exc:
        return None, f"unreadable {path}: {exc}"


def dump_json(path: Path, obj: Any) -> Path:
    """Write JSON atomically so a crash cannot leave a half-written result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, sort_keys=False, default=str)
        fh.write("\n")
    os.replace(tmp, path)
    return path


def dig(obj: Any, dotted: str, default: Any = None) -> Any:
    """``dig(run_json, "timings.first_tree_ms")`` without ever raising."""
    cur = obj
    for part in dotted.split("."):
        if isinstance(cur, dict):
            if part not in cur:
                return default
            cur = cur[part]
        elif isinstance(cur, list):
            try:
                cur = cur[int(part)]
            except (ValueError, IndexError):
                return default
        else:
            return default
    return default if cur is None else cur


@dataclass
class Table:
    """A TSV/CSV file read defensively."""

    path: Path
    header: List[str]
    rows: List[Dict[str, str]]
    problems: List[str] = field(default_factory=list)
    delimiter: str = "\t"

    def __len__(self) -> int:
        return len(self.rows)

    def columns(self) -> List[str]:
        return list(self.header)

    def find_column(self, *aliases: str) -> Optional[str]:
        """First header entry whose normalised name matches one of ``aliases``."""
        norm = {normalize_key(h): h for h in self.header}
        for alias in aliases:
            hit = norm.get(normalize_key(alias))
            if hit is not None:
                return hit
        return None

    def require_column(self, *aliases: str) -> str:
        hit = self.find_column(*aliases)
        if hit is None:
            raise PrereqError(
                f"{self.path}: no column matching {list(aliases)}; "
                f"the header is {self.header}"
            )
        return hit

    def values(self, column: str) -> List[str]:
        return [r.get(column, "") for r in self.rows]


def normalize_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def looks_numeric(text: str) -> bool:
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


def read_table(
    path: Path,
    *,
    has_header: Optional[bool] = None,
    delimiter: Optional[str] = None,
    max_rows: Optional[int] = None,
) -> Table:
    """Read a TSV (or CSV) file without assuming the header is complete.

    ``has_header=None`` auto-detects: a first line whose fields are not all
    numeric is taken as a header. Missing trailing fields become ``""`` rather
    than an exception, and the fact is recorded in ``problems``.
    """
    if not path.exists():
        raise PrereqError(f"missing table file: {path}")
    delim = delimiter
    problems: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        first = fh.readline()
        if not first.strip():
            return Table(path, [], [], ["file is empty"], delim or "\t")
        if delim is None:
            delim = "\t" if "\t" in first else ","
        fh.seek(0)
        raw = [line.rstrip("\r\n") for line in fh]
    raw = [line for line in raw if line.strip() != ""]
    if not raw:
        return Table(path, [], [], ["file has no data rows"], delim)
    first_fields = raw[0].split(delim)
    if has_header is None:
        has_header = not all(looks_numeric(f) for f in first_fields if f != "")
        if not first_fields or all(f == "" for f in first_fields):
            has_header = False
    if has_header:
        header = [h.strip() for h in first_fields]
        body = raw[1:]
    else:
        width = max(len(line.split(delim)) for line in raw)
        header = [f"col{i}" for i in range(width)]
        body = raw
    seen: Dict[str, int] = {}
    uniq_header: List[str] = []
    for h in header:
        name = h or f"col{len(uniq_header)}"
        if name in seen:
            seen[name] += 1
            name = f"{name}__{seen[name]}"
        else:
            seen[name] = 0
        uniq_header.append(name)
    width = len(uniq_header)
    rows: List[Dict[str, str]] = []
    short = 0
    long = 0
    first_bad_line = None
    for i, line in enumerate(body):
        if max_rows is not None and len(rows) >= max_rows:
            problems.append(f"stopped after {max_rows} rows (max_rows)")
            break
        fields = line.split(delim)
        if len(fields) < width:
            short += 1
            first_bad_line = first_bad_line or (i + 2)
            fields = fields + [""] * (width - len(fields))
        elif len(fields) > width:
            long += 1
            first_bad_line = first_bad_line or (i + 2)
            fields = fields[:width]  # extra columns are dropped, and reported
        rows.append(dict(zip(uniq_header, (f.strip() for f in fields))))
    if short or long:
        detail = []
        if short:
            detail.append(f"{short} row(s) had fewer than {width} fields (blank-filled)")
        if long:
            detail.append(f"{long} row(s) had more than {width} fields (extra dropped)")
        problems.append("; ".join(detail) + f", first at line {first_bad_line}")
    return Table(path, uniq_header, rows, problems, delim)


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(str(h) for h in header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")
    return path


def write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(",".join(str(h) for h in header) + "\n")
        for row in rows:
            fh.write(
                ",".join(
                    "" if v is None else (str(v) if _csv_safe(v) else json.dumps(str(v)))
                    for v in row
                )
                + "\n"
            )
    return path


def _csv_safe(value: Any) -> bool:
    s = str(value)
    return not any(c in s for c in (",", '"', "\n", "\r"))


def markdown_table(header: Sequence[str], rows: Iterable[Sequence[Any]]) -> str:
    """A GitHub-flavoured markdown table. ``None`` renders as ``n/a``."""
    body = [[("" if c is None else str(c)) for c in row] for row in rows]
    out = ["| " + " | ".join(header) + " |", "|" + "|".join(["---"] * len(header)) + "|"]
    for row in body:
        cells = list(row) + [""] * (len(header) - len(row))
        out.append("| " + " | ".join(c.replace("|", "\\|") for c in cells[: len(header)]) + " |")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# run.json
# ---------------------------------------------------------------------------

TIMING_KEYS = (
    "scan_ms",
    "cluster_ms",
    "backbone_ms",
    "refine_ms",
    "output_ms",
    "total_ms",
    "first_tree_ms",
    "exact_tree_ms",
)

COUNTER_KEYS = (
    "blocks_looked_up",
    "blocks_matched",
    "pieces_compared",
    "cache_hits",
    "alignments_run",
    "bases_aligned",
    "kmers_scanned",
    "genomes_scanned",
    "end_extensions",
    "stray_points_dropped",
    "breakpoints",
    "blocks_total",
)


@dataclass
class RunJson:
    """The contents of a ``run.json``, read defensively.

    Every accessor returns ``None`` for something the running build does not
    emit yet, and says so in :attr:`problems`. The report prints those problems
    rather than a zero.
    """

    path: Path
    raw: Dict[str, Any]
    problems: List[str] = field(default_factory=list)

    # -- scalars ---------------------------------------------------------
    @property
    def version(self) -> Optional[str]:
        return self._get("version")

    @property
    def cancelled(self) -> Optional[bool]:
        return self._get("cancelled")

    @property
    def genomes(self) -> Optional[int]:
        return self._int("genomes")

    @property
    def clusters(self) -> Optional[int]:
        return self._int("clusters")

    @property
    def total_length(self) -> Optional[int]:
        return self._int("total_length")

    @property
    def focus_counts(self) -> List[Optional[int]]:
        raw = self.raw.get("focus_counts")
        if not isinstance(raw, (list, tuple)):
            self._note("focus_counts missing or not a list")
            return [None, None, None]
        out: List[Optional[int]] = []
        for i in range(3):
            try:
                out.append(int(raw[i]))
            except (IndexError, TypeError, ValueError):
                out.append(None)
        return out

    @property
    def peak_memory_bytes(self) -> Optional[int]:
        value = self.raw.get("peak_memory_bytes", None)
        if value is None:
            self._note("peak_memory_bytes is null (platform did not report it)")
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            self._note(f"peak_memory_bytes unreadable: {value!r}")
            return None

    # -- timings ---------------------------------------------------------
    def timing_ms(self, key: str) -> Optional[int]:
        if key not in TIMING_KEYS:
            raise ValueError(f"unknown timing key {key!r}")
        timings = self.raw.get("timings")
        if not isinstance(timings, dict):
            self._note("timings object missing")
            return None
        if key not in timings:
            self._note(f"timings.{key} missing")
            return None
        value = timings[key]
        if value is None:
            self._note(f"timings.{key} is null (stage never completed)")
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            self._note(f"timings.{key} unreadable: {value!r}")
            return None

    def timing_s(self, key: str) -> Optional[float]:
        ms = self.timing_ms(key)
        return None if ms is None else ms / 1000.0

    @property
    def first_tree_s(self) -> Optional[float]:
        return self.timing_s("first_tree_ms")

    @property
    def exact_tree_s(self) -> Optional[float]:
        return self.timing_s("exact_tree_ms")

    @property
    def total_s(self) -> Optional[float]:
        return self.timing_s("total_ms")

    # -- counters --------------------------------------------------------
    def counter(self, key: str) -> Any:
        counters = self.raw.get("counters")
        if not isinstance(counters, dict):
            self._note("counters object missing")
            return None
        if key not in counters:
            self._note(f"counters.{key} missing")
            return None
        return counters[key]

    def counter_int(self, key: str) -> Optional[int]:
        value = self.counter(key)
        if value is None:
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            self._note(f"counters.{key} unreadable: {value!r}")
            return None

    @property
    def blocks_matched(self) -> List[Optional[int]]:
        raw = self.counter("blocks_matched")
        if not isinstance(raw, (list, tuple)):
            return [None, None, None]
        out: List[Optional[int]] = []
        for i in range(3):
            try:
                out.append(int(raw[i]))
            except (IndexError, TypeError, ValueError):
                out.append(None)
        return out

    def block_identity(self) -> Optional[float]:
        """Fraction of block lookups that hit, or ``None`` when unknown."""
        looked = self.counter_int("blocks_looked_up")
        matched = [m for m in self.blocks_matched if m is not None]
        if not looked or not matched:
            return None
        return sum(matched) / looked

    def param(self, key: str) -> Any:
        params = self.raw.get("params")
        return params.get(key) if isinstance(params, dict) else None

    # -- helpers ---------------------------------------------------------
    def _get(self, key: str) -> Any:
        if key not in self.raw:
            self._note(f"{key} missing from run.json")
            return None
        return self.raw[key]

    def _int(self, key: str) -> Optional[int]:
        value = self._get(key)
        if value is None:
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            self._note(f"{key} unreadable: {value!r}")
            return None

    def _note(self, message: str) -> None:
        if message not in self.problems:
            self.problems.append(message)

    def unique_problems(self) -> List[str]:
        out: List[str] = []
        for p in self.problems:
            if p not in out:
                out.append(p)
        return out


def load_run_json(path: Path, *, required: bool = True) -> Optional[RunJson]:
    """Read ``run.json`` from a run output folder.

    With ``required=False`` a missing or corrupt file yields ``None`` plus a
    problem string on the returned object's caller side; with ``required=True``
    it raises :class:`PrereqError`.
    """
    if isinstance(path, RunJson):
        return path
    p = Path(path)
    if p.is_dir():
        p = p / "run.json"
    if not p.exists():
        if required:
            raise PrereqError(
                f"no run.json at {p}\n"
                "       this output folder is incomplete: re-run the configuration "
                "(run_matrix.py --force) or point at a finished run"
            )
        return None
    obj, problem = load_json_optional(p)
    if problem is not None:
        if required:
            raise PrereqError(problem)
        return None
    if not isinstance(obj, dict):
        if required:
            raise PrereqError(f"{p}: expected a JSON object, got {type(obj).__name__}")
        return None
    return RunJson(path=p, raw=obj)


# ---------------------------------------------------------------------------
# events.jsonl
# ---------------------------------------------------------------------------

#: Event names as serde writes them (``rename_all = "snake_case"``).
EVENT_NAMES = (
    "genome_scanned",
    "clustered",
    "backbone_tree",
    "cluster_focus",
    "map_update",
    "counters",
    "progress",
    "warning",
    "finished",
)


@dataclass
class EventLog:
    """A parsed ``events.jsonl``."""

    path: Path
    events: List[Dict[str, Any]]
    problems: List[str] = field(default_factory=list)
    lines_read: int = 0
    lines_bad: int = 0

    def of_type(self, *names: str) -> List[Dict[str, Any]]:
        wanted = set(names)
        return [e for e in self.events if event_name(e) in wanted]

    def first(self, name: str) -> Optional[Dict[str, Any]]:
        for e in self.events:
            if event_name(e) == name:
                return e
        return None

    def last(self, name: str) -> Optional[Dict[str, Any]]:
        found = None
        for e in self.events:
            if event_name(e) == name:
                found = e
        return found


def event_name(envelope: Dict[str, Any]) -> Optional[str]:
    """The event name of one ``events.jsonl`` record, tolerating layout drift."""
    ev = envelope.get("event")
    if isinstance(ev, dict):
        name = ev.get("event")
        if isinstance(name, str):
            return name
        # Some builds may inline the payload; fall back to a single known key.
        for key in EVENT_NAMES:
            if key in ev:
                return key
    if isinstance(ev, str):
        return ev
    for key in ("type", "kind", "name"):
        if isinstance(envelope.get(key), str):
            return envelope[key]
    return None


def event_data(envelope: Dict[str, Any]) -> Dict[str, Any]:
    """The payload of one record (``event.data``), or the record itself."""
    ev = envelope.get("event")
    if isinstance(ev, dict):
        data = ev.get("data")
        if isinstance(data, dict):
            return data
        return {k: v for k, v in ev.items() if k not in ("event", "data")}
    return {k: v for k, v in envelope.items() if k not in ("seq", "t_ms", "wall_ms", "event")}


def event_t_ms(envelope: Dict[str, Any]) -> Optional[float]:
    for key in ("t_ms", "elapsed_ms", "t"):
        value = envelope.get(key)
        if value is None and isinstance(envelope.get("event"), dict):
            value = envelope["event"].get(key)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None
    return None


def read_event_log(path: Path, *, max_events: Optional[int] = None) -> EventLog:
    """Parse ``events.jsonl`` line by line, keeping every bad line as a problem."""
    p = Path(path)
    if p.is_dir():
        p = p / "events.jsonl"
    if not p.exists():
        raise PrereqError(
            f"no events.jsonl at {p}\n"
            "       the run did not produce an event stream; the accuracy-over-time "
            "figure needs one (PROJECT.md section 7)"
        )
    events: List[Dict[str, Any]] = []
    problems: List[str] = []
    bad = 0
    read = 0
    with open(p, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            read += 1
            if max_events is not None and len(events) >= max_events:
                problems.append(f"stopped after {max_events} events (max_events)")
                break
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                bad += 1
                if bad <= 3:
                    problems.append(f"line {lineno}: not JSON ({exc})")
                continue
            if not isinstance(obj, dict):
                bad += 1
                if bad <= 3:
                    problems.append(f"line {lineno}: JSON is not an object")
                continue
            if event_name(obj) is None:
                bad += 1
                if bad <= 3:
                    problems.append(
                        f"line {lineno}: no recognisable event name; keys were {sorted(obj)}"
                    )
                continue
            obj.setdefault("_line", lineno)
            events.append(obj)
    if bad > 3:
        problems.append(f"{bad} lines in total were unparsable and were skipped")
    if not events:
        problems.append("no usable events found")
    return EventLog(path=p, events=events, problems=problems, lines_read=read, lines_bad=bad)


# ---------------------------------------------------------------------------
# Newick parsing, branch collapsing, Robinson-Foulds
# ---------------------------------------------------------------------------


class NewickNode:
    """A node of a parsed Newick tree."""

    __slots__ = ("name", "length", "children", "parent")

    def __init__(
        self,
        name: str = "",
        length: Optional[float] = None,
        parent: Optional["NewickNode"] = None,
    ) -> None:
        self.name = name
        self.length = length
        self.children: List["NewickNode"] = []
        self.parent = parent

    def add(self, child: "NewickNode") -> "NewickNode":
        child.parent = self
        self.children.append(child)
        return child

    def is_leaf(self) -> bool:
        return not self.children

    def copy(self) -> "NewickNode":
        node = NewickNode(self.name, self.length)
        for child in self.children:
            node.add(child.copy())
        return node

    def iter_nodes(self) -> Iterator["NewickNode"]:
        yield self
        for child in self.children:
            yield from child.iter_nodes()

    def leaves(self) -> List["NewickNode"]:
        return [n for n in self.iter_nodes() if n.is_leaf()]

    def leaf_names(self) -> List[str]:
        return [n.name for n in self.leaves() if n.name]

    def total_length(self) -> float:
        return sum(
            n.length for n in self.iter_nodes() if n is not self and n.length is not None
        )

    def find(self, name: str) -> Optional["NewickNode"]:
        for node in self.iter_nodes():
            if node.name == name:
                return node
        return None

    def to_newick(self) -> str:
        if self.is_leaf():
            body = _quote_name(self.name)
        else:
            body = "(" + ",".join(c.to_newick() for c in self.children) + ")"
            if self.name:
                body += _quote_name(self.name)
        if self.length is not None:
            body += f":{self.length!r}"
        return body


def _quote_name(name: str) -> str:
    if not name:
        return ""
    if any(c in name for c in "()[],:;'\" "):
        return "'" + name.replace("'", "''") + "'"
    return name


_NEWICK_COMMENT = re.compile(r"\[[^\]]*\]")


def parse_newick(text: str) -> NewickNode:
    """Parse one Newick string.

    Tolerates internal node labels (``)100:0.01``), Nexus-style comments
    (``[&...]``), quoted names, a trailing ``;`` and surrounding whitespace.
    Raises :class:`PrereqError` with the offending position on a syntax error,
    because a malformed tree is a prerequisite problem, not a data problem.
    """
    if text is None:
        raise PrereqError("cannot parse an empty Newick string")
    cleaned = _NEWICK_COMMENT.sub("", str(text)).strip()
    if not cleaned:
        raise PrereqError("cannot parse an empty Newick string")
    # Take the first tree if the string holds several.
    depth = 0
    end = len(cleaned)
    for i, ch in enumerate(cleaned):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == ";" and depth <= 0:
            end = i
            break
    cleaned = cleaned[:end].strip().rstrip(";").strip()
    if not cleaned:
        raise PrereqError(f"Newick string has no tree: {text[:80]!r}")
    pos = 0
    length = len(cleaned)

    def parse_node() -> NewickNode:
        nonlocal pos
        node = NewickNode()
        if pos < length and cleaned[pos] == "(":
            pos += 1
            node.add(parse_node())
            while pos < length and cleaned[pos] == ",":
                pos += 1
                node.add(parse_node())
            if pos >= length or cleaned[pos] != ")":
                raise PrereqError(
                    f"Newick parse error at offset {pos}: expected ')' in {cleaned[:120]!r}"
                )
            pos += 1
        # label and/or length
        start = pos
        while pos < length and cleaned[pos] not in "(),:;":
            if cleaned[pos] == "'":
                pos += 1
                while pos < length and cleaned[pos] != "'":
                    pos += 1
            pos += 1
        label = cleaned[start:pos].strip()
        if label.startswith("'") and label.endswith("'") and len(label) >= 2:
            label = label[1:-1].replace("''", "'")
        if pos < length and cleaned[pos] == ":":
            pos += 1
            nstart = pos
            while pos < length and cleaned[pos] not in "(),:;":
                pos += 1
            token = cleaned[nstart:pos].strip()
            try:
                node.length = float(token)
            except ValueError:
                node.length = None
        node.name = label
        return node

    root = parse_node()
    if root.name and root.is_leaf() is False:
        pass  # internal label, keep as-is
    return root


def write_newick(node: NewickNode, semicolon: bool = True) -> str:
    return node.to_newick() + (";" if semicolon else "")


def read_newick_file(path: Path) -> NewickNode:
    p = Path(path)
    if p.is_dir():
        p = p / "tree.nwk"
    if not p.exists():
        raise PrereqError(f"no tree file at {p}")
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise PrereqError(f"cannot read tree {p}: {exc}") from exc
    if not text.strip():
        raise PrereqError(f"tree file {p} is empty")
    return parse_newick(text)


def contract_short_branches(node: NewickNode, threshold: float) -> NewickNode:
    """Return a copy of ``node`` with every branch shorter than ``threshold``
    contracted (its node suppressed, children re-attached to the parent).

    PROJECT.md section 11 asks for the Robinson-Foulds distance "after
    collapsing branches shorter than half a SNP": such a branch carries no
    evidence, so the split it implies must not count.
    """
    if threshold is None or threshold <= 0:
        return node.copy()

    def build(src: NewickNode, is_root: bool) -> NewickNode:
        out = NewickNode("" if is_root else src.name, src.length)
        kept: List[NewickNode] = []
        for child in src.children:
            if not child.is_leaf() and child.length is not None and child.length < threshold:
                # Contract: the child's subtree is spliced into this node and
                # the split the child defined disappears.
                spliced = build(child, is_root=True)
                for gc in spliced.children:
                    gc.parent = out
                    kept.append(gc)
            else:
                built = build(child, is_root=False)
                built.parent = out
                kept.append(built)
        out.children = kept
        return out

    contracted = build(node, is_root=True)
    # A root left with exactly two children is an unrooted artefact; keep it,
    # because splits() canonicalises complements anyway.
    return contracted


def _leaf_set(node: NewickNode) -> frozenset:
    return frozenset(n.name for n in node.leaves() if n.name)


def splits(node: NewickNode, leaves: Optional[frozenset] = None) -> set:
    """Non-trivial bipartitions of ``node`` (unrooted, canonicalised).

    Each bipartition is keyed by the side that does *not* hold the
    alphabetically first leaf. Keying by "the smaller side" instead would give
    a bipartition with two equal halves two different keys (one per side), so
    the two children of a root splitting the leaves 50/50 would be counted as
    two splits.
    """
    all_leaves = leaves if leaves is not None else _leaf_set(node)
    if len(all_leaves) < 4:
        return set()
    ref = min(all_leaves)
    out = set()
    for child_node in node.iter_nodes():
        if child_node is node:
            continue
        side = _leaf_set(child_node) & all_leaves
        if not side or side >= all_leaves:
            continue
        complement = all_leaves - side
        if len(side) < 2 or len(complement) < 2:
            continue  # trivial split (one leaf on one side)
        out.add(frozenset(complement if ref in side else side))
    return out


# -- bitmask splits ----------------------------------------------------------
#
# Splits as Python ints (one bit per leaf) instead of frozensets of names. The
# accuracy-over-time curve scores hundreds of trees of up to thousands of
# leaves per run, and the set-of-names representation is quadratic in the
# depth of the tree; with masks a whole tree is one post-order pass.

_INT_BIT_COUNT = getattr(int, "bit_count", None)


def popcount(value: int) -> int:
    """Number of set bits (``int.bit_count`` where Python has it)."""
    if _INT_BIT_COUNT is not None:
        return _INT_BIT_COUNT(value)
    return bin(value).count("1")


def leaf_index(names: Iterable[str]) -> Dict[str, int]:
    """``name -> bit`` for a leaf universe, in sorted name order."""
    return {name: 1 << i for i, name in enumerate(sorted(set(names)))}


def clade_masks(
    root: NewickNode, index: Dict[str, int], threshold: float = 0.0
) -> Tuple[int, List[int]]:
    """``(leaf mask of the whole tree, [clade mask of every kept edge])``.

    One entry per non-root internal node whose branch is *not* contracted by
    ``threshold`` (the same rule as :func:`contract_short_branches`: a branch
    shorter than the threshold is contracted, a missing length is kept, a leaf
    is never contracted). Leaves missing from ``index`` contribute no bit, which
    is how a tree is restricted to a leaf universe. Iterative, so a
    caterpillar tree with thousands of leaves cannot hit the recursion limit.
    """
    masks: List[int] = []
    value: Dict[int, int] = {}
    stack: List[Tuple[NewickNode, bool]] = [(root, False)]
    while stack:
        node, expanded = stack.pop()
        if not node.children:
            value[id(node)] = index.get(node.name, 0) if node.name else 0
            continue
        if not expanded:
            stack.append((node, True))
            for child in node.children:
                stack.append((child, False))
            continue
        mask = 0
        for child in node.children:
            mask |= value.pop(id(child), 0)
        value[id(node)] = mask
        if node is root:
            continue
        if threshold and node.length is not None and node.length < threshold:
            continue
        masks.append(mask)
    return value.get(id(root), 0), masks


def canonical_splits(masks: Iterable[int], universe: int) -> set:
    """Non-trivial bipartitions of ``universe`` among ``masks``, canonicalised.

    Each mask is restricted to ``universe``; a bipartition is kept when both
    sides hold at least two leaves, and keyed by the side without the lowest
    set bit of ``universe`` - so both sides of one bipartition map to one key.
    """
    n = popcount(universe)
    if n < 4:
        return set()
    ref = universe & -universe
    out = set()
    for mask in masks:
        mask &= universe
        k = popcount(mask)
        if k < 2 or n - k < 2:
            continue
        out.add(universe ^ mask if mask & ref else mask)
    return out


def prune_to(node: NewickNode, keep: set) -> NewickNode:
    """A copy of ``node`` restricted to the leaves in ``keep``."""
    def build(src: NewickNode) -> Optional[NewickNode]:
        if src.is_leaf():
            if src.name in keep:
                return NewickNode(src.name, src.length)
            return None
        children = [c for c in (build(ch) for ch in src.children) if c is not None]
        if not children:
            return None
        if len(children) == 1:
            only = children[0]
            if src is node:
                return only
            merged = NewickNode(only.name or src.name, None)
            for gc in only.children:
                merged.add(gc)
            return merged if merged.children else only
        out = NewickNode("" if src is node else src.name, None if src is node else src.length)
        for c in children:
            out.add(c)
        return out

    result = build(node)
    return result if result is not None else NewickNode()


@dataclass
class RfResult:
    """One Robinson-Foulds comparison."""

    rf: Optional[int] = None
    rf_normalized: Optional[float] = None
    splits_a: Optional[int] = None
    splits_b: Optional[int] = None
    shared_splits: Optional[int] = None
    leaves_a: Optional[int] = None
    leaves_b: Optional[int] = None
    leaves_shared: Optional[int] = None
    collapse_threshold: Optional[float] = None
    problems: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.rf is not None

    def as_dict(self) -> Dict[str, Any]:
        return {
            "rf": self.rf,
            "rf_normalized": self.rf_normalized,
            "splits_run": self.splits_a,
            "splits_truth": self.splits_b,
            "splits_shared": self.shared_splits,
            "leaves_run": self.leaves_a,
            "leaves_truth": self.leaves_b,
            "leaves_shared": self.leaves_shared,
            "collapse_threshold": self.collapse_threshold,
            "problems": list(self.problems),
        }


def rf_distance(
    run_tree: NewickNode,
    truth_tree: NewickNode,
    *,
    collapse: float = 0.0,
) -> RfResult:
    """Robinson-Foulds distance between two trees.

    Both trees are restricted to their shared leaf set first: a run that has
    only emitted part of the tree yet is compared on what it has, and
    ``leaves_shared`` says how much that was. ``rf_normalized`` is
    ``rf / (splits_run + splits_truth)`` in ``[0, 1]`` where 0 means identical.
    """
    result = RfResult(collapse_threshold=collapse or None)
    if run_tree is None or truth_tree is None:
        result.problems.append("missing tree")
        return result
    run_leaves = set(n.name for n in run_tree.leaves() if n.name)
    truth_leaves = set(n.name for n in truth_tree.leaves() if n.name)
    result.leaves_a = len(run_leaves)
    result.leaves_b = len(truth_leaves)
    shared = run_leaves & truth_leaves
    result.leaves_shared = len(shared)
    if not shared:
        result.problems.append(
            f"no shared leaf names (run has {len(run_leaves)}, truth has {len(truth_leaves)}); "
            "genome names must match between the run and the truth tree"
        )
        return result
    if len(shared) < len(truth_leaves):
        result.problems.append(
            f"only {len(shared)} of {len(truth_leaves)} truth leaves are in this tree"
        )
    if len(shared) < 4:
        result.problems.append(
            f"only {len(shared)} shared leaves: RF is undefined for fewer than 4 taxa"
        )
        return result
    # Contract short branches, restrict to the shared leaves, take the splits:
    # restricting every clade to the shared leaves gives exactly the splits of
    # the pruned tree, so no pruned copy is built.
    index = leaf_index(shared)
    universe = (1 << len(index)) - 1
    _, run_masks = clade_masks(run_tree, index, collapse or 0.0)
    _, truth_masks = clade_masks(truth_tree, index, collapse or 0.0)
    a = canonical_splits(run_masks, universe)
    b = canonical_splits(truth_masks, universe)
    return rf_from_split_sets(a, b, result)


def rf_from_split_sets(a: set, b: set, result: Optional[RfResult] = None) -> RfResult:
    """Fill an :class:`RfResult` from two canonical split sets.

    ``rf_normalized = rf / (splits_a + splits_b)``: 0 is identical, 1 means no
    split is shared (a star tree against any resolved truth scores 1).
    """
    result = result or RfResult()
    result.splits_a = len(a)
    result.splits_b = len(b)
    result.shared_splits = len(a & b)
    result.rf = len(a ^ b)
    denom = len(a) + len(b)
    result.rf_normalized = (result.rf / denom) if denom else 0.0
    return result


def guess_collapse_threshold(
    truth_tree: Optional[NewickNode],
    *,
    collapse_snps: Optional[float] = 0.5,
    comparable_bp: Optional[int] = None,
    mode: str = "auto",
    explicit_length: Optional[float] = None,
) -> Tuple[float, str]:
    """Work out the branch-length threshold for "collapse < half a SNP".

    Infinity-Win writes substitutions per site inside clusters (section 6.9) while the
    simulator's truth tree is in SNPs (section 10.1), so one number cannot serve
    both. Returns ``(threshold, explanation)``; ``threshold == 0`` means "do not
    collapse", which is recorded rather than silently applied.
    """
    if mode == "none":
        return 0.0, "collapsing disabled (--collapse-mode none)"
    if mode == "length":
        if explicit_length is None:
            raise PrereqError("--collapse-mode length needs --collapse-length <float>")
        return float(explicit_length), f"collapsing branches shorter than {explicit_length}"
    if collapse_snps is None:
        return 0.0, "no collapse: --collapse-snps not given"
    if comparable_bp and comparable_bp > 0:
        threshold = float(collapse_snps) / float(comparable_bp)
        return (
            threshold,
            f"collapsing branches shorter than {collapse_snps} SNP "
            f"= {threshold:.3g} substitutions/site over {comparable_bp} comparable bp",
        )
    if mode == "snps":
        raise PrereqError(
            "--collapse-mode snps needs a comparable site count "
            "(--comparable-bp N, or a run.json total_length, or cluster_*.dist.tsv)"
        )
    # auto: infer the unit from the truth tree's own branch lengths.
    if truth_tree is None:
        return 0.0, "no collapse: no truth tree to infer branch-length units from"
    lengths = [
        n.length
        for n in truth_tree.iter_nodes()
        if n.length is not None and n.length > 0
    ]
    if not lengths:
        return 0.0, "no collapse: truth tree has no branch lengths"
    total = sum(lengths)
    looks_like_snps = total >= 1.0 and all(abs(l - round(l)) < 0.02 for l in lengths[:200])
    if looks_like_snps:
        return (
            float(collapse_snps),
            f"truth branch lengths look like SNP counts (total {total:g}); "
            f"collapsing branches shorter than {collapse_snps} SNP",
        )
    return (
        0.0,
        "no collapse: truth branch lengths are not SNP counts and no comparable "
        "site count was supplied, so 'half a SNP' cannot be converted "
        "(pass --comparable-bp N or --collapse-mode length X)",
    )


# ---------------------------------------------------------------------------
# Neighbour-joining (used by the Mash baseline and by nothing else)
# ---------------------------------------------------------------------------


def neighbor_joining(labels: Sequence[str], distance: Callable[[str, str], float]) -> str:
    """Classical neighbour-joining, returning a Newick string.

    ``distance(a, b)`` must be symmetric. Negative branch lengths are clamped to
    zero, as PROJECT.md section 6.9 requires. Implemented here because the
    outside baselines must not depend on an external phylogenetics package.
    """
    names = list(labels)
    if not names:
        return "();"
    if len(names) == 1:
        return f"{_quote_name(names[0])};"
    d: Dict[Tuple[str, str], float] = {}

    def key(a: str, b: str) -> Tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            d[key(a, b)] = float(distance(a, b))
    active = list(names)
    newick: Dict[str, str] = {n: _quote_name(n) for n in names}
    next_id = 0
    while len(active) > 2:
        n = len(active)
        row_sum = {a: sum(d.get(key(a, b), 0.0) for b in active if b != a) for a in active}
        best_pair = None
        best_q = None
        for i, a in enumerate(active):
            for b in active[i + 1 :]:
                q = (n - 2) * d.get(key(a, b), 0.0) - row_sum[a] - row_sum[b]
                if best_q is None or q < best_q:
                    best_q = q
                    best_pair = (a, b)
        assert best_pair is not None
        a, b = best_pair
        dab = d.get(key(a, b), 0.0)
        la = 0.5 * dab + (row_sum[a] - row_sum[b]) / (2.0 * (n - 2))
        lb = dab - la
        node_name = f"__nj{next_id}"
        next_id += 1
        newick[node_name] = (
            f"({newick[a]}:{max(la, 0.0):.6g},{newick[b]}:{max(lb, 0.0):.6g})"
        )
        for other in active:
            if other in (a, b):
                continue
            d[key(node_name, other)] = 0.5 * (
                d.get(key(a, other), 0.0) + d.get(key(b, other), 0.0) - dab
            )
        for k in [k for k in d if a in k or b in k]:
            del d[k]
        active = [x for x in active if x not in (a, b)] + [node_name]
    a, b = active
    return f"({newick[a]}:{max(d.get(key(a, b), 0.0) / 2.0, 0.0):.6g}," \
           f"{newick[b]}:{max(d.get(key(a, b), 0.0) / 2.0, 0.0):.6g});"


# ---------------------------------------------------------------------------
# SKA2 distance files (the exact alignment-free outside baseline)
# ---------------------------------------------------------------------------


@dataclass
class SkaPair:
    """One row of a SKA2 ``*_distance.tsv`` file."""

    sample1: str
    sample2: str
    #: SKA2's own "Distance" column: SNP-like distance, can be fractional and
    #: can be ``NA`` for pairs with too few shared split k-mers.
    distance: Optional[float]
    mismatch_proportion: Optional[float]
    match_count: Optional[int]
    mismatch_count: Optional[int]
    raw: Dict[str, str] = field(default_factory=dict)

    def key(self) -> Tuple[str, str]:
        return pair_key(self.sample1, self.sample2)


SKA2_HEADER_ALIASES = {
    "sample1": ("sample1", "sample_1", "samplea", "query", "genome1", "a"),
    "sample2": ("sample2", "sample_2", "sampleb", "subject", "genome2", "b"),
    "distance": ("distance", "dist", "snp_distance", "ska_distance"),
    "mismatch_proportion": (
        "mismatches_proportion",
        "mismatch_proportion",
        "mismatches (proportion)",
        "mismatches",
        "proportion",
    ),
    "match_count": ("match_count", "matches", "match", "matched"),
    "mismatch_count": ("mismatch_count", "mismatches_count", "mismatch", "mismatched"),
}


def pair_key(a: str, b: str) -> Tuple[str, str]:
    """Canonical, order-independent key for a pair of genome names."""
    return (a, b) if str(a) <= str(b) else (b, a)


def _to_float(text: Optional[str]) -> Optional[float]:
    if text is None:
        return None
    t = str(text).strip()
    if t == "" or t.upper() in {"NA", "NAN", "N/A", "NULL", "INF", "-INF"}:
        return None
    try:
        value = float(t)
    except ValueError:
        return None
    return None if math.isnan(value) else value


def _to_int(text: Optional[str]) -> Optional[int]:
    value = _to_float(text)
    return None if value is None else int(round(value))


@dataclass
class SkaDistances:
    """A parsed SKA2 distance table."""

    path: Path
    pairs: List[SkaPair]
    problems: List[str] = field(default_factory=list)
    columns: Dict[str, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.pairs)

    def as_dict(self) -> Dict[Tuple[str, str], SkaPair]:
        return {p.key(): p for p in self.pairs}

    def samples(self) -> List[str]:
        out: List[str] = []
        for p in self.pairs:
            for s in (p.sample1, p.sample2):
                if s not in out:
                    out.append(s)
        return sorted(out)

    def distances(self) -> List[float]:
        return [p.distance for p in self.pairs if p.distance is not None]

    def stats(self) -> Dict[str, Any]:
        ds = self.distances()
        samples = self.samples()
        return {
            "rows": len(self.pairs),
            "rows_with_distance": len(ds),
            "rows_missing_distance": len(self.pairs) - len(ds),
            "sample_count": len(samples),
            "samples": samples,
            "min": min(ds) if ds else None,
            "max": max(ds) if ds else None,
            "mean": mean(ds),
            "median": median(ds),
            "problems": list(self.problems),
        }


def parse_ska2_distance_tsv(path: Path) -> SkaDistances:
    """Parse a SKA2 ``ska distance`` TSV.

    Expected header (as produced by SKA2 1.x on the real dataset)::

        Sample1  Sample2  Distance  Mismatches (proportion)  Match count  Mismatch count

    Column matching is by normalised alias, order does not matter, and unknown
    extra columns are kept in :attr:`SkaPair.raw`. Rows whose distance is ``NA``
    are kept with ``distance=None`` and counted in ``problems``, never dropped
    silently.
    """
    p = Path(path)
    if not p.exists():
        raise PrereqError(
            f"SKA2 distance file not found: {p}\n"
            "       SKA2 is not installed on this machine, so a pre-computed file "
            "must be supplied with --ska-dist <path>"
        )
    table = read_table(p)
    problems: List[str] = list(table.problems)
    columns: Dict[str, str] = {}
    for field_name, aliases in SKA2_HEADER_ALIASES.items():
        hit = table.find_column(*aliases)
        if hit is not None:
            columns[field_name] = hit
    for required in ("sample1", "sample2"):
        if required not in columns:
            raise PrereqError(
                f"{p}: cannot find the {required} column; header was {table.header}. "
                "Expected a SKA2 'ska distance' TSV."
            )
    if "distance" not in columns:
        problems.append(
            f"no distance column in {p.name} (header {table.header}); "
            "distances will be None"
        )
    pairs: List[SkaPair] = []
    missing_distance = 0
    for i, row in enumerate(table.rows, start=2):
        s1 = row.get(columns["sample1"], "").strip()
        s2 = row.get(columns.get("sample2", ""), "").strip()
        if not s1 or not s2:
            problems.append(f"line {i}: empty sample name, row skipped")
            continue
        if s1 == s2:
            problems.append(f"line {i}: self-pair {s1} skipped")
            continue
        dist = _to_float(row.get(columns.get("distance", ""), None)) if "distance" in columns else None
        if dist is None:
            missing_distance += 1
        pairs.append(
            SkaPair(
                sample1=s1,
                sample2=s2,
                distance=dist,
                mismatch_proportion=_to_float(row.get(columns.get("mismatch_proportion", "")))
                if "mismatch_proportion" in columns
                else None,
                match_count=_to_int(row.get(columns.get("match_count", "")))
                if "match_count" in columns
                else None,
                mismatch_count=_to_int(row.get(columns.get("mismatch_count", "")))
                if "mismatch_count" in columns
                else None,
                raw=dict(row),
            )
        )
    if missing_distance:
        problems.append(
            f"{missing_distance} of {len(pairs)} rows had no usable Distance value "
            "(SKA2 writes NA when too few split k-mers are shared)"
        )
    if not pairs:
        problems.append("no usable rows")
    return SkaDistances(path=p, pairs=pairs, problems=problems, columns=columns)


# ---------------------------------------------------------------------------
# Infinity-Win output folder readers
# ---------------------------------------------------------------------------

GENOME_NAME_ALIASES = ("name", "genome", "genome_name", "sample", "id", "genome_id")

#: Column aliases for ``cluster_<id>.dist.tsv``.
#:
#: The real file writes **both** a numeric id and a name per side::
#:
#:     genome_a  name_a  genome_b  name_b  snps  indels  comparable_bp  subs_per_site
#:
#: Joins against SKA2, a truth tree and the published expected values are all by
#: genome *name*, so the name columns are matched first and the id columns are
#: kept separately. Reading ``genome_a`` as though it were a name would produce a
#: comparison that silently matches nothing.
DIST_COLUMN_ALIASES = {
    "a": (
        "name_a", "genome_name_a", "sample_a", "sample1", "a", "left", "query",
        "genome_a", "genome1", "genome", "member",
    ),
    "b": (
        "name_b", "genome_name_b", "sample_b", "sample2", "b", "right", "subject",
        "genome_b", "genome2", "other", "member_b",
    ),
    "a_id": ("genome_a", "genome1", "id_a", "genome_id_a"),
    "b_id": ("genome_b", "genome2", "id_b", "genome_id_b"),
    "snps": ("snps", "snp", "snp_count", "snps_ab", "d_snps", "distance_snps", "snp_distance"),
    "comparable_bp": (
        "comparable_bp",
        "comparable",
        "comparable_bases",
        "sites",
        "n_sites",
        "compared_bp",
    ),
    "distance": (
        "subs_per_site",
        "distance",
        "dist",
        "d",
        "substitutions_per_site",
        "per_site",
    ),
    "focus": ("focus", "focus_level", "level"),
    "indels": ("indels", "indel", "indel_count"),
}


def load_cluster_distances(outdir: Path) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Read every ``cluster_<id>.dist.tsv`` in an output folder.

    Returns ``(rows, problems)`` where each row is a dict with keys ``cluster``,
    ``a``, ``b``, ``snps``, ``comparable_bp``, ``distance``, ``focus``. Column
    names are matched by alias; when the two genome columns cannot be identified
    the file is reported and skipped instead of guessed at.
    """
    outdir = Path(outdir)
    problems: List[str] = []
    rows: List[Dict[str, Any]] = []
    if not outdir.is_dir():
        return rows, [f"output folder does not exist: {outdir}"]
    files = sorted(outdir.glob("cluster_*.dist.tsv"))
    if not files:
        problems.append(f"no cluster_*.dist.tsv in {outdir}")
        return rows, problems
    for path in files:
        cluster = path.name[len("cluster_") : -len(".dist.tsv")]
        try:
            table = read_table(path)
        except PrereqError as exc:
            problems.append(str(exc))
            continue
        problems.extend(f"{path.name}: {p}" for p in table.problems)
        cols: Dict[str, Optional[str]] = {
            field: table.find_column(*aliases) for field, aliases in DIST_COLUMN_ALIASES.items()
        }
        if cols["a"] is None or cols["b"] is None:
            # A long-format file may instead be a matrix; report clearly.
            problems.append(
                f"{path.name}: cannot identify the two genome columns "
                f"(header {table.header}); file skipped"
            )
            continue
        if cols["snps"] is None and cols["distance"] is None:
            problems.append(
                f"{path.name}: neither a SNP-count nor a distance column found "
                f"(header {table.header}); file skipped"
            )
            continue
        for row in table.rows:
            a = row.get(cols["a"], "").strip()
            b = row.get(cols["b"], "").strip()
            if not a or not b:
                continue
            a_id = row.get(cols["a_id"]).strip() if cols.get("a_id") else None
            b_id = row.get(cols["b_id"]).strip() if cols.get("b_id") else None
            rows.append(
                {
                    "cluster": cluster,
                    "a": a,
                    "b": b,
                    # The numeric ids, kept apart from the names: a join on the
                    # wrong one matches nothing and looks like a clean zero.
                    "a_id": a_id if a_id != a else None,
                    "b_id": b_id if b_id != b else None,
                    "snps": _to_int(row.get(cols["snps"])) if cols["snps"] else None,
                    "comparable_bp": _to_int(row.get(cols["comparable_bp"]))
                    if cols["comparable_bp"]
                    else None,
                    "distance": _to_float(row.get(cols["distance"])) if cols["distance"] else None,
                    "focus": row.get(cols["focus"]) if cols["focus"] else None,
                    "indels": _to_int(row.get(cols["indels"])) if cols["indels"] else None,
                    "_source": path.name,
                }
            )
    return rows, problems


def distance_lookup(rows: Iterable[Dict[str, Any]], value: str = "snps") -> Dict[Tuple[str, str], Any]:
    """Index ``load_cluster_distances`` rows by canonical pair key."""
    out: Dict[Tuple[str, str], Any] = {}
    for row in rows:
        out[pair_key(row["a"], row["b"])] = row.get(value)
    return out


def load_genomes_tsv(outdir: Path) -> Tuple[Table, List[str]]:
    path = Path(outdir) / "genomes.tsv"
    problems: List[str] = []
    if not path.exists():
        return Table(path, [], [], [f"missing {path.name}"]), [f"missing {path.name}"]
    table = read_table(path)
    problems.extend(table.problems)
    return table, problems


def load_clusters_tsv(outdir: Path) -> Tuple[Table, List[str]]:
    path = Path(outdir) / "clusters.tsv"
    problems: List[str] = []
    if not path.exists():
        return Table(path, [], [], [f"missing {path.name}"]), [f"missing {path.name}"]
    table = read_table(path)
    problems.extend(table.problems)
    return table, problems


def load_variants_tsv(outdir: Path, cluster: Optional[str] = None) -> Tuple[Table, List[str]]:
    outdir = Path(outdir)
    if cluster is not None:
        path = outdir / f"cluster_{cluster}.variants.tsv"
    else:
        found = sorted(outdir.glob("cluster_*.variants.tsv"))
        if not found:
            return Table(outdir / "cluster_.variants.tsv", [], []), [
                f"no cluster_*.variants.tsv in {outdir}"
            ]
        path = found[0]
    if not path.exists():
        return Table(path, [], []), [f"missing {path.name}"]
    table = read_table(path)
    return table, list(table.problems)


CONTIG_LENGTH_IN_NAME = re.compile(r"_len(\d+)_")

VARIANT_COLUMN_ALIASES = {
    "genome_name": ("genome_name", "genome", "member", "name", "sample"),
    "pivot_contig": ("pivot_contig", "ref_contig", "contig"),
    "pivot_pos": ("pivot_pos", "pivot_position", "pos", "position"),
    "ref": ("ref", "reference", "ref_allele"),
    "alt": ("alt", "alternate", "alt_allele"),
    "type": ("type", "variant_type", "kind"),
    "flags": ("flags", "flag", "annotations"),
}


@dataclass
class RecombinedDistances:
    """Pairwise SNP counts recomputed from ``cluster_*.variants.tsv``."""

    pairs: Dict[Tuple[str, str], Optional[int]]
    problems: List[str] = field(default_factory=list)
    #: Variants dropped by the filter, so a filtered number is never silent.
    variants_total: int = 0
    variants_used: int = 0
    variants_excluded: int = 0
    filter_applied: bool = False
    filter_reason: str = ""


def recompute_pairwise_snps(
    outdir: Path,
    *,
    min_pivot_contig_bp: Optional[int] = None,
    exclude_flags: Sequence[str] = (),
    snp_types: Sequence[str] = ("snp", "SNP", "sub", "substitution"),
) -> RecombinedDistances:
    """Recompute pairwise SNP counts from the variant lists, optionally filtered.

    Infinity-Win does not filter its variants: section 2 lists no filtering among its
    goals and section 15 says repeats are flagged, not corrected. Published
    reference numbers, on the other hand, are always filtered (SKA2 restricts to
    the shared k-mer core; a mapping pipeline drops high-depth and low-identity
    regions). Comparing the two unfiltered is not apples to apples, so this
    recomputes distances on a stated basis:

    * ``min_pivot_contig_bp`` drops variants on pivot contigs shorter than this.
      The length is read from the contig name (``..._len442493_cov42.0``); when a
      name carries no length the filter cannot be applied and that is reported
      rather than silently ignored.
    * ``exclude_flags`` drops variants whose ``flags`` column contains any of
      them (for example ``low_identity``, ``repeat``, ``dup``).

    Distance follows section 6.9: a member without a variant at a pivot position
    carries the pivot base (the ``ref`` column), so two members differ where their
    alleles at that position differ.
    """
    outdir = Path(outdir)
    result = RecombinedDistances(pairs={})
    if not outdir.is_dir():
        result.problems.append(f"output folder does not exist: {outdir}")
        return result
    variant_files = sorted(outdir.glob("cluster_*.variants.tsv"))
    if not variant_files:
        result.problems.append(f"no cluster_*.variants.tsv in {outdir}")
        return result

    exclude = {f.strip().lower() for f in exclude_flags if f.strip()}
    snp_set = {t.lower() for t in snp_types}
    # genome -> {(contig, pos): allele}
    profiles: Dict[str, Dict[Tuple[str, int], str]] = {}
    contig_lengths: Dict[str, Optional[int]] = {}
    unparsable_contigs: set = set()
    # Seed every cluster member, including those with no variants at all and the
    # pivot (whose profile is empty by definition, section 6.9).
    members_by_cluster = _cluster_members(outdir)
    for members in members_by_cluster.values():
        for name in members:
            profiles.setdefault(name, {})
    if min_pivot_contig_bp:
        result.filter_applied = True
    if exclude:
        result.filter_applied = True

    for path in variant_files:
        cluster_id = path.name[len("cluster_") : -len(".variants.tsv")]
        try:
            table = read_table(path)
        except PrereqError as exc:
            result.problems.append(str(exc))
            continue
        cols = {f: table.find_column(*a) for f, a in VARIANT_COLUMN_ALIASES.items()}
        if cols["genome_name"] is None or cols["pivot_pos"] is None:
            result.problems.append(
                f"{path.name}: cannot identify the genome/position columns "
                f"(header {table.header}); file skipped"
            )
            continue
        # The pivot of this cluster has an empty profile but must still appear as
        # a member, otherwise its pairs are missing from the result.
        pivot_name = _pivot_name_for_cluster(outdir, cluster_id)
        if pivot_name:
            profiles.setdefault(pivot_name, {})
        for row in table.rows:
            result.variants_total += 1
            gname = (row.get(cols["genome_name"]) or "").strip()
            if not gname:
                continue
            contig = (row.get(cols["pivot_contig"]) or "").strip() if cols["pivot_contig"] else ""
            vtype = (row.get(cols["type"]) or "").strip().lower() if cols["type"] else "snp"
            flags_raw = (row.get(cols["flags"]) or "").strip().lower() if cols["flags"] else ""
            if vtype and vtype not in snp_set:
                result.variants_excluded += 1
                continue
            if exclude:
                flag_set = {f for f in re.split(r"[,;|\s]+", flags_raw) if f}
                if flag_set & exclude:
                    result.variants_excluded += 1
                    continue
            if min_pivot_contig_bp and contig:
                if contig not in contig_lengths:
                    match = CONTIG_LENGTH_IN_NAME.search(contig)
                    contig_lengths[contig] = int(match.group(1)) if match else None
                    if contig_lengths[contig] is None:
                        unparsable_contigs.add(contig)
                length = contig_lengths[contig]
                if length is not None and length < min_pivot_contig_bp:
                    result.variants_excluded += 1
                    continue
            pos_raw = (row.get(cols["pivot_pos"]) or "").strip()
            try:
                pos = int(float(pos_raw))
            except ValueError:
                result.variants_excluded += 1
                continue
            allele = (row.get(cols["alt"]) or "").strip() if cols["alt"] else ""
            ref = (row.get(cols["ref"]) or "").strip() if cols["ref"] else ""
            key = (contig, pos)
            profile = profiles.setdefault(gname, {})
            profile[key] = allele or ref or "?"
            profile.setdefault(("__ref__", contig, pos), ref)  # type: ignore[arg-type]
            result.variants_used += 1

    if unparsable_contigs and min_pivot_contig_bp:
        result.filter_reason = (
            f"{len(unparsable_contigs)} pivot contig name(s) carry no _len<N>_ field, "
            "so the contig-length filter could not be applied to their variants "
            f"(e.g. {sorted(unparsable_contigs)[0]}); those variants were KEPT"
        )
        result.problems.append(result.filter_reason)

    # ref allele lookup, so a member without a variant gets the pivot base
    ref_alleles: Dict[Tuple[str, int], str] = {}
    for profile in profiles.values():
        for key, value in list(profile.items()):
            if isinstance(key, tuple) and len(key) == 3 and key[0] == "__ref__":
                ref_alleles[(key[1], key[2])] = value
                del profile[key]  # type: ignore[arg-type]

    names = sorted(profiles)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            pa, pb = profiles[a], profiles[b]
            positions = set(pa) | set(pb)
            snps = 0
            for key in positions:
                ref = ref_alleles.get(key, "")
                allele_a = pa.get(key, ref)
                allele_b = pb.get(key, ref)
                if allele_a != allele_b:
                    snps += 1
            result.pairs[pair_key(a, b)] = snps
    if not result.pairs:
        result.problems.append("no variant rows could be read, so no pair was recomputed")
    return result


def _cluster_members(outdir: Path) -> Dict[str, List[str]]:
    """Every genome name per cluster id, from ``genomes.tsv``.

    Needed because a member identical to its pivot contributes no variant rows at
    all: without this it would be missing from the recomputed distances, and a
    true 0-SNP pair would look like "never compared" instead of "identical".
    """
    path = outdir / "genomes.tsv"
    out: Dict[str, List[str]] = {}
    if not path.exists():
        return out
    try:
        table = read_table(path)
    except PrereqError:
        return out
    name_col = table.find_column("name", "genome_name", "genome", "sample")
    cluster_col = table.find_column("cluster", "cluster_id")
    if name_col is None or cluster_col is None:
        return out
    for row in table.rows:
        name = (row.get(name_col) or "").strip()
        cluster = (row.get(cluster_col) or "").strip()
        if name and cluster:
            out.setdefault(cluster, []).append(name)
    return out


def _pivot_name_for_cluster(outdir: Path, cluster_id: str) -> Optional[str]:
    """The pivot genome name of one cluster, from ``clusters.tsv``."""
    path = outdir / "clusters.tsv"
    if not path.exists():
        return None
    try:
        table = read_table(path)
    except PrereqError:
        return None
    id_col = table.find_column("cluster_id", "cluster", "id")
    name_col = table.find_column("pivot_name", "pivot")
    if id_col is None or name_col is None:
        return None
    for row in table.rows:
        if (row.get(id_col) or "").strip() == str(cluster_id):
            value = (row.get(name_col) or "").strip()
            return value or None
    return None


def comparable_bp_for(outdir: Path, run_json: Optional[RunJson] = None) -> Optional[int]:
    """A comparable-site count to convert "half a SNP" into branch-length units.

    Prefers the largest ``comparable_bp`` seen in ``cluster_*.dist.tsv``, and
    falls back to the mean genome length (``run.json`` total_length / genomes).
    """
    rows, _ = load_cluster_distances(Path(outdir))
    values = [r["comparable_bp"] for r in rows if r.get("comparable_bp")]
    if values:
        return max(values)
    if run_json is None:
        run_json = load_run_json(Path(outdir), required=False)
    if run_json is not None and run_json.total_length:
        # total_length is summed over every genome; half a SNP is a fraction of
        # ONE genome's comparable sites, so use the mean genome length.
        genomes = run_json.genomes
        return int(run_json.total_length / genomes) if genomes else run_json.total_length
    return None


# ---------------------------------------------------------------------------
# Dataset descriptors written by real_datasets.py
# ---------------------------------------------------------------------------


def load_dataset(path: Path) -> Dict[str, Any]:
    """Read a ``dataset.json`` and check the assemblies are still there."""
    obj = load_json(path)
    if not isinstance(obj, dict):
        raise PrereqError(f"{path}: expected a JSON object")
    for key in ("name", "asm_dir"):
        if not obj.get(key):
            raise PrereqError(f"{path}: missing required key {key!r}")
    return obj


def dataset_asm_files(dataset: Dict[str, Any]) -> List[Path]:
    asm = Path(dataset["asm_dir"])
    if not asm.is_dir():
        raise PrereqError(
            f"dataset {dataset.get('name')!r}: assembly folder {asm} does not exist; "
            "re-run scripts/real_datasets.py"
        )
    files: List[Path] = []
    for pattern in ("*.fasta", "*.fa", "*.fna", "*.fasta.gz", "*.fa.gz", "*.fna.gz"):
        files.extend(sorted(asm.glob(pattern)))
    seen = set()
    out = []
    for f in files:
        if f not in seen:
            seen.add(f)
            out.append(f)
    if not out:
        raise PrereqError(f"dataset {dataset.get('name')!r}: no FASTA files in {asm}")
    return out


def find_datasets(work_dir: Path, names: Optional[Sequence[str]] = None) -> List[Dict[str, Any]]:
    """Every ``bench-work/<name>/dataset.json``, optionally filtered."""
    work_dir = Path(work_dir)
    if not work_dir.is_dir():
        raise PrereqError(
            f"benchmark work directory {work_dir} does not exist; "
            "run scripts/real_datasets.py first"
        )
    found = sorted(work_dir.glob("*/dataset.json"))
    out = []
    for path in found:
        obj, problem = load_json_optional(path)
        if problem is not None or not isinstance(obj, dict):
            print(f"WARNING: skipping {path}: {problem}", file=sys.stderr)
            continue
        obj["_dataset_json"] = str(path)
        obj["_dir"] = str(path.parent)
        out.append(obj)
    if names:
        wanted = set(names)
        selected = [d for d in out if d.get("name") in wanted]
        missing = wanted - {d.get("name") for d in selected}
        if missing:
            raise PrereqError(
                f"no dataset.json for {sorted(missing)} under {work_dir}; "
                f"available: {sorted(d.get('name') for d in out)}"
            )
        return selected
    if not out:
        raise PrereqError(
            f"no dataset.json found under {work_dir}; "
            "run scripts/real_datasets.py --datasets all first"
        )
    return out


def dataset_truth(
    work_dir: Optional[Path],
    name: str,
    summary: Optional[Dict[str, Any]] = None,
) -> Tuple[Optional[Path], Optional[Path], str]:
    """``(reference tree, truth folder, where they came from)`` for one dataset.

    Looks at ``<work>/<name>/dataset.json`` first, then at the ``datasets`` list
    of a matrix ``summary.json``: a simulated dataset can live outside the work
    folder (on tmpfs, say), and then the summary is the only place that still
    knows where its truth is. Paths that do not exist are dropped.
    """
    candidates: List[Tuple[Dict[str, Any], str]] = []
    if work_dir is not None:
        path = Path(work_dir) / name / "dataset.json"
        obj, _ = load_json_optional(path)
        if isinstance(obj, dict):
            candidates.append((obj, str(path)))
    for entry in (summary or {}).get("datasets") or []:
        if isinstance(entry, dict) and str(entry.get("name")) == str(name):
            candidates.append((entry, "matrix summary.json"))
    tree: Optional[Path] = None
    tdir: Optional[Path] = None
    origin = "no dataset descriptor found"
    for obj, where in candidates:
        truth = obj.get("truth") if isinstance(obj.get("truth"), dict) else {}
        t, d = truth.get("tree"), truth.get("dir")
        if tree is None and t and Path(t).exists():
            tree, origin = Path(t), where
        if tdir is None and d and Path(d).is_dir():
            tdir, origin = Path(d), where
        if tree is None and tdir is None:
            origin = f"{where} declares no existing truth"
    return tree, tdir, origin


# ---------------------------------------------------------------------------
# Optional plotting
# ---------------------------------------------------------------------------


def try_import_matplotlib() -> Tuple[Any, Optional[str]]:
    """Return ``(matplotlib.pyplot, None)`` or ``(None, reason)``.

    Figures are a nice-to-have; the numbers are the deliverable. Every caller
    writes its CSV/TSV first and only then tries to plot.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt, None
    except Exception as exc:  # noqa: BLE001 - any import failure is a "no plot"
        return None, f"matplotlib unavailable ({type(exc).__name__}: {exc})"


# ---------------------------------------------------------------------------
# Configuration matrix (PROJECT.md section 11)
# ---------------------------------------------------------------------------


@dataclass
class Config:
    """One cell of the configuration axis of the experiment."""

    name: str
    flags: List[str]
    #: Which part of the idea this configuration turns off, in the wording of
    #: section 11. ``None`` for the baseline (everything on).
    removes: Optional[str] = None
    #: The ``run.json`` ``params`` field that must be false (or true, for
    #: ``single_cluster``) when the switch took effect.
    param: Optional[str] = None
    param_expected: Optional[bool] = None

    @property
    def is_baseline(self) -> bool:
        return not self.flags


#: The six configurations of section 11: baseline plus one run per switch.
CONFIGS: Tuple[Config, ...] = (
    Config("baseline", [], removes=None),
    Config(
        "no_block_skip",
        ["--no-block-skip"],
        removes="block hashing altogether; every stretch between shared anchors is compared directly",
        param="block_skip",
        param_expected=False,
    ),
    Config(
        "flat",
        ["--flat"],
        removes="nested levels; only leaf blocks are used",
        param="hierarchical",
        param_expected=False,
    ),
    Config(
        "no_piece_cache",
        ["--no-piece-cache"],
        removes="resolving a shared difference piece once",
        param="piece_cache",
        param_expected=False,
    ),
    Config(
        "no_profile_grouping",
        ["--no-profile-grouping"],
        removes="handling identical genomes once",
        param="profile_grouping",
        param_expected=False,
    ),
    Config(
        "single_cluster",
        ["--single-cluster"],
        removes="cluster-first; all genomes are compared to one pivot",
        param="single_cluster",
        param_expected=True,
    ),
)

CONFIG_BY_NAME: Dict[str, Config] = {c.name: c for c in CONFIGS}

#: Switch names in the order the report presents them.
SWITCHES: Tuple[str, ...] = tuple(c.name for c in CONFIGS if not c.is_baseline)


def select_configs(names: Optional[Sequence[str]]) -> List[Config]:
    if not names:
        return list(CONFIGS)
    out: List[Config] = []
    for name in names:
        if name not in CONFIG_BY_NAME:
            raise PrereqError(
                f"unknown configuration {name!r}; available: {list(CONFIG_BY_NAME)}"
            )
        out.append(CONFIG_BY_NAME[name])
    if not any(c.is_baseline for c in out):
        # The verdict rule compares every switch against the baseline, so the
        # baseline is always part of a matrix run.
        out.insert(0, CONFIG_BY_NAME["baseline"])
    return out


# ---------------------------------------------------------------------------
# The verdict rule (PROJECT.md section 11)
# ---------------------------------------------------------------------------

#: "at least 20% worse"
VERDICT_THRESHOLD = 1.20
#: "on at least two of the three datasets"
VERDICT_MIN_DATASETS = 2

#: Metrics the rule looks at, in order. Section 11: "turning it off makes
#: time-to-exact or peak memory at least 20% worse".
VERDICT_METRICS = (
    ("exact_tree_s", "time to exact tree", "higher_is_worse"),
    ("peak_memory_bytes", "peak memory", "higher_is_worse"),
)


@dataclass
class VerdictMetric:
    dataset: str
    metric: str
    baseline: Optional[float]
    switched_off: Optional[float]
    ratio: Optional[float]
    worse: Optional[bool]
    note: str = ""

    def as_dict(self) -> Dict[str, Any]:
        return {
            "dataset": self.dataset,
            "metric": self.metric,
            "baseline": self.baseline,
            "switched_off": self.switched_off,
            "ratio": self.ratio,
            "worse": self.worse,
            "note": self.note,
        }


@dataclass
class Verdict:
    """The judgement for one switch, across every dataset."""

    switch: str
    removes: str
    metrics: List[VerdictMetric] = field(default_factory=list)
    pays: Optional[bool] = None
    #: One of ``pays``, ``does_not_pay``, ``undetermined`` (too few datasets were
    #: measurable to apply the rule), ``not_measurable`` (no metric at all).
    #: ``pays is None`` alone cannot tell the last two apart, and the report must.
    state: str = "not_measurable"
    #: Datasets where at least one verdict metric got >= 20% worse.
    datasets_worse: List[str] = field(default_factory=list)
    datasets_measurable: List[str] = field(default_factory=list)
    datasets_total: int = 0
    reason: str = ""
    threshold: float = VERDICT_THRESHOLD
    min_datasets: int = VERDICT_MIN_DATASETS

    def as_dict(self) -> Dict[str, Any]:
        return {
            "switch": self.switch,
            "removes": self.removes,
            "pays": self.pays,
            "state": self.state,
            "datasets_worse": list(self.datasets_worse),
            "datasets_measurable": list(self.datasets_measurable),
            "datasets_total": self.datasets_total,
            "reason": self.reason,
            "threshold": self.threshold,
            "min_datasets": self.min_datasets,
            "metrics": [m.as_dict() for m in self.metrics],
        }


def evaluate_verdict(
    switch: str,
    removes: str,
    baseline_by_dataset: Dict[str, Dict[str, Any]],
    switched_by_dataset: Dict[str, Dict[str, Any]],
    dataset_names: Optional[Sequence[str]] = None,
) -> Verdict:
    """Apply the section 11 rule to one switch.

    ``baseline_by_dataset`` / ``switched_by_dataset`` map a dataset name to a
    measurement dict that may contain ``exact_tree_s`` and
    ``peak_memory_bytes`` (either may be ``None``). ``dataset_names`` is the full
    set the rule should consider, so a dataset that was never run is named in the
    verdict instead of quietly shrinking the denominator.

    A part of the idea *pays* if turning it off makes time-to-exact **or** peak
    memory at least 20% worse on at least two of the three datasets. A dataset
    where neither metric is measurable counts as neither worse nor not-worse:
    the verdict then says it could not be decided, rather than reporting a win
    (section 2 anti-goal).
    """
    verdict = Verdict(switch=switch, removes=removes or "")
    considered = list(dataset_names) if dataset_names else []
    names = sorted(set(baseline_by_dataset) | set(switched_by_dataset) | set(considered))
    verdict.datasets_total = len(names)
    for name in names:
        base = baseline_by_dataset.get(name) or {}
        off = switched_by_dataset.get(name) or {}
        any_measurable = False
        any_worse = False
        for key, label, _ in VERDICT_METRICS:
            b = base.get(key)
            o = off.get(key)
            r = ratio(o, b)
            worse: Optional[bool] = None
            note = ""
            if r is None:
                if not base and not off:
                    note = "this dataset was never run"
                elif b is None and o is None:
                    note = "not measured in either run"
                elif b in (None, 0):
                    note = "baseline value missing or zero, ratio undefined"
                else:
                    note = "switched-off value missing"
            else:
                any_measurable = True
                worse = r >= VERDICT_THRESHOLD
                if worse:
                    any_worse = True
            verdict.metrics.append(
                VerdictMetric(
                    dataset=name,
                    metric=label,
                    baseline=None if b is None else float(b),
                    switched_off=None if o is None else float(o),
                    ratio=r,
                    worse=worse,
                    note=note,
                )
            )
        if any_measurable:
            verdict.datasets_measurable.append(name)
        if any_worse:
            verdict.datasets_worse.append(name)
    n_worse = len(verdict.datasets_worse)
    n_measurable = len(verdict.datasets_measurable)
    pct = int(round((VERDICT_THRESHOLD - 1) * 100))
    unmeasurable = [n for n in names if n not in verdict.datasets_measurable]
    if n_measurable == 0:
        verdict.pays = None
        verdict.state = "not_measurable"
        verdict.reason = (
            "not measurable: neither time-to-exact nor peak memory was available for "
            f"any dataset ({', '.join(unmeasurable) or 'none were run'}), so the rule "
            "cannot be applied"
        )
    elif n_measurable < VERDICT_MIN_DATASETS:
        # The rule needs "at least two of the three datasets". With fewer than two
        # measurable datasets it can decide neither way, and saying "does not pay"
        # would be a conclusion the data do not support.
        verdict.pays = None
        verdict.state = "undetermined"
        verdict.reason = (
            f"undetermined: only {n_measurable} of {verdict.datasets_total} datasets "
            f"were measurable ({', '.join(verdict.datasets_measurable)}), fewer than the "
            f"{VERDICT_MIN_DATASETS} the rule needs"
        )
        if unmeasurable:
            verdict.reason += f"; no usable numbers for {', '.join(unmeasurable)}"
    elif n_worse >= VERDICT_MIN_DATASETS:
        verdict.pays = True
        verdict.state = "pays"
        verdict.reason = (
            f"pays: turning it off made time-to-exact or peak memory >={pct}% worse "
            f"on {n_worse} of {verdict.datasets_total} datasets "
            f"({', '.join(verdict.datasets_worse)})"
        )
        if unmeasurable:
            verdict.reason += f"; not measurable on {', '.join(unmeasurable)}"
    else:
        verdict.pays = False
        verdict.state = "does_not_pay"
        detail = (
            f"only {n_worse} of {verdict.datasets_total} datasets got >={pct}% worse"
        )
        if unmeasurable:
            detail += f"; not measurable on {len(unmeasurable)} ({', '.join(unmeasurable)})"
        verdict.reason = f"does not pay by the section 11 rule: {detail}"
    return verdict


# ---------------------------------------------------------------------------
# Expected-value acceptance checks
# ---------------------------------------------------------------------------


@dataclass
class ExpectedCheck:
    """One comparison between an expected value and what a run produced."""

    name: str
    pair: Optional[Tuple[str, str]]
    expected: Optional[float]
    observed: Optional[float]
    tolerance: float
    unit: str
    passed: Optional[bool]
    source: str = ""
    note: str = ""

    def as_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "pair": list(self.pair) if self.pair else None,
            "expected": self.expected,
            "observed": self.observed,
            "tolerance": self.tolerance,
            "unit": self.unit,
            "passed": self.passed,
            "source": self.source,
            "note": self.note,
        }

    @property
    def status(self) -> str:
        if self.passed is None:
            return "not measured"
        return "PASS" if self.passed else "FAIL"


def check_expected_pairs(
    expected: Dict[str, Any], observed: Dict[Tuple[str, str], Optional[float]]
) -> List[ExpectedCheck]:
    """Compare published SNP distances against what Infinity-Win reported.

    ``observed`` keys are canonicalised, so the caller may index pairs in either
    order.
    """
    index: Dict[Tuple[str, str], Optional[float]] = {
        pair_key(*key): value for key, value in (observed or {}).items()
    }
    checks: List[ExpectedCheck] = []
    default_tolerance = float(expected.get("abs_tolerance_snps", 0))
    for entry in expected.get("pairs", []):
        a = entry.get("genome_a") or entry.get("a")
        b = entry.get("genome_b") or entry.get("b")
        value = entry.get("snps", entry.get("value"))
        tolerance = float(entry.get("abs_tolerance_snps", default_tolerance))
        unit = entry.get("unit", "SNPs")
        key = pair_key(a, b) if a and b else None
        obs = index.get(key) if key else None
        if obs is None:
            passed = None
            note = "Infinity-Win reported no distance for this pair (not in any cluster_*.dist.tsv)"
        else:
            try:
                passed = abs(float(obs) - float(value)) <= tolerance
            except (TypeError, ValueError):
                passed = None
            note = "" if passed else f"expected {value} +/- {tolerance}, observed {obs}"
        checks.append(
            ExpectedCheck(
                name=entry.get("name") or f"{a}-{b}",
                pair=key,
                expected=None if value is None else float(value),
                observed=None if obs is None else float(obs),
                tolerance=tolerance,
                unit=unit,
                passed=passed,
                source=entry.get("source", expected.get("source", "")),
                note=note or entry.get("note", ""),
            )
        )
    return checks


# ---------------------------------------------------------------------------
# bench.json reader
# ---------------------------------------------------------------------------


@dataclass
class BenchJson:
    """A ``infinity bench`` result, read defensively.

    The schema the CLI actually writes is::

        { "output", "truth",
          "run":      { ...the whole run.json... },
          "distance": { pairs_in_truth, pairs_compared, pairs_missing_from_output,
                        pairs_exact, max_abs_error, sum_abs_error, sum_true_snps,
                        sum_reported_snps, worst_pair: [a, b, true, reported] },
          "clusters": { genomes_in_truth, genomes_compared, true_clusters,
                        reported_clusters, recovered_exactly, pairwise_f1 },
          "trees":    [ {source, truth, leaves, rf, splits_emitted, splits_true,
                         rf_normalised, t_ms, approximate, focus, units} ],
          "timeline": [ ...same shape, one entry per emitted tree... ],
          "problems": [] }

    Note ``rf_normalised`` (British spelling) and that there is **no single
    full-tree RF**: the engine scores each cluster tree against its own truth
    tree. So the aggregate exposed here is the sum of the per-cluster RF values,
    and the mean of the per-cluster normalised values - both labelled as such by
    the report rather than presented as one tree's distance.
    """

    path: Path
    raw: Dict[str, Any]
    problems: List[str] = field(default_factory=list)

    def get(self, dotted: str, default: Any = None) -> Any:
        value = dig(self.raw, dotted, None)
        if value is None:
            self._note(f"{dotted} missing from {self.path.name}")
            return default
        return value

    def _note(self, message: str) -> None:
        if message not in self.problems:
            self.problems.append(message)

    # -- trees -----------------------------------------------------------
    def tree_entries(self) -> List[Dict[str, Any]]:
        entries = self.raw.get("trees")
        if isinstance(entries, list):
            return [e for e in entries if isinstance(e, dict)]
        self._note("no trees[] in bench output: no per-cluster RF was reported")
        return []

    @property
    def rf(self) -> Optional[float]:
        """Sum of the per-cluster RF distances (``trees[].rf``)."""
        entries = self.tree_entries()
        values = [_to_float(e.get("rf")) for e in entries]
        values = [v for v in values if v is not None]
        if values:
            return float(sum(values))
        for dotted in ("tree.rf", "rf", "tree.rf_distance", "accuracy.rf"):
            value = dig(self.raw, dotted, None)
            if value is not None:
                return _to_float(value)
        self._note("no RF value anywhere in the bench output")
        return None

    @property
    def rf_normalized(self) -> Optional[float]:
        """Mean of the per-cluster normalised RF values that are defined."""
        entries = self.tree_entries()
        values = [
            v
            for v in (
                _to_float(e.get("rf_normalised", e.get("rf_normalized"))) for e in entries
            )
            if v is not None
        ]
        if values:
            return sum(values) / len(values)
        for dotted in ("tree.rf_normalized", "tree.rf_normalised", "rf_normalized"):
            value = dig(self.raw, dotted, None)
            if value is not None:
                return _to_float(value)
        return None

    @property
    def n_trees_scored(self) -> Optional[int]:
        entries = self.tree_entries()
        return len(entries) if entries else None

    def rf_timeline(self) -> List[Dict[str, Any]]:
        """Per-emitted-tree RF values: the accuracy-over-time series."""
        for key in ("timeline", "rf_timeline", "tree.timeline", "trees"):
            value = self.raw.get(key)
            if isinstance(value, list) and value:
                return [v for v in value if isinstance(v, dict)]
        return []

    # -- distances -------------------------------------------------------
    @property
    def distance_mae(self) -> Optional[float]:
        """Mean absolute SNP error per compared pair."""
        total = _to_float(dig(self.raw, "distance.sum_abs_error", None))
        pairs = _to_float(dig(self.raw, "distance.pairs_compared", None))
        if total is not None and pairs:
            return total / pairs
        for dotted in (
            "distances.mean_abs_error_snps",
            "distances.mae_snps",
            "distance.mean_abs_error",
        ):
            value = dig(self.raw, dotted, None)
            if value is not None:
                return _to_float(value)
        return None

    @property
    def distance_max_error(self) -> Optional[float]:
        for dotted in ("distance.max_abs_error", "distances.max_abs_error_snps"):
            value = dig(self.raw, dotted, None)
            if value is not None:
                return _to_float(value)
        return None

    @property
    def pairs_compared(self) -> Optional[int]:
        return _to_int(
            dig(self.raw, "distance.pairs_compared", None)
            if dig(self.raw, "distance.pairs_compared", None) is not None
            else dig(self.raw, "distances.pairs_compared", None)
        )

    @property
    def pairs_exact(self) -> Optional[int]:
        return _to_int(
            dig(self.raw, "distance.pairs_exact", None)
            if dig(self.raw, "distance.pairs_exact", None) is not None
            else dig(self.raw, "distances.pairs_exact", None)
        )

    @property
    def pairs_in_truth(self) -> Optional[int]:
        return _to_int(dig(self.raw, "distance.pairs_in_truth", None))

    @property
    def pairs_missing_from_output(self) -> Optional[int]:
        return _to_int(dig(self.raw, "distance.pairs_missing_from_output", None))

    @property
    def worst_pair(self) -> Optional[List[Any]]:
        value = dig(self.raw, "distance.worst_pair", None)
        return list(value) if isinstance(value, (list, tuple)) else None

    @property
    def fraction_exact(self) -> Optional[float]:
        exact, compared = self.pairs_exact, self.pairs_compared
        if exact is None or not compared:
            return None
        return exact / compared

    # -- clusters --------------------------------------------------------
    @property
    def cluster_f1(self) -> Optional[float]:
        return _to_float(dig(self.raw, "clusters.pairwise_f1", None))

    @property
    def clusters_recovered(self) -> Optional[int]:
        return _to_int(dig(self.raw, "clusters.recovered_exactly", None))

    @property
    def clusters_true(self) -> Optional[int]:
        return _to_int(dig(self.raw, "clusters.true_clusters", None))

    # -- the embedded run ------------------------------------------------
    @property
    def run_summary(self) -> Optional[Dict[str, Any]]:
        """``bench.json`` embeds the whole run summary, so it is self-contained."""
        value = self.raw.get("run")
        return value if isinstance(value, dict) else None

    @property
    def bench_problems(self) -> List[str]:
        value = self.raw.get("problems")
        return [str(v) for v in value] if isinstance(value, list) else []


def load_bench_json(path: Path, *, required: bool = True) -> Optional[BenchJson]:
    p = Path(path)
    if p.is_dir():
        p = p / "bench.json"
    if not p.exists():
        if required:
            raise PrereqError(
                f"no bench result at {p}\n"
                "       run 'infinity bench --output <run dir> --truth <truth dir>' for this "
                "output folder first (run_matrix.py does this automatically when the "
                "dataset has a simulator truth folder)"
            )
        return None
    obj, problem = load_json_optional(p)
    if problem is not None or not isinstance(obj, dict):
        if required:
            raise PrereqError(problem or f"{p}: expected a JSON object")
        return None
    return BenchJson(path=p, raw=obj)


def rf_values_from_bench(bench: BenchJson) -> Dict[str, List[Dict[str, Any]]]:
    """Extract every RF-over-time series a bench JSON carries.

    Prefers the engine's own ``timeline`` (one entry per emitted tree), then the
    final per-cluster ``trees``, then a single aggregate ``tree`` block.
    """
    out: Dict[str, List[Dict[str, Any]]] = {}
    timeline = bench.rf_timeline()
    if timeline:
        out["timeline"] = timeline
    trees = bench.tree_entries()
    if trees and trees is not timeline:
        out["trees"] = trees
    if not out:
        rf = bench.rf
        if rf is not None:
            out["aggregate"] = [
                {
                    "t_ms": dig(bench.raw, "run.timings.exact_tree_ms", None),
                    "rf": rf,
                    "rf_normalised": bench.rf_normalized,
                    "source": "bench.json aggregate",
                }
            ]
    return out


def is_simulator_truth_dir(path: Optional[Path]) -> bool:
    """True when ``path`` looks like a ``infinity simulate`` truth folder.

    ``infinity bench --truth`` requires that layout (``tree_*.nwk``,
    ``pairwise_snps.tsv``, ``clusters.tsv``). A real dataset's reference folder --
    a SKA2 distance table and an IQ-TREE tree -- is *not* that, and passing it
    would produce a bench failure that looks like an Infinity-Win bug.
    """
    if not path:
        return False
    p = Path(path)
    if not p.is_dir():
        return False
    if list(p.glob("tree_*.nwk")):
        return True
    return (p / "pairwise_snps.tsv").exists() and (p / "clusters.tsv").exists()
