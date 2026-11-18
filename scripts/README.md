# Infinity-Win benchmark and reporting scripts

This directory is the measurement half of PROJECT.md section 11: *"the project is
an experiment as much as a tool, so measuring it honestly is part of the
deliverable."*

Everything here is **Python 3 standard library only**. `matplotlib` is used for
figures when it is importable; when it is not, every number is still written as
CSV/TSV, so a missing plotting library can never cost a measurement.

Three rules the scripts are written to obey:

| Rule | Where it comes from | How it shows up |
|---|---|---|
| A benchmark that only reports the cases where the idea wins is a failure | section 2 (anti-goal) | Failed configurations are written to `summary.json` with `status: "failed"`, their exit code and their stderr, and they are printed in `REPORT.md` section 8. Nothing is dropped. |
| Never present an approximate number as exact | section 0, rule 8 | A missing timing, counter or memory value prints as `n/a`, never `0`. Peak memory sampled from `/proc` is flagged `is_floor: true`. A wall-time substitute for `exact_tree_ms` is labelled "proxy". |
| The engine must never depend on the GUI; everything reachable from the CLI | section 0, rule 2 | These scripts drive `infinity run` / `infinity bench` only. The GUI is never involved. |

---

## 1. Prerequisites

```bash
export PATH="$HOME/.cargo/bin:$PATH"
cargo build --release          # -> /tmp/infinity-target/release/infinity   (see .cargo/config.toml)
python3 --version              # 3.8+ ; matplotlib optional
```

Scripts find the binary in this order: `--infinity <path>` → `$INFINITY_BIN` →
`$CARGO_TARGET_DIR/release/infinity` → `/tmp/infinity-target/release/infinity` →
`<repo>/target/release/infinity` → `<repo>/target/debug/infinity` → `infinity` on `PATH`.
If none exists you get one actionable message naming every path that was tried —
not a traceback.

Generated data lives under **`bench-work/`** (gitignored). The root partition of
the development machine is 99% full, so nothing is written anywhere else except
`REPORT.md` at the repository root and `dist/` for release binaries.

---

## 2. Run everything, in order

```bash
# 0. See what exists: datasets, sources, tools, binaries.
scripts/real_datasets.py --list
scripts/baselines.py status

# 1. Assemble the real datasets (copies + verifies; never moves or edits sources).
scripts/real_datasets.py --datasets all

# 2. Look at the exact command lines the matrix will run. Writes nothing.
scripts/run_matrix.py --dry-run

# 3. Run the matrix: 6 configurations x every dataset. Sequential on purpose, so
#    one configuration's wall time and peak memory are not contaminated by another.
scripts/run_matrix.py

# 4. Accuracy-over-time curves (run_matrix's outputs + events.jsonl).
scripts/accuracy_over_time.py --matrix-dir bench-work/matrix

# 5. Outside baselines, joined against one Infinity-Win run.
scripts/baselines.py all --dataset ncbi_pd_cluster \
    --run-dir bench-work/matrix/ncbi_pd_cluster/baseline

# 6. The report. Writes REPORT.md at the repository root.
scripts/make_report.py

# 7. Tests for all of the above (no infinity binary needed).
python3 -m unittest discover -s scripts/tests -t scripts -v

# 8. Release binaries.
scripts/release.sh --check     # what would happen
scripts/release.sh             # Linux x86-64 into ./dist
```

Re-running any step is safe. `run_matrix.py` skips a configuration whose output
folder already has `run.json` (marking it `cached`), and it **merges** with an
existing `summary.json`, so re-running one dataset never erases another
dataset's results — including its failures. Use `--force` / `--force-bench` to
redo work.

---

## 3. The scripts

### `real_datasets.py` — assemble the three real datasets

Copies source assemblies into `bench-work/<dataset>/asm/` and writes a
`dataset.json` recording provenance: source path, byte size and SHA-256 of both
source and copy, licence notes, caveats, and what truth (if any) exists.
Every copy is verified by checksum; the script proves it copied rather than
asserting that it did.

| Dataset | Source | Genomes | Truth |
|---|---|---|---|
| `salmonella_trio` | `.../analysis/asm/{isolate-A,isolate-B,isolate-C}.fasta` | 3 | published SNP distances in `scripts/expected/salmonella_trio.json`, plus a 3-row SKA2 extract |
| `ncbi_pd_cluster` | `.../analysis/ncbi_pd/asm/*.fasta` (PDS000037076) | 32 | SKA2 distance table (496 rows) + IQ-TREE tree, both copied into `truth/` |
| `s_aureus_cc` | `MORI/refstrains/closed/*/*.fna`, `MORI/refstrains/genomes/*.fna` | sampled | clonal-complex labels only (weak) |

```bash
scripts/real_datasets.py --list
scripts/real_datasets.py --datasets salmonella_trio ncbi_pd_cluster
scripts/real_datasets.py --datasets s_aureus_cc --max-per-cc 2 --max-extra-genomes 6
scripts/real_datasets.py --datasets all --clean          # rebuild from scratch
```

**`s_aureus_cc` is sampled, not copied wholesale**, because the MORI collection is
large and the disk is nearly full. The sample is deterministic (sorted order):
`--max-per-cc` closed genomes per clonal complex (default 2) and
`--max-extra-genomes` classic reference strains (default 6). `--all` takes
everything.

> **Caveat recorded in `dataset.json` and repeated in the report:** the
> *S. aureus* genomes are **different clonal complexes, i.e. far apart**. That
> makes them useful for checking that clustering separates lineages and useless
> for checking SNP-level accuracy — at several percent divergence Infinity-Win's block
> shortcut stops helping (PROJECT.md section 15, "Divergent clusters"). There is
> no ground-truth tree or distance for this set, so `infinity bench` cannot score
> accuracy on it.

### `run_matrix.py` — the configuration x dataset matrix

Six configurations, exactly as section 11 lists them: `baseline` (everything on)
plus one run per switch — `no_block_skip`, `flat`, `no_piece_cache`,
`no_profile_grouping`, `single_cluster`. Selecting any switch always drags the
baseline along, because the verdict rule compares against it.

```bash
scripts/run_matrix.py --dry-run                       # print command lines only
scripts/run_matrix.py --datasets ncbi_pd_cluster
scripts/run_matrix.py --configs baseline flat --force
scripts/run_matrix.py --threads 8 --k 31 --densities 16384,2048,256
scripts/run_matrix.py --extra-run-arg --verbose       # repeatable pass-through
```

Per cell it records: the exact `run` and `bench` command lines, wall time,
return code, `run.json` timings/counters/`focus_counts`/`peak_memory_bytes`,
whether the switch actually took effect, and the bench result.

**Switch verification.** A flag that the CLI accepts and then ignores would
silently invalidate the whole experiment, so each cell's `run.json.params` is
checked against what the switch should have set (`--flat` ⇒ `hierarchical=false`,
`--single-cluster` ⇒ `single_cluster=true`, baseline ⇒ all optimisations on). A
mismatch is recorded as `switch_verified: false`, printed as
`SWITCH NOT APPLIED`, and reported in `REPORT.md`.

**Peak memory** prefers `run.json.peak_memory_bytes` (exact — the engine reads
its own `VmHWM`). When that is null, the script samples `/proc/<pid>/status`
`VmHWM` every 50 ms while the run is in flight. `VmHWM` never decreases, so any
sample is a true *lower bound*; such values are flagged `is_floor: true` and
shown as `(floor)` rather than presented as the peak.

Output layout:

```
bench-work/matrix/
  summary.json                       every cell, including the failed ones
  <dataset>/<config>/
    genomes.tsv clusters.tsv tree.nwk cluster_<id>.{nwk,dist.tsv,variants.tsv}
    run.json events.jsonl            written by infinity run
    bench.json                       written by infinity bench
    infinity.run.log infinity.run.err      stdout/stderr of the run
    infinity.bench.log infinity.bench.err
  figures/                           accuracy curves (CSV always, PNG if possible)
```

`summary.json` is rewritten after **every** cell, so an interrupted or crashed
matrix still leaves a usable, honest summary behind.

### `accuracy_over_time.py` — the key figure

Parses `events.jsonl` and produces tree accuracy versus elapsed time, one curve
per configuration.

Accuracy at an instant = **Robinson-Foulds distance** between the tree emitted by
then and the truth tree, *after collapsing branches shorter than half a SNP*.
The full tree at each instant is re-assembled the way section 6.9 describes: take
the `backbone_tree` (leaves are clusters) and replace each pivot leaf with the
latest `cluster_focus` subtree.

```bash
# one run
scripts/accuracy_over_time.py --run-dir bench-work/matrix/ncbi_pd_cluster/baseline \
    --truth-tree bench-work/ncbi_pd_cluster/truth/cluster.treefile

# the whole matrix
scripts/accuracy_over_time.py --matrix-dir bench-work/matrix

# read RF values from a bench JSON instead of recomputing them
scripts/accuracy_over_time.py --run-dir <dir> --rf-from <dir>/bench.json --rf-mode bench

# force the pure-Python RF even when a bench.json exists
scripts/accuracy_over_time.py --matrix-dir bench-work/matrix --rf-mode internal

# let the CLI's own RF code do the scoring (one temp output folder per tree)
scripts/accuracy_over_time.py --matrix-dir bench-work/matrix --rf-mode infinity --infinity <path>
```

`--rf-mode` decides where RF values come from:

| mode | behaviour |
|---|---|
| `auto` (default) | the best source that can actually form a curve: a `bench.json` series with RF at two or more instants, else RF recomputed here, else the bench series anyway. The choice is recorded per row in `rf_source` and explained in the notes. |
| `bench` | require RF values from `--rf-from <bench.json>` (or `<run-dir>/bench.json`); no truth tree needed |
| `internal` | always compute RF in pure Python from the emitted Newick, against per-cluster truth trees (`--truth-dir`) or one reference tree (`--truth-tree`) |
| `infinity` | re-run `infinity bench --output <run dir> --truth <truth dir>` and read its series; needs a simulator truth folder |

`--rf-from`/`--rf-mode bench` is what makes the script usable **before** the CLI's
`bench` subcommand exists, and `internal` is what makes it usable when `bench`
exists but its series has no time resolution (see section 4).

**"Half a SNP" needs a unit.** Infinity-Win writes substitutions per site inside
clusters (section 6.9) while the simulator's truth tree is in SNPs (section 10.1),
so one number cannot serve both. `--collapse-mode`:

- `auto` (default) — use `--comparable-bp` if given (or the largest
  `comparable_bp` in `cluster_*.dist.tsv`, else `run.json.total_length`), giving
  `threshold = collapse_snps / comparable_bp`; otherwise, if the truth tree's own
  branch lengths look like integer SNP counts, use `collapse_snps` directly;
  otherwise **do not collapse and say so**.
- `snps` — require a comparable-site count; error out rather than guess.
- `length X` — an absolute threshold in tree units.
- `none` — no collapsing.

An unconverted "half a SNP" is never silently applied: the CSV carries the
threshold that was actually used, and the reason is printed.

Robinson-Foulds needs **four or more taxa**. Points where fewer than four truth
leaves are present (very early in a run, or the 3-genome outbreak trio) get an
empty RF and a note, never a zero. `leaves_shared` says how much of the tree
existed at that instant.

Outputs: `accuracy_over_time__<dataset>__<config>.csv` per run,
`accuracy_over_time.csv` combined, and `accuracy_over_time*.png` when matplotlib
is available. CSV columns:

```
dataset, config, t_ms, t_s, seq, event, cluster, focus, rf_source,
leaves_run, leaves_truth, leaves_shared, splits_run, splits_truth, splits_shared,
rf, rf_normalized, collapse_threshold, note
```

### `baselines.py` — the outside baselines (Linux only, not part of the program)

Mash + neighbour-joining is the fast blurry reference; SKA2 is the exact
alignment-free reference. **Neither is installed on this machine**, so every mode
detects the absence and degrades to a pre-computed result instead of failing.
A pre-computed SKA2 table for the real 32-genome cluster already exists and is
copied into `bench-work/ncbi_pd_cluster/truth/` by `real_datasets.py`.

```bash
scripts/baselines.py status                       # what is installed, what is cached
scripts/baselines.py ska2   --dataset ncbi_pd_cluster --tree
scripts/baselines.py mash   --dataset ncbi_pd_cluster
scripts/baselines.py compare --dataset ncbi_pd_cluster \
    --run-dir bench-work/matrix/ncbi_pd_cluster/baseline
scripts/baselines.py all --dataset ncbi_pd_cluster --run-dir <run dir>

# force pre-computed input, never invoke a tool
scripts/baselines.py ska2 --dataset X --ska-dist /path/to/cluster_reads_distance.tsv --no-run
```

- `ska2` parses a `ska distance` TSV, writes a normalised `ska2.dist.tsv`,
  reports row/sample counts and min/median/mean/max distance, and can build a
  neighbour-joining tree from the distances (NJ is implemented in `common.py`, so
  no external phylogenetics package is needed).
- `mash` runs `mash sketch` + `mash dist -t` when mash exists, otherwise reads a
  cached/pre-computed table (headerless `mash dist -t` layout or a headered one),
  then builds the NJ tree. With neither, it records `status: "unavailable"` and
  exits 0 — a missing outside tool is a reported gap, not an error.
- `compare` is the one that matters: it joins SKA2's pairwise distances against
  Infinity-Win's `cluster_<id>.dist.tsv` by genome-name pair and reports pairs in both /
  only in one, exact matches, within-1-SNP and within-2-SNP counts, mean, median,
  p95 and max absolute difference, and Pearson/Spearman correlation. It also
  computes the Robinson-Foulds distance between `tree.nwk` and the dataset's
  reference tree. Pairs that exist in only one table are written to their own TSV
  rather than dropped.

The neighbour-joining trees are the *baselines'* trees. They are built here
because section 2 forbids wrapping external tools inside the engine, and section
11 keeps the baselines outside the program.

### `make_report.py` — `REPORT.md`

```bash
scripts/make_report.py                          # bench-work/matrix/summary.json -> REPORT.md
scripts/make_report.py --stdout                 # also print it
scripts/make_report.py --no-figures             # skip regenerating the curves
scripts/make_report.py --out /tmp/REPORT.md     # write somewhere else
```

Sections: verdict summary → datasets → per-dataset results (time to first tree,
time to exact tree, peak memory, counters, and a separate accuracy table) →
verdict detail per switch and dataset → acceptance check against published values
→ **switch invariance** → the key figure → outside baselines → **failures, gaps
and data-quality notes** → performance targets (section 10.4) → reproduction
commands.

The verdict is exactly section 11's rule: *a part of the idea pays if turning it
off makes time-to-exact or peak memory at least 20% worse on at least two of the
three datasets.* Four outcomes are distinguished, because collapsing them would
hide missing data:

| Verdict | Meaning |
|---|---|
| **pays** | ≥ 2 datasets measurably ≥ 20% worse |
| does not pay | ≥ 2 datasets measurable, fewer than 2 got 20% worse |
| undetermined | only 1 dataset measurable — the rule cannot decide either way |
| not measurable | no dataset had a usable time-to-exact or peak memory |

"Undetermined" and "not measurable" are never reported as a pass.

**Switch invariance** (section 6 of the report) checks PROJECT.md section 10.2
property 6 — that flat mode, hierarchical mode and every cache switch give
*identical* variants, distances and trees. Each configuration's pairwise SNP
counts are compared with the baseline's; any difference is reported as a defect
with examples and the size of the disagreement, because a switch that is faster
but wrong has not been shown to pay. `single_cluster` is excluded from that
count and annotated instead: it changes which genome every member is compared to,
so a different comparable-pair set is expected by design, not a violation.

The report also carries an **acceptance check** for the Salmonella trio: the
published pairwise SNP distances in `scripts/expected/salmonella_trio.json`
(isolate-A–isolate-B = 2, isolate-A–isolate-C = 2, isolate-B–isolate-C = 0, each with
its provenance) against what Infinity-Win reported in `cluster_<id>.dist.tsv`. A pair
Infinity-Win never compared is "not measured", not a pass.

### `release.sh` — release binaries

```bash
scripts/release.sh --check      # report what would happen, change nothing
scripts/release.sh              # cargo build --release, stage into ./dist
scripts/release.sh --no-build   # stage binaries that are already built
scripts/release.sh --windows    # attempt a cross build; refuses with a clear
                                # message when mingw is absent (it is, here)
```

Copies `infinity` and `Infinity-Win` from `$CARGO_TARGET_DIR` (`/tmp/infinity-target` on
this machine, via `.cargo/config.toml`) into `./dist/<triple>/`, and writes
`dist/VERSION` and `dist/SHA256SUMS.txt`. Windows x86-64 binaries come from
`.github/workflows/release.yml`, which builds both platforms, uploads them as
artefacts, runs these script tests, and publishes a GitHub release on a `v*` tag.

### `common.py` — shared library

Not a script. Defensive readers (`run.json`, `events.jsonl`, every TSV), Newick
parsing, branch collapsing, Robinson-Foulds, neighbour-joining, the SKA2 parser,
statistics, subprocess running with wall-time and memory measurement, the
configuration matrix, the verdict rule and the acceptance-check logic. Every
reader returns `None` plus a recorded *problem* string instead of raising on a
field the still-being-built CLI does not emit yet.

### `expected/salmonella_trio.json` — published ground truth

The three pairwise SNP distances for the outbreak trio, each with the analysis
that produced it, the file and field it was read from, the raw and filtered
counts, and an independent cross-check. Tolerance is ±1 SNP with the reason
written down: Infinity-Win compares **assemblies**, while these numbers come from
**reads** (SKA2) and from a mapping-plus-filtering pipeline, so a one-SNP
difference is evidence about the assemblies rather than automatically a bug.
`make_report.py` prints expected, observed and the difference for every pair, so
a near-miss is visible.

---

## 4. The `infinity bench` interface these scripts use

This matches the **implemented** CLI (`infinity bench --help`), not an assumed one:

```
infinity bench --output <run dir>          # an Infinity-Win output folder (run.json + trees)
            --truth  <truth dir>        # a `infinity simulate` truth folder
            [-j, --json <path>]         # default <output>/bench.json
            [--collapse-below <float>]  # collapse branches shorter than this (default 0.5)
            [--quiet]                   # write the JSON, do not print the table
```

`--truth` is **required** and must be a simulator truth folder: `tree_*.nwk`
(one truth tree per cluster, branch lengths in SNPs), `pairwise_snps.tsv`
(`genome_a  genome_b  snps`) and `clusters.tsv` (`genome  cluster_id`).

**Consequence for the real datasets.** They have no simulator truth folder — they
have a SKA2 distance table and an IQ-TREE reference tree instead. So
`run_matrix.py` *skips* `bench` for them, records why, and the comparison is made
by the scripts instead:

| dataset | accuracy comes from |
|---|---|
| simulated (`truth.dir` is a simulator folder) | `infinity bench` → `bench.json` |
| `ncbi_pd_cluster` | `baselines.py compare` (SKA2 join + RF vs `cluster.treefile`) and `accuracy_over_time.py --rf-mode internal` |
| `salmonella_trio` | the acceptance check in `make_report.py` (published SNP distances) and `baselines.py compare` |
| `s_aureus_cc` | nothing — no truth exists; counters and timings only |

`common.is_simulator_truth_dir()` is what decides this, so a reference folder is
never handed to `--truth` where it would fail and look like an Infinity-Win bug.

### `bench.json` as written by the CLI

```jsonc
{
  "output": "...", "truth": "...",
  "run": { /* the whole run.json, so bench.json is self-contained */ },

  "distance": {
    "pairs_in_truth": 20, "pairs_compared": 20, "pairs_missing_from_output": 0,
    "pairs_exact": 16, "max_abs_error": 1, "sum_abs_error": 4,
    "sum_true_snps": 210, "sum_reported_snps": 212,
    "worst_pair": ["sim_c00_s000", "sim_c00_s004", 10, 11]   // a, b, true, reported
  },

  "clusters": {
    "genomes_in_truth": 10, "genomes_compared": 10, "true_clusters": 2,
    "reported_clusters": 2, "recovered_exactly": 2, "pairwise_f1": 1.0
  },

  // Final RF, one entry per cluster tree. There is NO single full-tree RF.
  "trees": [
    { "source": "cluster_1", "truth": "tree_c00", "leaves": 5, "rf": 1,
      "splits_emitted": 2, "splits_true": 1, "rf_normalised": 0.3333,
      "t_ms": 75.1, "approximate": false, "focus": 2, "units": "snps" }
  ],

  // One entry per emitted tree, i.e. the accuracy-over-time series.
  // NOTE: for live trees the engine currently reports rf_normalised = null
  // (it does not score them against truth), so this series has no usable time
  // resolution yet - see "Where the curve comes from" below.
  "timeline": [ { "source": "backbone (live)", "leaves": 2, "rf": 0,
                  "splits_emitted": 0, "splits_true": 0, "rf_normalised": null,
                  "t_ms": 24.9, "approximate": true, "focus": 0,
                  "units": "subs_per_site" } ],

  "problems": []
}
```

Note `rf_normalised` (British spelling) — `common.BenchJson` reads both
spellings. Because RF is per cluster, the aggregate the report shows is the
**sum** of `trees[].rf` and the **mean** of `trees[].rf_normalised`, and both are
labelled as such rather than presented as one tree's distance.

Exit status matters: `bench` exits non-zero when it cannot score the run (no
truth folder, no `tree.nwk`, unparsable tree). `run_matrix.py` treats "non-zero
exit and no `bench.json`" as `bench.status = "failed"` with the stderr tail
recorded, and keeps the run's own timings and counters, which are still valid.

### Where the accuracy-over-time curve comes from

`accuracy_over_time.py --rf-mode auto` picks the best available source, and says
which one it used in the CSV's `rf_source` column and in its notes:

1. a `bench.json` series **with RF at two or more distinct instants** (the
   engine's own numbers — preferred whenever it can actually form a curve);
2. otherwise RF **recomputed in Python** at every emitted tree: each emitted
   cluster tree against the matching `tree_*.nwk` truth tree (simulated data), or
   the assembled full tree against one reference tree (real data);
3. otherwise the `bench.json` series anyway, labelled, because one point is
   better than none.

Today the engine's `timeline` reports `rf_normalised: null` for live trees, so on
simulated data path 2 is what produces the curve. If `bench` later scores live
trees against truth, path 1 takes over automatically with no script change. Use
`--rf-mode bench` to force the bench values, `internal` to force the recomputed
ones, or `infinity` to have the CLI re-run `bench` first.

### Output files the scripts parse, and how leniently

| File | Parsed by | Leniency |
|---|---|---|
| `run.json` | `run_matrix`, `make_report` | every key optional; null timings stay `null`; unknown keys ignored |
| `events.jsonl` | `accuracy_over_time` | bad lines counted and skipped with their line numbers; envelope shape `{"seq","t_ms","wall_ms","event":{"event","data"}}` |
| `cluster_<id>.dist.tsv` | `baselines`, `make_report`, `common` | the real header is `genome_a name_a genome_b name_b snps indels comparable_bp subs_per_site`; **names are matched before ids**, because joining on the numeric id would match nothing and look like a clean zero. Alternatives (`sample1/sample2`, `snp_count`, `sites`, `distance`) are also accepted; a file whose two genome columns cannot be identified is reported and skipped, never guessed at |
| `tree.nwk`, `cluster_<id>.nwk` | `accuracy_over_time`, `baselines` | internal node labels, `[&...]` comments, quoted names, multiple trees |
| `genomes.tsv` | `run_matrix`, `baselines` | `id name path length contigs n50 gc cluster is_pivot focus snps_vs_pivot indels_vs_pivot uncovered_bp novel_bp` |
| `clusters.tsv`, `cluster_<id>.variants.tsv` | available via `common.load_*` | same alias matching; ragged rows blank-filled and reported |

A column may be missing while the tool is still being built. Readers report that
as a problem string that ends up in `REPORT.md` section 8 — they do not crash and
do not substitute a default.

---

## 5. Tests

```bash
python3 -m unittest discover -s scripts/tests -t scripts -v
python3 -m compileall -q scripts          # every script must compile
```

`scripts/tests/fixtures.py` writes an executable **fake `infinity`** that speaks the
documented CLI contract and emits every output file, so the whole pipeline is
tested end to end without the real binary. Environment variables steer it:

| Variable | Effect |
|---|---|
| `FAKE_INFINITY_FAIL=1` | exit non-zero after writing a partial folder |
| `FAKE_INFINITY_NO_RUN_JSON=1` | exit 0 without writing `run.json` |
| `FAKE_INFINITY_NO_TIMINGS=1` | null every timing |
| `FAKE_INFINITY_NO_PEAKMEM=1` | `peak_memory_bytes: null` |
| `FAKE_INFINITY_IGNORE_SWITCH=1` | accept the flags but echo default params |
| `FAKE_INFINITY_BENCH_FAIL=1` | `bench` exits non-zero, writes nothing |
| `FAKE_INFINITY_DIST_JSON=<path>` | take pairwise SNP counts from a JSON list |
| `FAKE_INFINITY_SLEEP=<float>` | make wall time visible |

Covered: the SKA2 parser (including the real 496-row table when present), the
`events.jsonl` parser, the `run.json` reader, defensive table reading, Newick /
branch collapsing / Robinson-Foulds / neighbour-joining, the `--rf-from` path,
the verdict rule (pays, does not pay, undetermined, not measurable, exactly 20%,
zero baselines), the acceptance check, matrix resumability, `--dry-run`, every
failure path, the Mash/SKA2 degradation, the report contents, and
`real_datasets.py` copying, verification and sampling. All fixtures are created
in `tempfile.mkdtemp()` directories under `/tmp` and removed afterwards.

---

## 6. Assumptions worth knowing

1. **Sequential runs.** Configurations run one at a time. Parallel runs would
   make wall time and peak memory meaningless.
2. **Genome names must match** between Infinity-Win's output, the SKA2 table and any
   truth tree, or the joins produce nothing. Infinity-Win uses the FASTA file stem; the
   real datasets were built so those agree. `baselines.py compare` reports a
   zero-overlap join explicitly instead of printing an empty comparison.
3. **Three taxa have no RF.** The Salmonella trio contributes distance accuracy
   and the acceptance check only.
4. **The IQ-TREE tree is a reference, not ground truth** — a maximum-likelihood
   tree from a SKA2 split-k-mer core alignment. Infinity-Win builds
   neighbour-joining trees, so exact agreement is not expected; the RF distance
   is reported as a number, not judged pass/fail. It also has **43 leaves while
   the dataset has 32 assemblies**: 11 reference/outgroup genomes (`GCF_*`) are
   in the tree but were not part of the cluster download. Comparisons are
   restricted to the shared leaves, and `real_datasets.py` records the exact
   mismatch in `dataset.json` under `truth.tree_leaf_check`.
5. **SKA2 distances came from reads, Infinity-Win's from assemblies.** Small
   disagreements are expected and are reported per pair.
6. **`s_aureus_cc` has no accuracy truth.** It tests clustering separation and
   shows where the block shortcut stops paying, nothing more.
7. **Simulated datasets** (section 11 items 1 and 2) are supported by the same
   code path: give `bench-work/<name>/dataset.json` a `truth.dir` (and optionally
   a `simulate.command` list, which `run_matrix.py` runs once to materialise the
   assemblies and truth) and the matrix, curves and report pick it up unchanged.
