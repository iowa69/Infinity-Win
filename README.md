# M-Tato

**Progressive, cluster-first, position-aware comparison of nearly identical
genomes, with a live tree that sharpens while the computation runs.**

M-Tato takes a set of assembled genomes and answers, in that order and as fast as
it can: *which genomes are near each other* → *where do they differ* → *what
exactly differs, base by base*. Every number it shows carries a **focus level**,
so you always know whether you are looking at a sketch estimate or an exact count.

```
Focus 0   sketch only          seconds    "these 41 genomes form a cluster"
Focus 1   blocks compared      ~10x more  "this member differs in 8 places"
Focus 2   exact, base level     minutes    "2 SNPs, at these positions, these bases"
```

It is designed for bacterial assemblies of 0.5–15 Mb, from tens up to a few
thousand genomes, that are **already close** — outbreak investigation,
transmission clustering, within-host evolution. It is not a phylogenetics package
and it is not a read mapper.

---

## Install

### From source (Linux and Windows)

You need a stable Rust toolchain. There is no C toolchain requirement, no
Python runtime, and no external bioinformatics tool to install: the aligner, the
neighbour-joining tree code, the Newick reader/writer and the Robinson–Foulds
distance are all implemented in the project.

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone <this repository> && cd M-Tato
cargo build --release
```

Two binaries land in `target/release/`:

| Binary | What it is |
|---|---|
| `mtato` | the command-line tool: `run`, `simulate`, `bench` |
| `mtato-gui` | the desktop program |

On Windows the same commands work in PowerShell or `cmd`. Prebuilt binaries for
Windows x86-64 and Linux x86-64 are produced by `.github/workflows/release.yml`
and attached to the release; see `DECISIONS.md` D4 for why cross-compiling from
Linux is not part of the local build.

---

## Quick start

Point it at a folder of assemblies. One file is one genome; `.fa`, `.fasta` and
`.fna`, plain or gzipped, are accepted.

```sh
mtato run ./assemblies -o ./results
```

```
[    0.00s] discover  312/312
[    2.14s] scan      312/312
[    2.31s] clustered 312 genomes into 12 clusters at d<=0.005 (largest 88)
[    2.34s] backbone tree ready over 12 clusters  <-- first blurry tree
[   41.87s] refine    12/12
genomes 312 | clusters 12 | focus 0/1/2 = 0/0/312 | first tree 2.34s | exact 41.90s
output written to ./results
```

Then look at the results:

```sh
column -t -s$'\t' results/clusters.tsv      # which genomes cluster together
less results/cluster_3.variants.tsv         # every difference, with positions
cat results/tree.nwk                        # the full tree
```

Or open the folder in the GUI, which reads a finished run without recomputing it:

```sh
mtato-gui
```

### A worked example with known truth

The simulator produces genomes together with the answer, which is the fastest way
to see whether the tool is behaving:

```sh
mtato simulate -o /tmp/sim --seed 42 --length 2000000 -m 1 -n 50 --mean-snps 3
mtato run /tmp/sim/genomes -o /tmp/out
mtato bench --output /tmp/out --truth /tmp/sim/truth
```

```
  distances:
    pairs in truth 1225 | compared 1225 | missing 0
    exact 1225 | max |error| 0 | mean |error| 0.000
    fraction exact 1.0000 | sum reported / sum true 1.0000
  clusters:
    true 1 | reported 1 | recovered exactly 1 | pairwise F1 1.0000
  trees:
    cluster_0  vs tree_c00  leaves 50  RF 0  normalised 0.0000
```

---

## `mtato run`

```
mtato run <files or folders...> -o <outdir> [options]
```

### Parameters (defaults from section 6.1 of the design outline)

| Flag | Default | Meaning |
|---|---|---|
| `--k` | 31 | Anchor k-mer length. Must be odd, at most 31. |
| `--densities` | `16384,2048,256` | One anchor per this many k-mers at levels 1, 2, 3. Each must divide the previous. |
| `--sketch-level` | 2 | Level whose anchors form the sketch used for clustering. |
| `--cluster-dist` | 0.005 | Sketch distance threshold for joining a cluster (≈99.5% identity). |
| `--max-span-ratio` | 10 | Collinearity limit on span ratio. |
| `--max-span-extra` | 50000 | Collinearity limit slack, in bases. |
| `--dp-max-len` | 5000 | Longest trimmed piece aligned by dynamic programming. |
| `--min-identity` | 0.70 | Below this, a piece is reported as one `complex` variant. |
| `--xdrop` | 15 | Stop rule for ungapped extension at contig ends. |
| `--threads` | all cores | Worker threads. `0` means all. |
| `--memory-budget` | ¼ of RAM | Bytes of non-pivot sequence kept before re-reading from disk. |

If the median genome is under 500 kb the densities are automatically switched to
`1024,128,16` and a warning is logged — at 256-base spacing a small genome would
produce too few anchors to map.

### Experiment switches

Each of these removes one part of the idea so that its cost can be measured. They
exist for the benchmark in section 11 and are not normally used.

| Flag | Removes |
|---|---|
| `--flat` | The nested levels; only leaf blocks are used. |
| `--no-block-skip` | Block hashing entirely; every stretch between shared anchors is compared directly. |
| `--no-piece-cache` | Resolving a shared difference once. |
| `--no-profile-grouping` | Handling identical genomes once. |
| `--single-cluster` | Cluster-first; every genome is compared to one pivot. |
| `--focus0-only` | Refinement; stop at the blurry tree. |

`--flat`, `--no-piece-cache` and `--no-profile-grouping` leave the variants,
distances and trees unchanged, which is asserted by a test. `--no-block-skip`
does not: see `DECISIONS.md` D16.

### Output folder

| File | Content |
|---|---|
| `genomes.tsv` | One row per genome: length, contigs, N50, GC, cluster, is_pivot, focus, SNPs vs pivot, indels vs pivot, uncovered bp, novel bp. |
| `clusters.tsv` | Cluster id, pivot, pivot name, size, focus. |
| `tree.nwk` | The full tree. Branch lengths inside clusters are substitutions per site, so they are on the same scale as the backbone. |
| `cluster_<id>.nwk` | One cluster's tree, branch lengths in SNP counts. |
| `cluster_<id>.dist.tsv` | Pairwise SNP counts, indel counts and comparable bp. |
| `cluster_<id>.variants.tsv` | Every difference from the pivot: genome, pivot contig, pivot position (1-based), ref, alt, type, genome contig, genome position, strand, flags. |
| `run.json` | Parameters, version, per-stage timings, counters, peak memory. |
| `events.jsonl` | The timestamped event stream, one JSON object per line. |

**Coordinates.** Inside the program everything is 0-based and half-open. In the
files it is 1-based, converted in one place. For an insertion the reported
position is that of the base the insertion sits *after*, so an insertion and the
matching deletion report the same coordinate. That convention applies to **both**
coordinate columns of a variant row: `pivot_pos` and `genome_pos` are converted
the same way, so for a forward-strand insertion the two are the same number and
the row is self-consistent.

**Variant types** are `snp`, `ins`, `del` and `complex`. A `complex` variant means
the piece was too divergent or too long to resolve into individual events; it is
excluded from coverage rather than being guessed at.

**Flags** on a variant: `repeat` (the block matched a pivot repeat, so the copy
chosen may be wrong), `low_identity`, `align_failed`, `end_ext` (found by
end-of-contig extension rather than between two anchors), `dup` (the interval is
covered by two or more places of the member).

**`run.json` timings can be `null`.** A `null` means the stage did not complete —
for example `exact_tree_ms` on a cancelled run. It is never written as `0`,
because "instant" and "never" are different claims.

---

## `mtato simulate`

Deterministic for a given seed: the same seed and configuration produce
byte-identical FASTA files and truth files.

```sh
mtato simulate -o <outdir> [--seed 1] [--reference REF.fasta | --length 2000000 --gc 0.5]
               [-m CLUSTERS] [-n LEAVES] [--divergence 0.02]
               [--mean-snps 3] [--mean-indels 0] [--mean-segments 0] [-c CONTIGS]
```

It writes one FASTA per leaf genome under `genomes/`, the cluster roots under
`roots/`, and the truth under `truth/`:

| File | Content |
|---|---|
| `truth/tree_c<NN>.nwk` | The true tree per cluster, branch lengths in SNPs. |
| `truth/pairwise_snps.tsv` | True pairwise SNP counts. |
| `truth/clusters.tsv` | True cluster of each genome. |
| `truth/mutations_<genome>.tsv` | Every mutation on each lineage. |
| `truth/manifest.json` | The configuration and seed, so a dataset can be replayed exactly. |

Every SNP is placed at a position never mutated before on that lineage, so the
true distance between two leaves is exactly the path length between them in the
truth tree. That is what makes the truth tests meaningful.

---

## `mtato bench`

Compares an output folder with a truth folder and writes the measurements as
JSON. This is what the experiment in section 11 is built on.

```sh
mtato bench --output <outdir> --truth <truthdir> [--json bench.json] [--collapse-below 0.5] [--quiet]
```

It reports: distance accuracy (pairs exact, maximum absolute error, reported/true
total), cluster recovery (exact clusters and pairwise co-clustering F1), and tree
accuracy as a Robinson–Foulds distance after collapsing branches shorter than
`--collapse-below` SNPs — the resolution section 11 asks for, so that a branch of
less than half a SNP cannot count as a topological difference. It also reads
`events.jsonl` to give every emitted tree a timestamp, which is the x axis of the
accuracy-over-time figure.

Trees are matched by leaf set, not by cluster number, because the reported
cluster ids need not agree with the truth's. A `tree.nwk` comparison is reported
in substitutions-per-site units and a `cluster_<id>.nwk` comparison in SNP units;
the JSON records which, because the collapse threshold is only meaningful in SNP
units.

---

## The desktop program

One window, four panels: genomes on the left, the tree in the middle, parameters
and live counters on the right, and the map along the bottom.

- **Live.** The window never blocks while computing. Every stage streams events,
  and the tree redraws as focus improves: Focus 0 branches are dashed and pale,
  Focus 1 solid and thin with an `≈` marker, Focus 2 solid and dark. Nodes move
  smoothly to their new places when a cluster is refined.
- **Clusters start collapsed** as a triangle with the member count. Click to open.
  Clicking a cluster also moves it to the front of the refinement queue, so you
  can look at the cluster you care about first.
- **The cluster-distance slider takes effect immediately.** It re-clusters from
  the stored sketches; clusters whose pivot and members are unchanged keep their
  results, and only the others are refined again.
- **The map panel** shows the pivot as a horizontal axis with one row per selected
  member: grey line where identical, orange box for a difference piece not yet
  resolved, red tick for a SNP, blue tick for an indel, purple box for a complex
  variant, triangle for a breakpoint, gap where not covered. In hierarchical mode
  a region that failed at level 1 is shown as a wide pale box that narrows as
  levels 2 and 3 are checked — the zoom made visible.
- **Select two leaves** to see their pairwise distance, SNP count and focus.
- **Export** the tree as Newick or PNG, plus the distance matrices, variants and
  clusters.
- **Open a finished run** from its output folder and browse it without
  recomputing.

Genome folders are added by dragging them onto the window or by typing a path,
with a list of recently used folders. A native file-picker dialog is deliberately
not used; see `DECISIONS.md` D3.

---

## How it works, briefly

1. **One hash pass feeds every level.** Sliding a window of `k` bases along each
   contig produces canonical k-mer hashes. The same hashes serve as the sketch
   (blurry comparison), as anchors (mapping) and as cut points for nested blocks
   (zoom).
2. **Nested content-defined blocks.** Contigs are cut at anchors into blocks of
   about 16 kb, 2 kb and 256 bp. A block whose hash the pivot already knows is
   identical and is skipped without looking inside. Block hashes are
   strand-independent, so a reverse-complemented contig still matches.
3. **Cluster first.** Genomes are clustered greedily against pivots only, through
   an inverted index of sketch hashes — never all-versus-all. The most contiguous
   assembly of each cluster becomes its pivot. Fine comparison happens inside
   clusters.
4. **Position is never lost.** Shared points anchor the member onto the pivot;
   the stretches between them are compared base by base. Every difference maps
   back to a position on the pivot, so variants can be reported, flagged and
   drawn.

Full detail is in `PROJECT.md` sections 5–6.

---

## Development

```sh
cargo test --workspace                      # 152 core unit tests + 19 integration tests
cargo clippy --workspace --all-targets      # must be warning-free
cargo fmt --all -- --check
cargo check --target x86_64-pc-windows-msvc # portability gate
```

The integration tests in `tests/` are the ones that matter most:

- `tests/tests/properties.rs` — the six properties of section 10.2 that must
  always hold, including "same input gives byte-identical output with 1 thread and
  with many" and "every optimisation switch gives the same answer".
- `tests/tests/truth.rs` — the truth tests T1–T5 of section 10.3, run against
  simulated genomes with known answers.

### Benchmark

```sh
scripts/run_matrix.py --dry-run     # show the configuration x dataset matrix
scripts/run_matrix.py               # run it
scripts/make_report.py              # write REPORT.md
scripts/README.md                   # the full procedure
```

---

## Limitations

These are real and are described in `PROJECT.md` section 15. `DECISIONS.md`
records where the implementation deviates from the design, with measurements.

- **Pivot bias.** Sequence the pivot lacks is not compared between members. It is
  reported as novel bp, not as distance.
- **Fragmented assemblies cost coverage.** A member contig that runs across a
  junction between two pivot contigs cannot be placed. At 100 contigs over a 2 Mb
  genome this loses about 6% of the true SNP total; at 20 contigs about 2%; at one
  contig nothing. Measured and documented in `DECISIONS.md` D13.
- **Repeats.** A block that occurs more than once in the pivot can be mapped to
  the wrong copy. Variants there are flagged `repeat`, not corrected. On real
  data this matters: on three *Salmonella* Heidelberg assemblies the tool reports
  2 chromosomal SNPs, matching the published figure, plus 13 differences on a
  1.4 kb element at 139× coverage that the reference pipeline filters out. Both
  are in the variant file, and it is up to you to decide which you want.
- **Trees.** Neighbour joining on SNP counts ignores recombination and is not a
  substitute for a full phylogenetic analysis.
- **Divergent clusters.** The block shortcut stops helping around 0.5% divergence
  and anchors become too sparse above a few percent. Keep `--cluster-dist` small.

---

## Licence

MIT. See `LICENSE`.
