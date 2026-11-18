# Infinity-Win: project outline

**Progressive, cluster-first, position-aware comparison of nearly identical genomes, with a live tree that sharpens while the computation runs.**

| | |
|---|---|
| Status | Design outline, version 1, 2026-11-18. Ready for a builder. |
| Name | Infinity-Win (owner, 2026-11-18). Earlier working name: M-Tato. |
| Reader | An AI coding assistant (or a human developer) who builds this without access to the design conversation. |

---

## 0. How to use this document

- Sections 1-4 explain **why** the project exists and what was considered.
- Sections 5-9 say **what to build**.
- Sections 10-11 say **how the result is judged**.
- Sections 12-13 say **how to build it and in which order**.
- Section 14 lists **who decided what** and the questions still open.

Two labels are used throughout:

- **(owner)**: decided by the project owner. Do not change.
- **(author's call)**: chosen by the author of this outline because the owner left it open. The owner may override it. The builder may not, except through the rule below.

Rules for the builder:

1. Follow the build order in section 13. Each milestone has acceptance checks. Do not start the next milestone until they pass.
2. The engine must never depend on the GUI. Everything the GUI shows must also be reachable from the command line.
3. Same inputs and parameters must give byte-identical output files, whatever the thread count.
4. If something here turns out to be wrong, impossible, or ambiguous, do not silently change the algorithm. Pick the closest working alternative, and record what changed and why in `DECISIONS.md` at the repository root.
5. Do not add features that are not in this document.
6. The program needs no network access at run time.
7. Library names in section 12 are suggestions from memory. Check that each exists and is maintained before depending on it.
8. Never present an approximate number as exact. Every distance shown or written carries its focus level (section 5).

---

## 1. What this is and why

Tools like Mash compare genomes very fast by reducing each genome to a small sample of k-mer hashes. That speed has a price: the sample forgets *which* k-mers were kept and *where* they sit, and it is far too sparse to see single mutations. A single SNP changes 21 of roughly 5 million k-mers in a bacterial genome, so a 1,000-hash sketch notices it about 0.4% of the time.

Exact tools (alignment, split k-mers) see every SNP but give nothing until they finish, and many of them compare everything against everything.

Infinity-Win tests one idea **(owner)**:

> Treat genome comparison like a photo coming into focus. Start with a fast blurry answer, then raise the resolution step by step, and only where genomes actually differ. Keep the identity and position of the k-mers at every step so that each level can build on the previous one and every difference can be mapped back onto the genome.

The owner's requirements:

- **Target (owner):** measure divergence between nearly identical genomes, down to single SNPs, in minutes.
- **Coarse to fine (owner):** contigs, then shared content, then exact positions. The tree gets "higher definition" as the run proceeds.
- **Position-aware (owner):** like Mash and MeShClust in speed and spirit, but k-mers carry what they are and where they are. Mapping matters.
- **Cluster first (owner):** genomes are clustered early. Fine divergence is measured inside clusters, never all-versus-all.
- **Visual (owner):** a desktop program with a live interface that shows the tree and the parameters. A usable picture early is worth more than waiting for the exact one.
- **Platforms (owner):** Windows and Linux.
- **Purpose (owner):** find out whether this multidimensional approach actually pays. The project is an experiment as much as a tool, so measuring it honestly is part of the deliverable (section 11).

---

## 2. Goals, non-goals, anti-goals

### Goals

1. Take a set of assembled genomes (FASTA) and produce clusters, a tree, pairwise distances inside each cluster, and the list of differences with their positions.
2. Show results progressively: a blurry tree within seconds, exact SNP-level branches within minutes.
3. Do work in proportion to what differs. Identical stretches are recognised once and skipped.
4. Run as a desktop app on Windows and Linux, and as a command-line tool with the same engine.
5. Include a simulator and a benchmark that say, with numbers, whether each part of the idea pays.

### Non-goals (version 1)

- Raw reads (FASTQ) as input. Assemblies only.
- Large eukaryotic genomes. Design target is bacterial-size genomes (about 0.5-15 Mb), from tens up to a few thousand genomes **(author's call, see section 14)**.
- Maximum-likelihood phylogenetics, recombination detection, dating, annotation, gene calling.
- macOS, web, or mobile versions.
- Wrapping or bundling external tools (Mash, SKA2, and so on) inside the engine. They are used only as outside baselines in the benchmark.

### Anti-goals (these would make the project a failure even if it runs)

- A window that freezes while computing.
- A tree that reshuffles which genomes belong to which cluster after the user has seen them (without the user changing a parameter).
- Approximate numbers shown as if they were exact.
- Base-level comparison of all genomes against all genomes.
- Requiring the user to supply a reference genome.
- A benchmark that only reports the cases where the idea wins.

---

## 3. Prior art, and what is actually new

This table comes from the author's knowledge plus a short search. It is not a systematic review.

| Piece of the idea | Existing work |
|---|---|
| Fast k-mer sketches | Mash, sourmash, Dashing |
| Sequence clustering without all-vs-all | MeShClust, CD-HIT, Linclust |
| Several k sizes at once | PopPUNK (separates SNP divergence from gene content by comparing across k), CMash and MOS (multi-resolution sketches), variable-order de Bruijn graphs |
| SNPs straight from k-mers | SKA2 (split k-mers), kSNP, KmerAperture |
| Cluster, then refine inside clusters | PopPUNK followed by SKA, packaged as PopPIPE (separate tools chained) |
| Anchors plus local alignment | MUMmer, Parsnp, minimap2, skani |
| Storing repeated sequence once | prefix-free parsing, compressed pangenome indexes |

What Infinity-Win tries that the search did not find in one tool:

1. **One hash pass feeds every level.** The same sampled k-mers serve as the sketch (blurry), as position anchors (mapping), and as cut points for nested blocks (zoom).
2. **Nested content-defined blocks.** Genomes are cut at the anchors into blocks of about 16 kb, 2 kb, and 256 bp. A block whose hash is already known is identical and is skipped without looking inside.
3. **A staged, live result.** Every distance carries a focus level, and the tree is redrawn as focus improves.

Whether this beats existing tools is the hypothesis under test, not a claim.

---

## 4. Approaches considered

### 4.1 Comparison engine

| | Approach | Trade-off |
|---|---|---|
| **A (chosen)** | **Anchor pyramid with a pivot genome per cluster.** Sampled k-mers give sketch, anchors, and nested blocks. Each cluster member is compared to the cluster's pivot. Identical blocks are skipped by hash, the rest is compared base by base between shared anchors. | Exact for SNPs and indels, including SNPs close together. Work scales with differences. Weakness: sequence missing from the pivot is not compared between members (section 15). |
| B | **Nested split k-mers.** One SKA-style table per cluster, filled progressively by hash threshold. | Reference-free and statistically clean. But work at full resolution scales with genome size times genome count, SNPs closer than half a k-mer are missed, indels are invisible, and positions need extra storage. |
| C | **Chain existing tools behind a GUI** (Mash, a clusterer, SKA2 or Parsnp). | Fastest to build. But it cannot carry k-mer identity and position between stages, so it does not test the idea. Those tools are also poorly supported on Windows. |
| D | **Full all-k index** (suffix array, BWT, maximal exact matches). | Exact for every k at once. Heavy to build, nothing to show until the index exists. |

**(author's call)** A is chosen because it is the only one that keeps position at every level, gives an early answer, and handles indels. B stays relevant as the main outside baseline (SKA2).

### 4.2 Application shell

| | Approach | Trade-off |
|---|---|---|
| **A (chosen)** | **Rust engine + native Rust GUI (egui/eframe).** | One language, one build, one executable per platform, no runtime to install. Tree drawing is custom code. |
| B | Rust engine + web front end in a desktop wrapper (Tauri). | Richer ready-made visual libraries. Two toolchains and a system web-view dependency on Linux. |
| C | Python + Qt, calling native code. | Familiar to bioinformaticians. Packaging for Windows is fragile and the hot loops still need native code. |

**(author's call)** A is chosen to keep the build simple for an AI builder and the result easy to install.

---

## 5. Vocabulary

| Term | Meaning |
|---|---|
| **Genome** | One FASTA file. Each record is a contig. |
| **Anchor** | One occurrence of a sampled k-mer: its hash, contig, position, strand. |
| **Level** | Sampling density of anchors. Level 1 is the sparsest (coarse), level 3 the densest (fine). Every level-1 anchor is also a level-2 and level-3 anchor. |
| **Block** | The stretch of a contig from the start of one anchor to the end of the next anchor of the same level. Level-1 blocks are about 16 kb, level-2 about 2 kb, level-3 about 256 bp. Level-3 blocks are also called **leaf blocks**. |
| **Sketch** | The set of level-2 anchor hashes of a genome. Used for clustering. |
| **Cluster** | A group of genomes within a sketch distance threshold of the same pivot. |
| **Pivot** | The representative genome of a cluster. Differences inside the cluster are expressed in pivot coordinates. |
| **Shared point** | An anchor of a member genome whose position on the pivot is known. |
| **Piece** | The stretch between two consecutive shared points, in the member and in the pivot. |
| **Variant** | One difference between a member and its pivot: SNP, insertion, deletion, or complex. |
| **Profile** | Everything that distinguishes a member from its pivot: its variants plus the pivot intervals it does not cover. |
| **Focus** | How sharp a result is. **Focus 0** = sketch only. **Focus 1** = blocks compared, differing places known but not resolved. **Focus 2** = exact, base level. |

The "dimensions" of the structure are: **genome** (who), **position** (where), **level** (how fine), and, experimentally, **k** (how specific). A cell is "this block of this genome at this level", and its value is a hash. Only cells that differ from what is already known get expanded.

Coordinates are 0-based, half-open `[start, end)` inside the program. Output files use 1-based positions.

---

## 6. The algorithm

### 6.1 Parameters

| Name | Default | Meaning |
|---|---|---|
| `k` | 31 | Anchor k-mer length. Must be odd and at most 31 (fits in 64 bits). |
| `densities` | 16384, 2048, 256 | One anchor per this many k-mers at levels 1, 2, 3. Each must divide the previous. |
| `sketch_level` | 2 | Level whose anchor hashes form the sketch. |
| `cluster_dist` | 0.005 | Sketch distance threshold for joining a cluster (about 99.5% identity). |
| `max_span_ratio` | 10 | Collinearity limit, see 6.6. |
| `max_span_extra` | 50000 | Collinearity limit, see 6.6. |
| `dp_max_len` | 5000 | Longest piece (after trimming) that is aligned with dynamic programming. |
| `min_identity` | 0.70 | Below this, an aligned piece is reported as one complex variant. |
| `xdrop` | 15 | Stop rule for ungapped extension at contig ends. |
| `threads` | all cores | |

If the median genome length is under 500 kb, the default densities are too sparse. In that case use 1024, 128, 16 automatically and log it.

### 6.2 Reading and scanning a genome

1. Read the FASTA file (plain or gzip). Convert to upper case. Replace every character that is not A, C, G, or T with N.
2. Record per genome: name, total length, number of contigs, N50, GC fraction.
3. Slide a window of `k` bases along each contig. Skip windows containing N.
4. For each window compute the anchor hash exactly like this, so that results are reproducible everywhere:

```
A=0, C=1, G=2, T=3
fwd    = the k-mer as a 2k-bit unsigned integer, first base in the highest bits
rc     = its reverse complement, encoded the same way
canon  = min(fwd, rc)
strand = '+' if fwd <= rc else '-'
h      = fmix64(canon XOR 0x9E3779B97F4A7C15)

fmix64(x):
    x ^= x >> 33;  x *= 0xff51afd7ed558ccd
    x ^= x >> 33;  x *= 0xc4ceb9fe1a85ec53
    x ^= x >> 33
    (all arithmetic wraps at 64 bits)
```

5. The window is an **anchor at level j** if `h < (2^64 - 1) / densities[j]` (integer division). Because the thresholds are nested, a level-1 anchor is automatically a level-2 and level-3 anchor.
6. Keep every level-3 anchor as `(h, contig, position, strand)`.

Using the canonical k-mer makes the anchors identical whether a contig is read forward or reverse-complemented. This is what "bidirectional" means in this project.

### 6.3 Blocks and block hashes

Within one contig of length `L`, let the level-j anchors start at `p1 < p2 < ... < pn`.

- Interior blocks: `[p_i, p_(i+1) + k)` for `i = 1 .. n-1`. Each block begins with one anchor and ends with the next, so neighbouring blocks overlap by exactly `k` bases.
- Head block: `[0, p1 + k)`, omitted if `p1 = 0`.
- Tail block: `[pn, L)`, omitted if `pn + k = L`.
- No anchors at all: the whole contig is one block.

The hash of a block must not depend on strand:

```
hf      = xxh3_64(bytes of the block)
hr      = xxh3_64(bytes of the reverse complement of the block)
hash    = min(hf, hr)
flipped = (hr < hf)
```

Two blocks with the same hash are treated as identical sequence (possibly on opposite strands).

### 6.4 Focus 0: sketch, clusters, blurry tree

**Sketch.** The set of hashes of a genome's anchors at `sketch_level`.

**Distance.** For sketches A and B with `s` shared hashes:

```
J = s / (|A| + |B| - s)
D = 1                              if J = 0
D = -(1/k) * ln( 2J / (1 + J) )    otherwise
```

**Clustering without all-vs-all.** Greedy, against pivots only:

```
order genomes by (N50 descending, total length descending, name ascending)
pivots = empty;  index = map from hash to list of pivots containing it
for each genome g in that order:
    count shared hashes between g and every pivot, using the index
    best = pivot with the smallest D(g, pivot)      (ties: the earliest pivot)
    if best exists and D(g, best) <= cluster_dist:
        g joins best's cluster
    else:
        g becomes a new pivot; add its hashes to the index
```

The best assembly of each cluster therefore becomes its pivot, and the result is deterministic.

**Backbone tree.** Compute D between pivots only (there are few), build a neighbour-joining tree, root it at the midpoint. Its leaves are clusters.

**Blurry cluster view.** Until a cluster is refined, its members hang from the pivot as a star, each at its sketch distance. This is the Focus 0 picture.

### 6.5 Pivot index

For the pivot of a cluster being refined, hold in memory:

- its sequence;
- its level-3 anchors, sorted by position;
- **anchor map**: canonical k-mer hash to `(contig, position, strand)`, only for anchors whose k-mer occurs exactly once among the pivot's anchors;
- **block dictionaries**, one per level: block hash to list of `(contig, start, end, flipped)`.

Clusters with a single genome need no index.

### 6.6 Comparing a member to its pivot

For each contig of member `G`:

**Step A, block shortcut.** Look up each leaf block of `G` in the pivot's level-3 dictionary. A hit means the block is identical to a pivot block; mark it *matched* and remember the mapping. Everything else is *variable*. If the hash occurs several times in the pivot (a repeat), prefer the occurrence that continues directly from the previous matched block; otherwise take the first.

For a matched pair, member block `[gs, ge)` with `flipped_G` and pivot block `[ps, pe)` with `flipped_P`, let `rel = flipped_G XOR flipped_P`. A member position `x` maps to the pivot as:

```
rel = false:  P(x) = ps + (x - gs)
rel = true :  P(x) = pe - 1 - (x - gs)
```

**Step B, shared points.** Give every anchor of the contig a pivot position if one can be found:

- an anchor at the start or end of a matched block gets its pivot position from that block's mapping;
- any other anchor is looked up in the pivot's anchor map (unique anchors only), with `rel = strand_G XOR strand_P`.

A shared point is `(g, p_contig, p, rel)` where `g` and `p` are the start of the k-base anchor on the member and on the pivot's forward strand.

Drop a shared point that came from the anchor map if it is not collinear (defined below) with its nearest shared neighbours on both sides while those two neighbours are collinear with each other. This removes stray matches.

**Step C, pieces.** Walk the remaining shared points in member order. For consecutive points `x` and `y` (`gx < gy`):

- If the only thing between them is one matched block, there is nothing to do.
- Otherwise check **collinearity**: same pivot contig, same `rel`, pivot order consistent with `rel` (`py >= px` if `rel` is false, `py <= px` if true), and neither span longer than `max_span_ratio` times the other plus `max_span_extra`.
- Collinear: this is a **piece**.

```
member piece:  [gx, gy + k)
pivot piece :  [px, py + k)   if rel = false
               [py, px + k)   if rel = true
ref = pivot piece
alt = member piece, reverse-complemented if rel = true
```

  `ref` and `alt` begin with the same k bases and end with the same k bases. If `ref` equals `alt`, the piece is identical. Otherwise it is a **difference piece**.
- Not collinear: record a **breakpoint**; the member bases between the two points are *unplaced*.

**Step D, contig ends.** Between a contig end and its nearest shared point there is no second anchor. Extend outward from the shared point without gaps, comparing member and pivot base by base (score +1 for a match, -3 for a mismatch, stop when the score falls `xdrop` below its maximum, keep the extension up to the maximum). Mismatches inside the kept part are SNPs. Member bases beyond it are unplaced. A contig with no shared point at all is entirely unplaced.

**Focus 1 result** for the member: the list of difference pieces as `(p_contig, p_start, p_end, hash of alt)`, the breakpoints, and the unplaced bases. No base-level work has been done yet.

### 6.7 Focus 2: resolving a difference piece

```
diff(ref, alt):
    trim the longest common prefix, then the longest common suffix of what remains
    r, a = the remainders
    if both empty            -> nothing
    if r empty               -> insertion of a
    if a empty               -> deletion of r
    if len(r) = len(a) and (mismatches <= 3 or mismatches / len(r) <= 0.10)
                             -> one SNP per mismatching position
    else if max(len(r), len(a)) <= dp_max_len
                             -> global alignment with affine gaps
                                (match +1, mismatch -3, gap open -5, gap extend -1);
                                read SNPs, insertions, deletions from the alignment;
                                if identity < min_identity: one COMPLEX variant for r
    else                     -> one COMPLEX variant for r
```

Additional rules:

- A position where either base is N is never a SNP. It counts as not covered.
- Shift every insertion and deletion as far left as it will go (standard left-alignment), so the same event always gets the same coordinates.
- Each variant stores: pivot contig, pivot position, ref, alt, type, member contig, member position, strand, flags.

### 6.8 Coverage

For each member, record the pivot intervals it does **not** cover. A pivot base is covered if it lies in a matched block or in a compared piece, and is not inside a COMPLEX variant, an N, or an interval covered by two or more places of the member (flag `dup`). Member bases that are unplaced or inserted are counted as **novel bp**.

### 6.9 Profiles, distances, trees

**Profile.** A member's sorted variant list plus its uncovered pivot intervals. The pivot's own profile is empty.

**Identical profiles are handled once.** Group members with identical profiles, compute distances between groups only, and attach the members of a group to the tree as zero-length siblings.

**Pairwise distance inside a cluster**, for members A and B:

```
comparable(A,B) = pivot length - | uncovered(A) UNION uncovered(B) |
snps(A,B)       = number of comparable pivot positions where the base of A differs
                  from the base of B (a member without a SNP there has the pivot base)
indels(A,B)     = number of insertion/deletion variants, inside comparable intervals,
                  present in one profile and not identically in the other
dist(A,B)       = snps(A,B) / comparable(A,B)        (substitutions per site)
```

**Focus 1 estimate.** Before pieces are resolved, the distance between A and B is estimated as the number of difference-piece records that are in one list and not in the other. This is shown as approximate ("≈"). The true SNP count is at least half that number, rounded up.

**Cluster tree.** Neighbour-joining on `snps` (Focus 2) or on the Focus 1 estimate, negative branch lengths set to zero, rooted at the midpoint.

**Full tree.** Take the backbone tree and replace each pivot leaf with its cluster tree. Use `dist` (substitutions per site) for branch lengths inside clusters so that units match the backbone. Each cluster subtree is labelled with its focus level.

### 6.10 The zoom: nested levels and caches

The steps above work with leaf blocks only ("flat mode"). Build that first. Then add the following, each behind a switch, because each one is a part of the idea that the benchmark must be able to turn off.

**Top-down matching (hierarchical mode).** Start with level-1 blocks. A level-1 block found in the pivot's level-1 dictionary is identical: mark it matched and do not look inside. Otherwise look at the level-2 blocks inside it, and for those that miss, at the level-3 blocks. Only the missed leaf blocks are variable. Head and tail blocks are treated the same way. The final variants must be identical to flat mode.

**Piece cache.** Key: `(hash of alt, p_contig, p_start, p_end)`. Value: the variants of that piece. A mutation shared by 500 members of a clade is then resolved once.

**Profile grouping.** As in 6.9.

**Memory.** Sequences of non-pivot genomes are not kept after use unless they fit a memory budget (default: a quarter of physical RAM); otherwise the file is read again when needed. A member whose difference pieces are all in the piece cache does not need to be re-read for Focus 2.

### 6.11 Order of work

1. Scan all genomes in parallel. Emit each genome's statistics as it finishes.
2. Cluster, build the backbone tree, emit the Focus 0 picture.
3. Refine clusters one at a time, largest first, except that a cluster the user selects jumps to the front. For each cluster: build the pivot index, run 6.6 for all members in parallel and emit the Focus 1 tree, then run 6.7-6.9 and emit the Focus 2 tree.
4. Emit a final summary.

---

## 7. Engine interface

The engine is a library. It runs in background threads and communicates only through messages. The command-line tool and the GUI are two consumers of the same messages.

**Events (engine to consumer).** Field names are indicative.

| Event | Content |
|---|---|
| `GenomeScanned` | genome id, name, length, contigs, N50, GC, sketch size |
| `Clustered` | threshold, list of clusters with pivot and members |
| `BackboneTree` | Newick, leaves are clusters |
| `ClusterFocus` | cluster id, focus (0, 1, or 2), Newick, per-member numbers (difference pieces, SNPs, indels, uncovered bp, novel bp) |
| `MapUpdate` | cluster id, genome id, difference pieces and variants in pivot coordinates |
| `Counters` | blocks looked up, blocks matched per level, pieces compared, cache hits, alignments run, bases aligned |
| `Progress` | stage, done, total, elapsed time |
| `Warning` | text |
| `Finished` | summary and timings |

**Commands (consumer to engine).** `Start {inputs, parameters, output folder}`, `Prioritize {cluster}`, `SetClusterDist {value}`, `Cancel`.

Changing `cluster_dist` re-clusters from the stored sketches. Clusters whose pivot and members are unchanged keep their results; the others are refined again.

Every event is also appended, with a timestamp, to `events.jsonl` in the output folder. The benchmark uses this file to measure how accuracy grows over time.

---

## 8. The desktop program

One window. The layout below is the author's proposal **(author's call)**; the owner should adjust it to taste.

```
+--------------------------------------------------------------------------------+
| [Add genomes] [Run] [Stop]    cluster distance [--|-----] 0.005    threads [8] |
+----------------+-------------------------------------------+-------------------+
| GENOMES (312)  |  TREE                                     | PARAMETERS        |
| # c1 strainA01 |        +------ c3  (41 genomes)  blurry   | k            31   |
| # c1 strainA02 |   +----+                                  | levels 16k/2k/256 |
| # c2 strainB01 |   |    |   +-- strainA07        sharp     | cluster dist 0.005|
| # c2 strainB02 | --+    +---+                              +-------------------+
| ...            |   |        +-- strainA09        sharp     | LIVE COUNTERS     |
| filter: [    ] |   +----------- c2  (88 genomes)  medium   | blocks identical  |
|                |                                           |   99.7 %          |
|                |                    |-----| 5 SNPs         | pieces compared   |
|                |                                           | SNPs found        |
+----------------+-------------------------------------------+-------------------+
| MAP  (cluster c1, pivot strainA01)                                             |
| pivot     ================================================================     |
| strainA07 --------|--------------[==]------------------|------------------     |
| strainA09 --------|------------------------------------------   (gap)   --     |
+--------------------------------------------------------------------------------+
| Scan done 2.1 s | Clusters done 0.3 s | Blocks 7/12 clusters | Exact 3/12      |
+--------------------------------------------------------------------------------+
```

**Tree panel.**
- Rectangular tree, pan and zoom with the mouse, scale bar.
- Clusters start collapsed (a triangle with the member count) and open on click. Clicking a cluster also sends `Prioritize`.
- Focus is always visible: Focus 0 branches are dashed and pale, Focus 1 solid and thin with an "≈" marker, Focus 2 solid and dark. A legend explains this.
- When a cluster's tree changes, nodes move smoothly to their new place (about 300 ms).
- Hovering a leaf shows its numbers. Selecting two leaves shows their pairwise distance, SNP count, and focus.
- Search box to find a genome by name.

**Genome list.** Name, cluster colour, length, contigs, focus. Sortable, filterable. Selection is linked to the tree and the map.

**Parameters and counters.** Shows every parameter in 6.1 and the live counters. The cluster distance is a slider and takes effect immediately. Other parameters apply at the next run.

**Map panel.** The pivot of the selected cluster as a horizontal axis (contigs side by side). One row per selected member:
- identical stretches: thin grey line;
- difference piece not yet resolved: orange box spanning the piece;
- SNP: red tick; insertion or deletion: blue tick; complex: purple box; breakpoint: triangle; not covered: gap.

In hierarchical mode, show a region that failed at level 1 as a wide pale box that narrows as levels 2 and 3 are checked. This is the "zoom" made visible. Zoom and pan along the axis; hovering shows position and alleles.

**Status bar.** One segment per stage with progress and elapsed time.

**Export.** Tree (Newick), distance matrices, variants, clusters, and the tree picture (PNG; SVG if practical).

**Open results.** A finished output folder can be opened and viewed without recomputing.

**Responsiveness.** The window must never stall for more than about 50 ms. All computation happens off the interface thread.

---

## 9. Command line and output files

```
infinity run <files or folders...> -o <outdir> [--threads N] [--k 31]
          [--densities 16384,2048,256] [--cluster-dist 0.005]
          [--flat] [--no-block-skip] [--no-piece-cache] [--no-profile-grouping]
          [--single-cluster]
infinity simulate ...          (section 10)
infinity bench ...             (section 11)
Infinity-Win               (the desktop program)
```

Input files: `.fa`, `.fasta`, `.fna`, optionally `.gz`. One file is one genome. A folder means all such files inside it.

Output folder:

| File | Content |
|---|---|
| `genomes.tsv` | id, name, path, length, contigs, N50, GC, cluster, is_pivot, focus, SNPs vs pivot, indels vs pivot, uncovered bp, novel bp |
| `clusters.tsv` | cluster id, pivot, size, focus |
| `tree.nwk` | full tree |
| `cluster_<id>.nwk` | tree of one cluster |
| `cluster_<id>.dist.tsv` | pairwise SNP counts and comparable bp |
| `cluster_<id>.variants.tsv` | genome, pivot contig, pivot position (1-based), ref, alt, type, genome contig, genome position, strand, flags |
| `run.json` | parameters, software version, timings per stage, counters, peak memory |
| `events.jsonl` | timestamped event stream |

---

## 10. Testing

### 10.1 Simulator (`infinity simulate`)

Deterministic for a given seed. It produces genomes with known truth.

- **Reference:** a FASTA supplied by the user, or a random sequence of chosen length.
- **Clusters:** `m` cluster roots, each derived from the reference by a chosen divergence (default 2%).
- **Tree per cluster:** a random binary tree with `n` leaves. SNPs per branch drawn from a Poisson distribution with a chosen mean. Every SNP hits a position not mutated before, so true distances are simply path lengths.
- **Optional:** small indels per branch (length 1-10); gain and loss of segments (500-5000 bp); cutting each genome into `c` contigs in random order and orientation.
- **Outputs:** one FASTA per leaf; the true tree (Newick, branch lengths in SNPs); the true pairwise SNP matrix; the true cluster of each genome; the true list of mutations per leaf.

### 10.2 Properties that must always hold (automated tests)

1. A genome and its reverse complement give the same set of anchor hashes and the same set of block hashes.
2. A genome compared with itself, or with its reverse complement, or with a copy whose contigs are shuffled and flipped, gives zero variants and 100% matched blocks.
3. One SNP introduced into a genome gives exactly one SNP, at the right position, with the right bases.
4. On a random sequence, the number of level-j anchors is within 10% of `length / densities[j]`.
5. Same input and parameters give byte-identical output with 1 thread and with many.
6. Flat mode, hierarchical mode, and every cache switch give identical variants, distances, and trees.

### 10.3 Truth tests

| Test | Simulation | Must hold |
|---|---|---|
| T1 | 2 Mb random reference, 1 cluster, 50 leaves, mean 3 SNPs per branch | Pairwise SNP matrix equals truth exactly. Differing positions equal truth exactly. |
| T2 | T1 plus mean 0.5 indels per branch | At least 99% of pairs exact; no pair off by more than 2. |
| T3 | T2, each genome cut into 100 contigs, shuffled and flipped | No pair's distance exceeds truth by more than 2. Sum of reported distances at least 99% of the true sum. |
| T4 | T2 plus segment gain and loss | T2 thresholds still hold. Each gained segment appears as an insertion or novel bp, each lost one as a deletion or uncovered interval, with length within 5%. |
| T5 | 5 clusters of 20, 2% apart | Clusters recovered exactly. |
| T6 (optional) | As T3 on a real bacterial genome supplied by the owner, to exercise repeats | At least 99% of true SNP differences found; at most 0.1% false. |

### 10.4 Performance targets

These are targets to report against, not gates. Machine: 8 cores, SSD. Data: 1,000 simulated genomes of 5 Mb in 20 clusters.

- First blurry tree: within 60 s.
- All clusters exact: within 5 min.
- Peak memory: under 4 GB.

---

## 11. The experiment: does it pay?

This section is a deliverable, not an appendix **(owner's purpose)**.

**Switches.** Each one turns off one part of the idea while leaving results identical.

| Switch | What it removes |
|---|---|
| `--no-block-skip` | Block hashing altogether; every stretch between shared anchors is compared directly. |
| `--flat` | Nested levels; only leaf blocks are used. |
| `--no-piece-cache` | Resolving a shared difference once. |
| `--no-profile-grouping` | Handling identical genomes once. |
| `--single-cluster` | Cluster-first; all genomes are compared to one pivot. |

**Measurements**, from `run.json` and `events.jsonl`: time to first tree; time to exact tree; peak memory; counters; accuracy against truth (distance error, and Robinson-Foulds distance between each emitted tree and the true tree, after collapsing branches shorter than half a SNP).

**The key figure:** tree accuracy plotted against elapsed time, one curve per configuration. It shows the "photo coming into focus" as a measurement.

**Datasets.**
1. Simulated, nearly identical (T1-style, scaled to 100, 1,000, and 5,000 genomes).
2. Simulated with indels, segment gain and loss, and fragmentation.
3. A real outbreak dataset chosen by the owner. A suggestion to check: the public benchmark datasets for phylogenomic pipeline validation by Timme and colleagues (2017).

**Outside baselines** (a separate script, Linux only, not part of the program): Mash plus neighbour-joining as the fast blurry reference; SKA2 as the exact alignment-free reference.

**`infinity bench`** compares an output folder with a truth folder and writes the measurements as JSON. A script runs the full matrix and writes `REPORT.md` with tables and the accuracy-over-time figure.

**Proposed verdict rule (author's call).** A part of the idea "pays" if turning it off makes time-to-exact or peak memory at least 20% worse on at least two of the three datasets. Report every result, including the ones where it does not pay.

### Experimental extension: the k dimension

Build this last. Its payoff is the least certain part of the idea, and the report should say plainly what was found.

- **Sketch level.** Compute sketches at several k (for example 15, 21, 31). For a genome and a candidate pivot, with `p_k = 2J_k / (1 + J_k)`, fit `ln p_k = ln(1 - a) + k * ln(1 - pi)` by least squares. `pi` estimates SNP divergence and `a` gene-content difference. Test whether clustering on `pi` matches the truth better than single-k clustering.
- **Anchor level.** Run the comparison at more than one k (for example 21 and 31). For pieces that end as breakpoints or unplaced at one k, retry with the other. Measure whether the comparable fraction or the accuracy in repeats improves on real data.

---

## 12. Technology and layout **(author's call)**

- **Language:** Rust (stable).
- **Workspace:**

```
infinity/
  Cargo.toml              workspace
  crates/
    infinity-core/           engine library: scan, sketch, cluster, compare, tree, events
    infinity-sim/            simulator and truth comparison
    infinity-cli/            binary "infinity"
    infinity-gui/            binary "Infinity-Win"
  scripts/                benchmark matrix, outside baselines, report
  tests/                  integration tests
  DECISIONS.md            deviations from this outline
  README.md               user guide
```

- **Suggested libraries (verify first):** `needletail` (FASTA, gzip), `rayon` (parallelism), `xxhash-rust` (xxh3), `rustc-hash` or `ahash`, `dashmap`, `crossbeam-channel`, `clap`, `serde` and `serde_json`, `eframe`/`egui` (GUI), `rfd` (file dialogs), `bio` (pairwise alignment) or a small in-house aligner.
- Prefer pure-Rust dependencies so that the Windows build needs no extra toolchain.
- Neighbour-joining, midpoint rooting, Newick reading and writing, and Robinson-Foulds distance are small enough to implement in the project.
- Targets: Windows x86-64 and Linux x86-64. Deliver one executable per program per platform.

---

## 13. Build order

| # | Milestone | Done when |
|---|---|---|
| M0 | Workspace skeleton, formatting and lint checks, test runner, `infinity --version`. | Builds and tests run on Linux and Windows. |
| M1 | Scan (6.2), blocks (6.3), sketch, and the simulator (10.1). | Properties 1, 4, 5 pass. |
| M2 | Clustering and backbone tree (6.4), event stream (7), `infinity run` writing `genomes.tsv`, `clusters.tsv`, `tree.nwk`, `events.jsonl`. | T5 passes. |
| M3 | GUI skeleton: add genomes, run, live genome list, clusters, backbone tree drawn and updating. | A 100-genome simulated run shows a live blurry tree; the window never stalls. |
| M4 | Exact engine in flat mode (6.5-6.9): variants, coverage, profiles, distances, cluster trees, full tree. | Properties 2, 3 and tests T1-T4 pass. |
| M5 | GUI complete (section 8): focus styling, smooth updates, map panel, parameters and counters, export, open results. | Every item in section 8 can be demonstrated on a simulated run. |
| M6 | Staged Focus 1 then Focus 2, hierarchical mode, caches, memory budget, all switches (6.10, 6.11), `infinity bench`. | Property 6 passes. Performance numbers reported against 10.4. |
| M7 | Experiment (section 11): matrix script, outside baselines, `REPORT.md`. | Report exists with all configurations and datasets available. |
| M8 | Experimental k dimension (section 11). | Findings added to `REPORT.md`. |
| M9 | Release builds for both platforms, user guide in `README.md`. | A new user can install, load genomes, and read the tree from the guide alone. |

---

## 14. Who decided what

**Decided by the owner**

- The core idea: coarse-to-fine, multi-level, position-carrying comparison (section 1).
- Target: divergence of nearly identical genomes, to single-SNP level, in minutes.
- Cluster first; fine comparison inside clusters only.
- A visual, live interface showing the tree and the parameters.
- A desktop program for Windows and Linux.
- The project must show whether the approach pays.

**Decided by the author of this outline (owner may override)**

1. Design target is bacterial-size assemblies, tens to a few thousand genomes.
2. Engine approach A: anchors, nested blocks, one pivot per cluster (section 4.1).
3. Rust with a native GUI (section 4.2).
4. Default parameters in 6.1, including k = 31 and cluster distance 0.005.
5. The pivot is the most contiguous assembly of its cluster.
6. Neighbour-joining trees, midpoint rooted.
7. The window layout in section 8.
8. The verdict rule in section 11.
9. The k dimension is built last, as an experiment.
10. The working name "M-Tato" (since renamed to Infinity-Win by the owner).

**Open questions for the owner**

1. Which organisms and how many genomes are the first real test? This fixes the real dataset of section 11 and confirms or changes decision 1.
2. Is the window layout in section 8 what you have in mind?
3. Under which licence should the code be released?
4. Should reads (FASTQ) ever be accepted, or are assemblies enough?

---

## 15. Known limitations of this design

- **Pivot bias.** Sequence that the pivot lacks is not compared between members. It is reported as novel bp, not as distance. A later version could add such sequence to a growing cluster reference.
- **Repeats.** Blocks that occur more than once in the pivot can be mapped to the wrong copy. Variants there are flagged, not corrected.
- **Divergent clusters.** The block shortcut stops helping once most 256-bp blocks differ (around 0.5% divergence). The anchor path still works but more slowly, and above a few percent anchors become too sparse. Keep `cluster_dist` small.
- **Fragmented assemblies.** Every contig end costs a little coverage, even with end extension.
- **Trees.** Neighbour-joining on SNP counts ignores recombination and is not a substitute for a full phylogenetic analysis.
- **Novelty.** The prior-art check in section 3 was brief. A proper literature review should precede any publication.
