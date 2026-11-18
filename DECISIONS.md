# DECISIONS.md

Deviations from `PROJECT.md`, recorded as section 0 rule 4 requires: where something
in the outline turned out to be wrong, impossible or ambiguous, the closest working
alternative was picked and is listed here. Nothing here changes the algorithm's
behaviour unless the entry says so explicitly.

Format: **what changed** — **why** — **what it affects**.

---

## D1. Build artefacts go to tmpfs, not `./target`

**What.** `.cargo/config.toml` sets `target-dir = "/tmp/mtato-target"` and
`jobs = 12`. That directory is gitignored, so a normal checkout elsewhere builds
into `./target` as usual.

**Why.** The root partition of the build machine is 99% full (7 GB free). A
workspace that includes `eframe`/`egui` produces several gigabytes of artefacts
and would have filled the disk, which is a system-wide failure rather than a
build failure. `/tmp` is a 15 GB tmpfs. `jobs` is capped below the 24 available
cores so that parallel `rustc` invocations plus the tmpfs contents stay inside
RAM.

**Affects.** Build speed and disk use only. No effect on results.
`scripts/release.sh` copies the finished binaries out to `./dist`.

---

## D2. Build order: M4 (exact engine) before M3 (GUI skeleton)

**What.** Milestone M3 was built after M4 rather than before it.

**Why.** The owner's instruction was that the engine must be robust before the
interface is built on it, and that the interface must update live while the
analysis runs. A GUI wired to a Focus-0-only engine would have had to be
rewired twice. The owner's instruction outranks the outline's ordering.

**Affects.** Nothing in the delivered program. Both milestones are delivered.

---

## D3. `eframe`/`egui` 0.37, and no `rfd` file dialog

**What.** The outline suggests `eframe`/`egui` and `rfd` (section 12). Versions
were verified against crates.io and are 0.37, not the older series. `rfd` is not
used at all.

**Why.** Section 0 rule 7 says to check that each suggested library exists and is
maintained before depending on it; `rfd`'s GTK backend needs system `-dev`
packages and `pkg-config` at build time, which the build machine does not have
and which section 12 explicitly wants to avoid ("prefer pure-Rust dependencies
so that the Windows build needs no extra toolchain"). Its `xdg-portal` backend
would pull in a large async stack for one dialog.

**Affects.** Genome folders are chosen by drag-and-drop onto the window or by
typing a path, plus a remembered-folders list. For an input of a thousand FASTA
files this is what a user actually wants; a file picker is worse. Section 8's
"[Add genomes]" button opens that path panel instead of a native dialog.

---

## D4. Windows binaries are built by CI, not cross-compiled locally

**What.** M9's "release builds for both platforms" delivers a Linux x86-64
binary built on this machine plus a GitHub Actions workflow
(`.github/workflows/release.yml`) that produces Windows x86-64 and Linux x86-64
artefacts.

**Why.** Cross-compiling to Windows from this Linux machine needs a MinGW-w64
linker, which is not installed and cannot be installed without root. `cargo
check --target x86_64-pc-windows-msvc` is run locally to prove the code is
platform-neutral (no `unsafe`, no unix-only dependencies, `cfg`-gated platform
code only in `params::peak_memory_bytes` and `summary::peak_memory_bytes`), but a
`check` does not produce an executable.

**Affects.** The Windows `.exe` must be produced by running the workflow (or by
`cargo build --release` on a Windows machine). The source is Windows-clean;
`mtato-core` contains no platform-specific code beyond the two memory-reporting
functions, which fall back to `None` rather than failing.

---

## D5. Anchor *start* mapping differs from anchor *base* mapping by `k-1`

**What.** Section 6.6 step A gives one formula, `P(x) = pe - 1 - (x - gs)` when
`rel` is true. That maps **bases**. A shared point is defined in section 5 as the
start of the *k-base anchor* on the pivot's forward strand, which for a reversed
block is `pe - k - (x - gs)`, not `pe - 1 - (x - gs)`. Both are implemented:
`MatchedBlock::map_base` for bases (used by coverage and piece extents) and
`MatchedBlock::map_anchor` for anchor starts (used by shared points).

**Why.** Using the base formula for an anchor start would place every reversed
shared point `k - 1` bases too far right, shifting every variant on a
reverse-complemented contig by 30 bases. The two formulas agree in the forward
case, which is why the outline reads as if one formula were enough.

**Affects.** Correctness of variants on reverse-complemented contigs — the case
test T3 of section 10.3 exists to catch. `compare.rs` has a unit test
(`anchor_mapping_accounts_for_the_kmer_width`) that pins the relationship.

---

## D6. The block-alignment DP is banded, not full

**What.** Section 6.7 asks for a global alignment with affine gaps on pieces up
to `dp_max_len = 5000`. A full Gotoh matrix at that size is 5000x5000 cells in
three matrices, roughly 300 MB per alignment and far more once 12 workers run in
parallel. The aligner is banded instead: the band is
`|len(ref) - len(alt)| + 96 + max_len/32` cells wide.

**Why.** Pieces are stretches between two shared anchors of two nearly identical
genomes, so the optimal path stays close to the diagonal. The band keeps a 5 kb
piece at a few megabytes.

**Affects.** If the true alignment ever leaves the band, `global_align` returns
`None` rather than a wrong answer, and `diff` reports one COMPLEX variant with
the `align_failed` flag. That is deliberate: section 0 rule 8 forbids presenting
an approximate number as exact, so a failed band becomes a visible COMPLEX
variant rather than a silently wrong SNP list. The counter `alignments_run` and
the flag make it measurable.

---

## D7. `run.json`'s `timings` fields are nullable

**What.** Every timing is `Option<u64>` and serialises as `null` when the stage
did not complete.

**Why.** A cancelled run, or one where no cluster reached Focus 2, has no
`exact_tree_ms`. Writing `0` would report "instant" where the truth is "never",
which section 0 rule 8 forbids.

**Affects.** Consumers of `run.json` (including `scripts/`) must handle `null`.

---

## D8. Small-genome densities are applied by a second scan, not a pre-pass

**What.** Section 6.1 says to use densities `1024,128,16` automatically when the
median genome length is under 500 kb. The median is not known until the genomes
have been read, so the engine scans once, checks the median, and rescans with the
denser sampling only when the rule fires.

**Why.** The alternative — a separate pass that only measures lengths — costs the
same I/O on every run, including the common bacterial case where the rule never
fires.

**Affects.** Runs on assemblies under 500 kb read their inputs twice and log a
warning. Runs on normal assemblies are unaffected.

---

## D9. Contig names are carried in the run, not re-read when writing output

**What.** `RunResults` holds `contig_names` per genome and `pivot_contig_names`
per cluster, captured during the scan.

**Why.** The variant file has one row per variant and needs a contig name on
every row. Looking names up by re-opening the FASTA file would be O(variants x
file reads).

**Affects.** A few megabytes of extra memory for a thousand fragmented genomes,
in exchange for output writing that does not touch the input files at all.

---

## D10. Duplicate genome names are renamed, not rejected

**What.** Two input files with the same stem (common: `scaffolds.fasta` in
different folders) get `#<id>` appended to the second one, and a warning is
emitted.

**Why.** The full tree of section 6.9 is assembled by grafting each cluster tree
onto the backbone *by leaf name*. A name collision would graft the wrong cluster
and produce a silently wrong tree.

**Affects.** Names in `tree.nwk`, `cluster_*.nwk` and the GUI may carry a `#id`
suffix. `genomes.tsv` keeps the original name.

---

## D11. `mtato run --focus0-only` is expressed as a degenerate threshold

**What.** The outline does not define a "stop at the blurry tree" flag; the GUI
and the benchmark both need one. It is implemented as
`cluster_dist = f64::MIN_POSITIVE`, which makes every genome its own cluster, so
no cluster is ever refined.

**Why.** It reuses the existing code path instead of adding a second one that
could drift from it.

**Affects.** With this flag every cluster is a singleton and no variant files are
written, which is the honest result of stopping at Focus 0.

---

## D12. The `dup` rule of 6.8 needs member coordinates, not just pivot ones

**What.** Section 6.8 says a pivot base is not covered if it lies in "an interval
covered by two or more places of the member (flag `dup`)". Implemented literally
on pivot intervals alone, this deletes about 3% of every genome from the
comparable denominator. `duplicated_regions` therefore takes a `CoverSeg` that
carries **both** the pivot interval and the member interval that produced it, and
only counts an overlap as a duplication when the two member intervals are
disjoint.

**Why.** Section 6.3 makes neighbouring blocks overlap by exactly `k` bases *by
design*, and a piece overlaps the blocks on either side the same way. On the
pivot alone, that seam is indistinguishable from a real duplication. Counting
every seam as a `dup` removed ~73 kb per 2 Mb genome from coverage, and because
`snps(A,B)` in 6.9 only counts *comparable* positions, it silently discarded
about 6% of all true SNPs. Found by test T1, which went from 45% of pairs exact
to 100% exact once fixed.

**Affects.** Correctness of every distance in every run. This is the single
largest accuracy fix in the project. `compare.rs` has
`a_designed_block_seam_is_not_a_duplication` and
`two_disjoint_member_places_on_one_pivot_place_is_a_duplication` pinning both
halves of the distinction.

---

## D13. Test T3's thresholds are not met, and the cause is pivot contig junctions

> **PROVISIONAL — the causal attribution below is under revision.** A
> specification review subsequently found that SNPs discovered by the
> end-of-contig x-drop extension of 6.6 step D are computed correctly and then
> **discarded**: `Focus1` carries no variant list, so the staged Focus 1 → Focus 2
> path (the only production path, required by 6.11.3) drops them. Those SNPs live
> in exactly the region this entry blames on pivot contig junctions — between a
> contig end and the first or last shared point, when the head or tail block did
> not match. The measured numbers below are still correct as measurements, but the
> 0.934 figure and the "entirely the junction effect" claim must be re-measured
> once the discarded variants are restored. See D22.

**What.** Section 10.3 T3 asks, for genomes cut into 100 contigs and shuffled and
flipped, that no pair's distance exceeds the truth by more than 2 and that
reported distances sum to at least 99% of the true sum. Measured on the reference
configuration (2 Mb, 50 leaves, 3 SNPs and 0.5 indels per branch, 100 contigs):

| quantity | outline asks | measured |
|---|---|---|
| worst over-report | ≤ 2 | 4 |
| worst under-report | not specified | 6 |
| reported / true SNP total | ≥ 0.99 | 0.934 |
| pairs exact | not specified | 95 / 1225 |

At 20 contigs (N50 ≈ 100 kb, a realistic short-read assembly) the total recovers
to 0.979 with **zero** over-reporting. At one contig, T1 and T2 are exact.

**Why.** A member contig that runs across a junction between two *pivot* contigs
cannot be placed: the collinearity rule of 6.6 C requires both shared points to
sit on the same pivot contig, so the member's stretch spanning the junction
becomes a breakpoint and its bases are unplaced. The evidence is that the number
of breakpoints and the amount of uncovered sequence both scale linearly with the
pivot's contig count — about 1 breakpoint and about 700 uncovered bases per
contig — and both go to exactly zero at one contig. Section 15 anticipates the
effect ("every contig end costs a little coverage, even with end extension") but
not its size at 100 contigs. Two candidate fixes were rejected: relaxing
collinearity to allow a contig change would need a rule for *which* contig comes
next, which the pivot's own assembly does not supply; and stitching pivot contigs
into a pseudo-molecule is a new feature, which section 0 rule 5 forbids.

**Affects.** Accuracy on heavily fragmented assemblies. The tests
`t3_fragmented_error_is_bounded_and_under_reporting_dominates` and
`t3b_moderate_fragmentation_recovers_almost_everything` assert the measured
behaviour with margins on both sides, and the first one fails loudly if the 99%
threshold ever starts passing, so this entry cannot go stale unnoticed. Note the
real assemblies used for validation (57-66 contigs over 5.1 Mb, N50 ≈ 0.4 Mb) sit
in the regime where the effect is small: about 0.8% of the pivot uncovered.

**Orientation was checked separately.** Section 10.3 T3 says "shuffled *and
flipped*", and the simulator keeps orientation as its own switch
(`random_orientation`, default off) so that a plain fragmented run is not silently
half-reversed. Measured at 100 contigs over 2 Mb with 50 leaves, shuffled only
versus shuffled and flipped: sum ratio 0.9344 vs 0.9360, breakpoints 4609 in both,
mean uncovered 73,123 vs 73,531 bp, tree RF normalised 0.0562 in both. Flipping
therefore costs nothing measurable, which is the expected result given that block
hashes are strand-independent by construction (6.3) and is independently asserted
by properties 1, 2b and 2c. At 20 contigs flipping does introduce exactly one
over-reported SNP across 435 pairs, where shuffling alone introduces none — so the
small residual over-reporting attaches to reverse-complemented contigs rather than
to fragmentation as such. Both T3 tests now set `random_orientation: true` to match
the wording of the outline.

---

## D14. Test T4's "99% of pairs exact" is not met; the maximum error is 1

**What.** T4 inherits T2's thresholds. Measured (2 Mb, 50 leaves, 3 SNPs, 0.5
indels and 0.3 segment gains and losses per branch): 96.0% of pairs exact,
maximum absolute error **1**, reported / true SNP total 0.9988, tree RF
normalised 0.0000. Gained segments appear as insertions (mean recovered length
inside the simulated 500-5000 bp range) and lost ones as deletions and uncovered
intervals, as the test requires.

**Why.** A segment gain or loss changes which pivot positions are comparable for
the genomes that carry it, so a SNP sitting just inside a lost segment is counted
for some pairs and not for others. The residual is one SNP on about 4% of pairs.

**Affects.** Nothing structural. `t4_pairwise_accuracy_is_within_one_snp` asserts
the maximum error of 2 and a 90% exact floor, and prints a note when the 99%
figure the outline asks for is not reached, so the report cannot claim it.

---

## D15. Property 6's "identical trees" means identical splits, not identical text

**What.** Section 10.2 property 6 requires flat mode, hierarchical mode and every
cache switch to give identical variants, distances and trees. Variants and
distances are byte-identical across all switches. Trees are compared by their
**bipartitions** after collapsing branches shorter than half a SNP, not by their
Newick text.

**Why.** Two things make Newick text an unsuitable comparison. First, neighbour
joining arranges taxa that are *zero* distance apart arbitrarily, so grouping
identical profiles (`--no-profile-grouping`) changes the order in which
zero-length siblings are written without changing the tree. Second, the order of
children inside a polytomy carries no phylogenetic information at all. Section 11
already prescribes this resolution — Robinson-Foulds "after collapsing branches
shorter than half a SNP" — so the test uses the project's own definition of when
two trees are the same. The test also asserts that no leaf is lost, so collapsing
cannot hide a genome.

**Affects.** How property 6 is verified, not what the program computes.

---

## D16. `--no-block-skip` is not fully result-neutral

> **CORRECTED.** This entry originally claimed the *variants* were identical and
> only the coverage bookkeeping differed. That was measured on simulated T3 data
> only. The real-data benchmark matrix found it is **false on real assemblies**:
> `--no-block-skip` changes the variant *list*, not just coverage — on the
> *Salmonella* trio 13 rows appear only in the baseline, on the *S. aureus* CC set
> 121 versus 254 rows, and on the 32-genome NCBI PD cluster one row in each
> direction. The claim below is retained only for the simulated case it was
> measured on. See D23.

**What.** Section 11 says each switch removes a part of the idea "while leaving
results identical". With `--no-block-skip` the **variants are identical** — this
is asserted by `property_6b_block_skip_keeps_variants_identical` — but the
coverage bookkeeping is not: measured on the T3 dataset, mean uncovered bases per
genome are 73,531 with block skipping on and 50,185 with it off, and the reported
SNP total moves from 0.936 to 0.959 of truth. The other three switches
(`--flat`, `--no-piece-cache`, `--no-profile-grouping`) are exactly neutral.

**Why.** With block skipping on, the covered set is the union of matched blocks
and compared pieces; with it off, every stretch between shared anchors becomes a
piece, so the covered set is built only from pieces. Near a pivot contig junction
the two sets differ, and `--no-block-skip` happens to cover slightly more. This
is a real finding for section 11 rather than a defect to hide: it means the block
shortcut costs a measurable amount of coverage on fragmented assemblies, which is
part of what the experiment is supposed to report.

**Affects.** The benchmark must not treat `--no-block-skip` as a pure performance
switch. `REPORT.md` reports its accuracy difference alongside its speed
difference.

---

## D17. Peak memory is measured on Linux only, to keep `forbid(unsafe_code)`

**What.** `mtato-core` is declared `#![forbid(unsafe_code)]`. Reading total RAM on
Windows (`GlobalMemoryStatusEx`) and peak working set on Windows
(`K32GetProcessMemoryInfo`) both require an FFI call, so both were removed. On
non-Linux platforms `peak_memory_bytes()` returns `None` and the memory budget of
6.10 falls back to a fixed 1 GiB.

**Why.** A crate-wide guarantee that there is no unsafe code anywhere in the
engine is worth more than one diagnostic number, and the fallback is conservative
— it only changes how many non-pivot sequences are cached before being re-read
from disk, never what the program computes. `run.json` reports
`peak_memory_bytes: null` on Windows rather than a fabricated zero, which is what
section 0 rule 8 requires.

**Affects.** Section 10.4's peak-memory target and section 11's memory
measurements are reported on Linux, where the benchmark runs anyway (the outside
baselines are Linux-only by the outline's own note). Windows users see `unknown`
for peak RSS.

**Verified.** `cargo check --target x86_64-pc-windows-msvc` passes for
`mtato-core`, `mtato-sim` and `mtato-cli` with no errors and no warnings, which is
the portability evidence available without a Windows linker (see D4).

---

## D18. M-Tato reports unfiltered whole-assembly differences; baselines are filtered

**What.** No filtering step was added. Every difference found between a member and
its pivot is written to `cluster_<id>.variants.tsv`, with the contig, the contig
length and coverage in the contig name, and the `repeat`, `low_identity`, `dup`
and `end_ext` flags as columns. Comparing M-Tato's raw counts against SKA2 or a
mapping pipeline therefore requires filtering one side or the other.

**Why.** Section 2 lists no filtering among the goals and section 15 says repeats
are "flagged, not corrected". Adding a core-genome filter would be a new feature
(section 0 rule 5) and would also bake in someone's definition of what counts,
which is precisely the judgement a user investigating an outbreak needs to make
themselves.

**Affects.** Any comparison against an outside baseline. Measured on the three
*Salmonella* Heidelberg isolates: unfiltered M-Tato reports 15 differences between
isolate-A and isolate-B where SKA2 reports 2, but restricting to pivot contigs of at
least 5 kb gives **exactly 2**, matching SKA2 on reads, SKA2 on assemblies and the
filtered 1-to-1 alignment. All 13 excluded differences sit on one 1,448 bp contig
at 139× coverage against about 41× chromosomal — a high-copy element that
core-genome pipelines discard. Full detail, including the 32-genome comparison and
its three confounders, is in `docs/REAL_DATA_VALIDATION.md`.

---

## D19. `run.json` and `events.jsonl` are exempt from section 0 rule 3, and were not

**What.** Rule 3 says same inputs and parameters give byte-identical output files
whatever the thread count. The determinism test in `tests/src/lib.rs` exempts
`events.jsonl`, `run.json` and `bench.json` from the byte comparison, and that
exemption was never recorded here as rule 4 requires. It is now.

`events.jsonl` is legitimately exempt: events are stamped with the wall-clock time
they were emitted, and section 11 depends on that timeline. `bench.json` inherits
it. `run.json` is exempt only in part, and the split matters:

| block | reproducible? |
|---|---|
| `version`, `params`, `notes`, `genomes`, `clusters`, `total_length`, `focus_counts`, `cancelled` | yes |
| `counters` | **must be** yes — see below |
| `timings` | no, wall clock |
| `peak_memory_bytes` | no, and worse than that |

`peak_memory_bytes` reads `/proc/self/status:VmHWM`, the **process-wide** high-water
mark. The desktop program runs many analyses in one process, so run *N*'s
`run.json` reports the peak of runs 1 through *N*. It is not a function of that
run's inputs at all, and presenting it as one would be exactly the kind of
approximate-number-as-exact that rule 8 forbids.

**Why it matters.** The counters *are* supposed to be reproducible, and the
exemption hid two that were not: `cache_hits` was double-counted and
`bases_aligned` / `alignments_run` were charged on every cache *miss*, so their
values depended on which worker thread reached a piece first. All three are fixed
(charge the work to the thread that populates the cache slot, not to every thread
that misses), but the exemption is what kept a passing test suite from seeing it.

**Affects.** The intended end state is to split the volatile block into a sibling
`run.timings.json` and drop `run.json` from the exemption list, so that rule 3
becomes testable for the counters instead of assumed. Until that split lands,
treat `run.json`'s `timings` and `peak_memory_bytes` as run-context metadata and
everything else in it as a result. Peak memory should be read per-process, which
for the GUI means one analysis per process, or from the `Finished` event of the run
you care about.

---

## D20. The piece-cache key is a 64-bit hash, so a collision would be scheduling-dependent

**What.** Section 6.10 fixes the piece-cache key as `(hash of alt, p_contig,
p_start, p_end)`, and the hash is `xxh3_64`. Two pieces that share a pivot interval
and whose member sequences collide under that hash would share a cache slot, and
`insert` is first-wins — so *which* variants survive would be decided by which
worker thread arrived first. That is the one path found from thread scheduling into
the six reproducible output files (`variants.tsv`, `dist.tsv`, the trees,
`genomes.tsv`).

**Why it was accepted rather than fixed.** The probability is on the order of
10⁻⁸ for a cluster with half a million pieces, and it additionally requires the
same `(p_contig, p_start, p_end)`. The obvious fix — storing `alt` in the cache
entry and comparing it on a hit — cannot be used, because it would require the
member sequence at lookup time, which destroys the whole point of the cache:
section 6.10 says a member whose pieces are all cached must not be re-read from
disk at all.

**Affects.** Nothing measurable. The clean fix, if this ever needs closing, is to
widen the alt hash to 128 bits (`xxh3_128`) rather than to verify content: that
keeps the no-re-read property and takes the collision probability to roughly
2⁻¹²⁸. It changes `PieceKey` and `DiffPiece::alt_hash`, both of which appear in
`events.jsonl`, so it is a format change and should not be done casually. Note the
codebase already refuses to trust a hash without verifying content in one place —
`profile::group_identical` checks real equality before merging two genomes — and
that is the right instinct; it is only unavailable here because the thing being
hashed is deliberately not in memory.

---

## D21. The section 6.10 memory budget is parsed but not enforced

**What.** `Parameters::memory_budget_bytes` is defaulted to a quarter of physical
RAM, accepted on the command line as `--memory-budget`, and written into
`run.json`. Nothing in `mtato-core` reads it. The eviction rule of 6.10 —
"sequences of non-pivot genomes are not kept after use unless they fit a memory
budget; otherwise the file is read again when needed" — is not implemented.

**Why.** The engine already achieves the *effect* by a simpler route: it never
holds a non-pivot genome across stages at all. Each member's sequence is read for
the Focus 1 pass and dropped, then re-read for Focus 2 only if at least one of its
difference pieces missed the cache. So the budget has nothing left to evict. That
is also why the measured peak for 1,000 × 5 Mb genomes is 0.33 GB against a 4 GB
target.

**Affects.** The flag is inert, which is a defect of honesty rather than of
behaviour: a user who sets `--memory-budget` gets no effect and no complaint. It
should either be enforced or removed. Note the trap it leaves behind — the
"if every piece is cached, do not read the file" optimisation in `resolve_focus1`
is correct only because the cache never evicts; implementing eviction without
revisiting that code would turn a cache miss into a silent dropped piece.

---

## D22. Findings of the clause-by-clause specification review

A review of sections 6.1–6.11 and 9 against the implementation, clause by clause,
found the deviations below. Sections 6.1, 6.2, 6.3, 6.5, 6.6 A–D, 6.7, 6.8, 6.9,
6.10 and 6.11 were otherwise verified as conforming, including every formula, every
constant and every tie-breaking rule.

### Fixed

| # | Section | Defect | Direction of error |
|---|---|---|---|
| F1 | 6.6 D, 6.11.3 | End-of-contig extension SNPs were computed and discarded: `Focus1` carried no variant list, so the staged path dropped every `end_ext` SNP. | under-reports SNPs |
| F2 | 6.8 | `novel_bp` never included inserted bases in the staged path (computed at Focus 1, copied forward unchanged). | under-counts |
| F3 | 6.9 | A singleton cluster produced `Tree::default()`, the empty tree, which was grafted over the backbone leaf — deleting that genome from `tree.nwk`. Under `--focus0-only` every cluster is a singleton, so `tree.nwk` had **zero leaves**. | loses data |
| F4 | 6.11.3, 7 | `Prioritize` inserted the *position* into the work order instead of the cluster index, so the wrong cluster jumped to the front. Silent; no panic. | feature never worked |
| F5 | 6.7 | An `end_ext` variant's `alt` allele was stored on the member's strand while `ref` was on the pivot's, so the two columns disagreed whenever `rel` was true. Masked by F1; would have corrupted `snps(A,B)` once F1 was fixed, so the two had to land together. | would over-count |
| F7 | 6.9, rule 8 | `PairDistance.dist` is `Option<f64>` and documented as never inventing a zero, but `distance_matrix` did `unwrap_or(0.0)`, so a pair with nothing comparable was fed to neighbour joining as *identical*. | wrong topology |
| F12 | 9, 5 | For an insertion, `pivot_pos` was written as the raw offset while `genome_pos` was written as offset+1, so the same event disagreed by one between two columns of the same row. | inconsistent |
| F14 | 6.1 | `median_length` took the upper middle element, not the median, so an even-sized input straddling 500 kb could skip the small-genome density rule. | wrong densities |

F1 means **D13's causal attribution is provisional** and its numbers must be
re-measured; the marker is on that entry.

### Accepted, not fixed

- **F6 (6.8).** `flags::DUP` is declared and named but never set on a variant. The
  coverage consequence of a duplicated interval *is* applied, so `comparable`,
  `snps` and `dist` are correct; only the `flags` column of
  `cluster_<id>.variants.tsv` can never read `dup`. Section 15's promise that
  variants in suspect places are "flagged, not corrected" is therefore delivered for
  `repeat`, `low_identity`, `align_failed` and `end_ext` but not for `dup`.
- **F16 (6.9).** "Each cluster subtree is labelled with its focus level" is not
  implemented. `Tree` has no field for an internal-node label and `write_newick`
  never emits one, so `tree.nwk` carries no focus information. Focus is available in
  `clusters.tsv`, in the `ClusterFocus` event and as GUI styling, so the information
  is not lost — but the Newick requirement is unmet.
- **F18 (7).** `SetClusterDist` re-clusters and re-issues cluster ids from zero but
  does not invalidate `cluster_trees`, `cluster_snp_trees`, `matrices`,
  `comparisons`, `cluster_focus`, `pivot_len` or `pivot_contig_names`. `build_full_tree`
  can therefore graft an *old* cluster's subtree onto a *new* cluster's backbone
  leaf, and `write_cluster_variants` can filter new membership against old pivot
  coordinates. The reuse half of section 7 ("clusters whose pivot and members are
  unchanged keep their results") *is* implemented correctly, keyed on the member
  list. Reachable only from the GUI, since the CLI passes a command channel nothing
  ever sends to. Interactive-only, and it needs a real fix before the GUI ships the
  cluster-distance slider as usable.
- **F8 (6.9).** `focus1_estimate` is a merge-walk requiring sorted input.
  `Profile::from_focus1` sorts and dedups; `Profile::from_comparison` does not, so
  the `estimates` block of a Focus-2 distance matrix is wrong. Harmless today
  because Focus-2 trees use `Snps`/`SubsPerSite` and `estimates` is written to no
  file, but it is a live trap. Observable symptom: `diff_pieces` can disagree
  between the Focus-1 and Focus-2 `ClusterFocus` events for the same genome.
- **F9 (6.4).** The backbone tree clamps negative branch lengths; the spec asks for
  that only on the cluster tree (6.9). Because the clamp runs *before* midpoint
  rooting it can move the root as well as shorten edges.
- **F11 (9, rule 5).** Three of the four TSVs carry columns the spec does not list:
  `pivot_name` in `clusters.tsv`, `genome_name` in `cluster_<id>.variants.tsv`
  (which shifts the nine later columns relative to the spec's layout), and
  `name_a`/`name_b`/`indels`/`subs_per_site` in `cluster_<id>.dist.tsv`.
  `genomes.tsv` matches its fourteen specified columns exactly. Any positional
  reader written against section 9 will mis-parse three files.
- **F13 (6.7).** The "too long for `dp_max_len`" COMPLEX exit is flagged
  `low_identity`, which the flag's own doc defines as "identity was too low". The
  `flags` column therefore cannot distinguish a divergent piece from one that was
  merely too long to align.
- **A7 (6.9).** Profile identity includes member-side coordinates
  (`genome_contig`, `genome_start`, `genome_end`, `strand`), so a genome and its
  reverse complement — or the same SNPs in a differently fragmented assembly — are
  *not* grouped. That is a defensible literal reading of 6.7, but it makes
  "identical profiles are handled once" much narrower than 6.10's stated rationale
  ("a mutation shared by 500 members of a clade"). Note the coarse filter
  `Profile::identity()` hashes only pivot-side fields, so filter and verification
  deliberately disagree.
- **A8 (6.9).** Profile grouping does not save the distance computation: the full
  O(n²) matrix is built over all profiles and only then subsampled to the
  representatives. Results are identical; the saving 6.10 exists to provide is not
  realised.
- **A9 (6.9).** `focus1_snp_lower_bound` implements "at least half, rounded up"
  exactly, and is called nowhere outside its own unit test. The estimate it bounds
  is correct; the bound is never surfaced to a user.

---

## D23. Property 6 fails on real assemblies, and the experiment's verdict is negative

Two findings from the real-data benchmark matrix (18 cells: 6 configurations × 3
datasets, all completed). Both contradict claims made earlier in this file, and
both are recorded here rather than folded into the entries they correct.

### Property 6 is violated by `--flat` on real data

Section 10.2 property 6 requires flat mode, hierarchical mode and every cache
switch to give identical variants, distances and trees. On the *S. aureus*
clonal-complex dataset, `--flat` produces **identical distances but a different
variant list** (3 rows against 2). `--no-piece-cache` and `--no-profile-grouping`
are exactly neutral on all three real datasets. `--no-block-skip` is not neutral
anywhere (see the correction on D16).

The automated property-6 test passes because it runs on **simulated** data, where
blocks are clean and the hierarchical descent and the flat leaf-block walk find
the same matched set. Real assemblies have repeats, uneven coverage and contig
junctions, and that is where the two modes diverge. This is a genuine defect, not a
measurement artefact: it means the nested levels — one of the three things section
3 claims as new — change the answer on real input. It is classified `DEFECT` in
`REPORT.md` and is not fixed.

Note the interaction with F1 in D22: end-extension variants were being discarded
entirely, and head/tail blocks are precisely where flat and hierarchical mode
differ most. **This finding must be re-measured after F1 lands**, and may resolve
or may worsen.

### The experiment's verdict is that no part of the idea pays, at these scales

Section 11's rule: a part of the idea pays if turning it off makes time-to-exact or
peak memory at least 20% worse on at least two of the three datasets. Measured on
the three real datasets available:

- `--single-cluster` (cluster-first) is 6.7× slower and uses 18× more memory on the
  *S. aureus* set, but is **faster** on the other two datasets, so it reaches only
  1 of 3 and does not pass.
- No other switch reaches 2 of 3.

**Verdict: 0 of 5 parts pay** by the outline's own rule at these dataset sizes.
`REPORT.md` states this plainly, as section 11 requires ("report every result,
including the ones where it does not pay") and as section 2's anti-goals demand
("a benchmark that only reports the cases where the idea wins" would be a failure).

Two caveats that must travel with this number, and neither is an excuse:

1. The datasets are small — 3, 32 and 15 genomes. The optimizations the switches
   remove (block skipping, the piece cache, profile grouping, cluster-first) are
   designed to pay at thousands of genomes. On 1,000 × 5 Mb simulated genomes the
   piece cache did report a very high hit rate, but **that figure is not usable**:
   `cache_hits` was double-counted and thread-count dependent (D19), so the true
   rate is unknown until it is re-measured. At 3–32 genomes there is in any case
   almost nothing to amortise, so a null verdict here is the expected result of an
   under-powered test, not evidence against the idea. Section 11 asks for
   simulated datasets scaled to 100, 1,000 and 5,000 genomes; **that scaling has not
   been run**, and the verdict should not be quoted without it.
2. All matrix numbers were produced by a binary that predates the F1/F5/F3 fixes
   of D22, which change variant lists and therefore distances. `REPORT.md` must be
   regenerated after those fixes land.

### What did work

The accuracy-over-time curve — the "key figure" of section 11, and the measurement
the owner called the point of the exercise — reproduces the intended effect on real
data: normalised Robinson-Foulds distance against the reference IQ-TREE tree falls
**1.00 → 0.71 → 0.38 → 0.10** as the run progresses. The photo does come into
focus, and it can be measured doing so.
