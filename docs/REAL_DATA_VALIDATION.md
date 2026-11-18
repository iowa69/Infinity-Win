# Validation on real data: three *Salmonella* Heidelberg outbreak isolates

This is the acceptance check the project owner asked for: three similar
*Salmonella* genomes where two are close and one carries 2 SNPs.

**Data.** `<salmonella-data>/asm/`
— `isolate-A.fasta`, `isolate-B.fasta`, `isolate-C.fasta`. Each about 5.13 Mb,
57–66 contigs, 51.8% GC. Assemblies and the reference distances below come from
the existing analysis in that folder (`analysis/results.json`), not from Infinity-Win.

**Command.**

```sh
infinity run $ASM/isolate-A.fasta $ASM/isolate-B.fasta $ASM/isolate-C.fasta \
          -o out --threads 8
```

**Run.** 3 genomes → 1 cluster, all at Focus 2. First blurry tree 0.20 s, exact
tree 0.57 s, peak RSS 38 MB. Pivot chosen automatically: `isolate-B`
(best N50, 442,493).

---

## Result

| Pair | SKA2 (reads) | SKA2 (assemblies) | 1-to-1 alignment, raw | 1-to-1, filtered | **Infinity-Win, all variants** | **Infinity-Win, contigs ≥ 5 kb** |
|---|---|---|---|---|---|---|
| isolate-A – isolate-B | 2 | 2 | 32 | 2 | 15 | **2** |
| isolate-A – isolate-C | 2 | 2 | 16 | 2 | 6 | **2** |
| isolate-B – isolate-C | 0 | 2 | 17 | 0 | 0 | **0** |

**Filtering Infinity-Win's variants to pivot contigs of at least 5 kb reproduces the
published chromosomal SNP counts exactly.** The two SNPs it reports are:

```
isolate-B_contig4_len442493_cov42.0   pos 175934   G>T
isolate-B_contig8_len231110_cov40.4   pos 156784   C>A
```

## Where the unfiltered extra differences are

All 13 extra differences between isolate-A and isolate-B sit on a single pivot
contig:

```
isolate-B_contig37_len1448_cov139.5    1448 bp    14 variants (13 snp + 1 complex)
```

That contig is **1,448 bp at 139.5× coverage**, against about 40–42× for the
chromosomal contigs — roughly 3.4× the chromosomal depth. It is a high-copy
element, and its variants include one `complex` flagged `low_identity`, i.e. a
stretch too divergent to resolve into individual events. Every reference
pipeline in the comparison excludes this kind of region: SKA2 by restricting to
the shared k-mer core, the alignment pipeline by its explicit filtering step
(32 raw → 2 filtered).

Infinity-Win does not filter, by design: section 2 lists no filtering among its goals,
and section 15 states that repeats are flagged rather than corrected. What it
does instead is give you everything needed to filter it yourself — the contig,
the contig length and coverage are in the contig name, and the `repeat`,
`low_identity`, `dup` and `end_ext` flags are columns in
`cluster_<id>.variants.tsv`.

## What this does and does not show

**Shows.** On real assemblies the engine finds the real chromosomal SNPs at the
right positions with the right bases, agrees with two independent methods once
the same filtering is applied, and costs well under a second for three 5 Mb
genomes. It also shows something the reference numbers hide: the two isolates
differ by 13 additional variants on a high-copy element, which a core-genome
filter discards silently.

**Does not show.** Infinity-Win is not a drop-in replacement for a filtered
core-genome SNP pipeline. Its raw counts are whole-assembly counts and are
larger than a filtered core-genome number. Anyone comparing Infinity-Win output to
SKA2 or to a mapping pipeline must filter both, or filter neither.

---

## Larger real dataset: the 32-genome NCBI Pathogen Detection cluster

Cluster PDS000037076, 32 assemblies at
`<salmonella-data>/ncbi_pd/asm/`.

**Run.** 32 genomes → 1 cluster at `d ≤ 0.005`, all at Focus 2. First tree
0.46 s, exact 1.04 s, peak RSS 272 MB, 92% of block lookups identical, piece
cache hit 508 of 754 times.

**Against SKA2 on reads** (496 pairs, all joinable):

| measure | value |
|---|---|
| exact agreement | 22 / 496 (4.4%) |
| mean absolute difference | 11.6 SNPs |
| median absolute difference | 6 SNPs |
| SKA2 range | 0–46 |
| Infinity-Win range | 0–145 |
| sum ratio (Infinity-Win / SKA2) | 1.618 |
| Pearson r | 0.44 |

This disagreement is **not** a like-for-like comparison and should not be read as
an accuracy figure. Three confounders, in order of importance:

1. **Reads versus assemblies.** The SKA2 baseline is `cluster_reads_distance.tsv`
   — computed from reads. Infinity-Win compares assemblies. The trio above shows the
   same effect in miniature: SKA2 on reads gives 2/2/0 while SKA2 on the
   assemblies gives 2/2/2.
2. **Mixed assemblers.** The three outbreak isolates were assembled locally; the
   other 29 are NCBI SKESA assemblies. Assembly method affects the accessory
   content and the error profile, and a whole-assembly comparison sees both.
   QUAST in the source analysis reports about 41 mismatches per 100 kbp and
   95.5% genome fraction for these assemblies, i.e. they are not error-free.
3. **No filtering on either side of the comparison.** Infinity-Win reports 776 SNPs,
   48 deletions, 26 insertions and 11 complex variants over 31 genomes, with 849
   of 861 variants carrying no flag at all — so this is not predominantly a
   repeat artefact. It is a whole-assembly count against a core-genome count.

A fair baseline for section 11 needs SKA2 run **on these same assemblies**, which
is what `scripts/baselines.py` is meant to produce when SKA2 is installed. Until
then this number is reported as an unadjusted, confounded comparison.

---

*Generated during the M6/M7 milestone work. Reproduce with the commands above;
the engine version is in `run.json` of each output folder.*
