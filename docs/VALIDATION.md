# Validation of Infinity-Win 0.2 on real data

All numbers below were produced on 2026-11-18 with Infinity-Win 0.2.0 at default
settings (`infinity run <folder> -o <out>`), on a 24-core Linux workstation. The
reference for real data is **SKA2 0.5.1 run on the sequencing reads**, which is
the closest thing to truth available: it does not depend on how the genomes were
assembled. Where Infinity-Win and SKA2 disagree on a close pair, every disputed
SNP was checked directly in the raw reads (count of reads carrying each allele).

## 1. *Salmonella* Heidelberg outbreak trio (owner's acceptance test)

Three assemblies of 57–66 contigs. Published SNP distances (SKA2 on reads): 2, 2, 0.

| pair | published | 0.1.0 | **0.2.0** | 0.2.0 raw (no rules) |
|---|---|---|---|---|
| isolate-A – isolate-B | 2 | 17 | **2** | 20 |
| isolate-A – isolate-C | 2 | 2 | **2** | 4 |
| isolate-B – isolate-C | 0 | 4 | **0** | 13 |

The raw differences sit on a 1.4 kb contig at 3.4x the chromosomal coverage (a
collapsed high-copy element) and at contig ends; the SNP rules remove both.

## 2. *Salmonella* NCBI Pathogen Detection cluster (32 genomes)

SPAdes and SKESA assemblies with 30–96 contigs each; 496 pairs; SKA2-on-reads
distances 0–46 SNPs.

| | 0.1.0 | **0.2.0** | SKA2 on the same *assemblies* |
|---|---|---|---|
| pairs exact | 12 | **109** | 63 vs reads |
| within 2 SNPs | 35 | **381** | 145 |
| within 5 SNPs | ~120 | **494** | – |
| worst error | 130 | **6** | 521 |
| mean bias | +13.7 | **−0.2** | +68.9 |
| correlation | 0.41 | **0.977** | 0.07 |
| run time | 1.7 s | **1.8 s** | – |

SKA2 itself, run on these assemblies instead of the reads, is far noisier than
Infinity-Win: the assemblies carry errors that a whole-assembly comparison sees.
MUMmer `dnadiff` against the pivot agrees with Infinity-Win's raw counts to within
2 SNPs for 17 of 31 genomes and is never lower, so the engine reports the real
differences between the assemblies; the rules decide which of them count.

## 3. *S. aureus* bloodstream cohort (137 genomes, MORI)

137 real assemblies (about 60 contigs each), 24 clusters at the default 99.5%
sketch identity. First tree after 0.6 s, every cluster exact after 27 s, peak
memory 1.7 GB. SKA2 was run on the reads of every multi-genome cluster (126
isolates, 16 GB of reads, 20 min).

**The four transmission pairs** (isolate identifiers omitted):

| pair | Infinity-Win | SKA2 on reads | checked in the raw reads |
|---|---|---|---|
| A | 8 | 6 | **8 real** (SKA2 misses two adjacent SNPs: split k-mers cannot see SNPs that close) |
| B | 8 | 8 | not checked (methods agree) |
| C | 14 | 8 | **8 real**, 4 collapsed-repeat artefacts, 2 unclear |
| D | 15 | 15 | **14 real**, 1 collapsed-repeat artefact |

Before the four-gamete test, pair A was 11: three SNPs inside one kilobase whose
alleles were present in *both* isolates' reads. The test recognises such sites
from the assemblies alone, because neighbouring sites are carried by unrelated
subsets of genomes, which a clonal tree cannot produce.

**What remains:** isolated SNPs inside collapsed multi-copy regions (pairs C and
D). Each assembly holds a single copy there, so no rule based on assemblies can
tell; read depth can (5–7x the genome average at those sites). This is the
strongest argument for adding reads (FASTQ) as an input.

**Distant pairs.** Inside the 24 clusters most pairs are different strains of one
lineage, 160–2,300 SNPs apart by SKA2. There the rules remove the SNP-dense
windows that recombination leaves, so Infinity-Win reports fewer SNPs than SKA2
(correlation 0.77). For questions at that distance, use the `raw_snps` columns or
run with `--no-snp-filter`.

## 4. Real genomes, fragmented, with known truth

Three closed *S. aureus* genomes from different lineages (2.8–2.9 Mb, with
plasmids) as references for the simulator: 2 clusters x 25 descendants each,
3 SNPs, 0.5 indels and 0.3 segment gains or losses per branch, every genome cut
into 60 contigs in random order and orientation. 600 pairs per reference.

| | before the contig-end fix | **0.2.0** |
|---|---|---|
| uncovered bases per contig | ~800 | **~70** |
| pairs over-reported | 0 | **0** |
| true SNPs counted (strict core) | 82% | **93%** |
| true SNPs found before the rules (raw) | 95% | **99%** |

The 7% not counted lie outside the strict core by design: lost segments, contig
ends, and bases next to contig junctions.

## How to reproduce

```sh
infinity run <assembly folder> -o out --threads 16
# SKA2 on reads, per cluster (sample<TAB>R1<TAB>R2 per line):
ska build -f cluster.tsv -k 31 --threads 12 -o cluster && ska distance cluster.skf > ska.tsv
# known-truth fragmented runs:
infinity simulate -o sim --reference closed_genome.fna -m 2 -n 25 --divergence 0.003 \
  --mean-snps 3 --mean-indels 0.5 --mean-segments 0.3 -c 60 --random-orientation
infinity run sim/genomes -o sim_out && infinity bench --output sim_out --truth sim/truth
```
