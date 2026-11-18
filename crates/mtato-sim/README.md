# mtato-sim

Deterministic genome simulator with known truth, for M-Tato (section 10.1 of
`PROJECT.md`). It grows a pedigree of nearly identical bacterial genomes from a
reference, records every mutation it makes, and writes both the genomes and the
record out, so that the engine's accuracy can be *measured* instead of argued
about. Tests T1-T5 of section 10.3 and the accuracy columns of section 11 all
read the files this crate produces.

Library only: the CLI exposes it as `mtato simulate`.

```rust
use mtato_sim::{simulate, RefSpec, SimConfig};

let cfg = SimConfig {
    seed: 20260118,
    reference: RefSpec::Random { length: 2_000_000, gc: 0.5 },
    clusters: 1,
    leaves_per_cluster: 50,
    mean_snps_per_branch: 3.0,
    ..SimConfig::default()
};
let out = simulate(&cfg, std::path::Path::new("/tmp/t1")).unwrap();
assert_eq!(out.genomes.len(), 50);
let names = out.truth.genome_names();
let snps = out.truth.pairwise_snps(&names[0], &names[1]);
```

## What it builds

1. **A reference.** Either a FASTA the user supplies - its contigs are
   concatenated into one linear sequence - or an i.i.d. random sequence of a
   chosen length and GC fraction.
2. **`m` cluster roots.** Each is the reference with `divergence` (default 2%) of
   its positions substituted.
3. **A pedigree per cluster.** A random binary tree with `n` leaves. Every branch
   draws its SNP count from a Poisson distribution with mean
   `mean_snps_per_branch`, and that one number is both how many substitutions the
   branch applies and what the truth Newick says the branch length is, so the two
   cannot drift apart.
4. **Optional branch events.** Small indels of 1-10 bases; gains and losses of
   500-5000 base segments.
5. **Optional fragmentation.** Cut each genome into `contigs` pieces, shuffle
   them, and - with `random_orientation` - reverse-complement each independently.

## The invariant that makes the truth exact

Inside a cluster, **no two events ever touch the same reference position**: not
two SNPs, not a SNP and a deletion, not two insertion anchors, and not any event
and any event of any *other* branch of the same cluster. A branch SNP also never
lands on a position an ancestor deleted, or on one of the substitutions that made
the cluster root.

Two consequences, both of which the outline depends on:

- The true SNP distance between two genomes is exactly the sum of the branch
  lengths on the path between them in the truth tree. Nothing cancels out and
  nothing is overwritten, which is what lets T1 demand an exact match.
- A leaf's mutation log is a set of edits at pairwise distinct reference
  coordinates, so two leaves differ at exactly the symmetric difference of their
  private SNP positions - and the log is a replay script that reproduces the
  genome from its cluster root.

Segment losses are the one place where truth and observation legitimately
diverge: a lost stretch is not a SNP difference, it is missing coverage. The
substitutions the *cluster root* carried inside it disappear with it, but those
are shared by every member of the cluster, so they cancel out of every
within-cluster distance. This is why T2 and T4 allow a pair to be off by up to
two rather than demanding exactness.

## Determinism

Same seed, same configuration, same reference file: byte-identical FASTA and
byte-identical truth files. Nothing is read from the clock, the generator is this
crate's own (xoshiro256\*\* seeded from SplitMix64) rather than an external one
whose version could silently change a benchmark, no hash map is ever iterated on
the way to a file, and every derived stream is keyed by what it is for. Two
useful consequences:

- Changing `contigs` or `random_orientation` changes how genomes are cut, not one
  base of what they contain.
- Changing one cluster's seed does not move another cluster.

The one thing the seed does not fix is the reference file: a different FASTA
gives a different dataset, as it must. Float draws (`poisson`, GC content) go
through the platform libm, so bit-identical output across *different* platforms
is not promised; the truth files are the record, and they travel with the data.

## Output layout

Written under the output directory given to `simulate`:

```text
reference.fasta                  the sequence the cluster roots were derived from
genomes/<genome>.fasta           one file per leaf genome
roots/<root>.fasta               one file per cluster root (an ancestor, not a leaf)
truth/manifest.json              format version, config and seed, reference, genome list
truth/clusters.tsv               genome, cluster_id, file, length, contigs
truth/pairwise_snps.tsv          genome_a, genome_b, snps
truth/tree_c<NN>.nwk             true tree of cluster NN, branch lengths in SNPs
truth/mutations_<name>.tsv       kind, position, ref, alt - one file per leaf and per root
```

Point the engine at `genomes/`, not at the output directory: `reference.fasta`
and `roots/` are truth artefacts, not simulated genomes.

Every TSV has a header row, is UTF-8, uses `\n` line endings on every platform,
and ends with a single newline. `manifest.json` is pretty-printed JSON.

Names are `<name_prefix>_c<cluster>_s<leaf>` with both numbers zero-padded to a
fixed width, so name order and cluster/leaf order agree and a lexicographic sort
is also a numeric one. Cluster roots are `root_c<NN>`. A FASTA record is named
after its genome when the file holds one contig, and `<genome>_ctg<i>` (1-based,
in output order) when it holds more.

### `truth/clusters.tsv`

| Column | Meaning |
|---|---|
| `genome` | Leaf genome name, also the FASTA file stem |
| `cluster_id` | 0-based cluster the genome belongs to |
| `file` | Its FASTA path, relative to the output directory |
| `length` | Total bases over all contigs |
| `contigs` | Number of records in that file |

The first two columns are the assignment section 10.1 asks for; the rest lets a
script sanity-check a run without opening the FASTA files. The same information
is in `manifest.json`, and `Truth::read_from` refuses a directory where the two
disagree.

### `truth/pairwise_snps.tsv`

| Column | Meaning |
|---|---|
| `genome_a`, `genome_b` | Two leaves of the same cluster, `genome_a < genome_b` |
| `snps` | True SNP count between them: the path length in the truth tree |

One row per unordered pair *inside* a cluster, rows sorted by the pair.
Cross-cluster pairs are deliberately absent: the engine never measures them
(section 6.4), so a truth value would only invite a comparison that cannot be
made. `Truth::pairwise_snps` returns `None` for them and `Some(0)` for a genome
against itself.

### `truth/tree_c<NN>.nwk`

Newick, written with `mtato_core::tree::to_newick`, rooted at the cluster root,
with integer branch lengths in SNPs. `NN` is the cluster id padded to at least
two digits. `mtato_sim::pairwise_from_newick` recomputes the whole pairwise
matrix from this file alone, which is the independent check truth test 3 is made
of.

### `truth/mutations_<name>.tsv`

| Column | Meaning |
|---|---|
| `kind` | `SNP`, `INS`, `DEL`, `GAIN` or `LOSS` |
| `position` | 1-based reference coordinate, see below |
| `ref` | Bases before the event, or `-` when there were none |
| `alt` | Bases after the event, or `-` when there are none |

For a leaf this is the log of the path from its cluster root to it; for a
`root_c<NN>` file it is the path from the reference to that root, which is all
substitutions. Rows are sorted by position.

Positions are coordinates on the reference, the one sequence every genome of a
cluster descends from. Section 5's rule applies: 0-based in the program
(`Mutation::pos`), 1-based in the files. For `SNP`, `DEL` and `LOSS` the position
is the first affected base. For `INS` and `GAIN` it is the base immediately
*after* the inserted sequence, which is what anchors an insertion to a coordinate
it does not itself occupy.

`ref` for a `LOSS` is the sequence as it stood in that lineage, so it carries any
substitution the cluster root introduced inside it - that is the string a search
of the root's FASTA will find. The length change an event caused is
`len(alt) - len(ref)`, so a genome's total length is its root's length plus the
sum of the deltas in its log.

### `truth/manifest.json`

Carries `format_version`, the crate version, the whole `SimConfig` with its seed,
the reference's provenance and length, and one entry per leaf genome and per
cluster root with its file path, length and contig count. Everything except the
mutation rows, the Newick strings and the pairwise table is in this file, so a
dataset can be regenerated exactly from it.

## Consuming the truth

From Rust, `simulate` returns the truth in memory as well as on disk, so an
integration test does not have to read back what it just wrote:

```rust
let out = simulate(&SimConfig::default(), std::path::Path::new("/tmp/x")).unwrap();
let names = out.truth.genome_names();
let d = out.truth.pairwise_snps(&names[0], &names[1]);
let log = out.truth.mutations_of(&names[0]).unwrap();
```

`Truth::write_to(outdir)` and `Truth::read_from(outdir)` are exact inverses -
`outdir` is the simulation output directory in both cases, and the files live in
`outdir/truth`. Useful accessors: `cluster_of`, `clusters`, `cluster_ids`,
`members`, `genome`, `root`, `mutations_of`, `tree`, `pairwise_snps`,
`genome_names`.

From Python, the TSV files are plain `csv.reader(..., delimiter='\t')` material
with a header row, `manifest.json` is `json.load`, and the Newick files parse
with any reader (ete3, Bio.Phylo); branch lengths are integers.

## Configuration

`SimConfig::default()` is the section 10.1 default: a 2 Mb random reference at
50% GC, one cluster, ten leaves, 2% divergence, three SNPs per branch, no indels,
no segment events, one forward contig per genome. Every field is
serde-serialisable, so the CLI can drop the whole configuration into `run.json`.

The section 10.3 datasets are:

| Test | Changes from the default |
|---|---|
| T1 | `leaves_per_cluster: 50` |
| T2 | T1 plus `mean_indels_per_branch: 0.5` |
| T3 | T2 plus `contigs: 100`, `random_orientation: true` |
| T4 | T2 plus `segment_events` with non-zero gain and loss means |
| T5 | `clusters: 5`, `leaves_per_cluster: 20`, `divergence: 0.02` |

`SimConfig::validate` checks the numbers without touching the filesystem;
problems that depend on the reference (it does not exist, it is shorter than the
contig count) are reported by `simulate`.

## Tests

```
cargo test -p mtato-sim
```

Unit tests cover the generator, the occupancy bitsets, the spine, the pedigree
builder and the truth file round-trip. `tests/simulation.rs` covers the
properties that make the datasets usable as ground truth:

1. Determinism - the same seed gives byte-identical files, a different seed does
   not, and rerunning into the same directory changes nothing.
2. Mutation exclusivity - no two distinct events in a cluster touch one
   position, and the symmetric difference of two leaves' SNP sets is exactly
   their truth distance.
3. Truth self-consistency - the published matrix equals the path sums read back
   out of the Newick with `mtato_core::tree::parse_newick`.
4. Reference round-trip - a supplied FASTA gives roots at the requested
   divergence, no substitution on an N, and unchanged genome lengths.
5. Indels - every base of length change is accounted for by a recorded indel,
   and replaying the log rebuilds the genome.
6. Segment gain and loss - gained segments are in the leaf and not in the root,
   lost ones the reverse, at the recorded coordinate and length.
7. Fragmentation - exactly `c` contigs, total length preserved, and the pieces
   tile the unfragmented genome up to orientation.
8. Pinned small-scale numbers - 3 leaves at 2 SNPs per branch, checked by hand
   against the topology, the logs and the sequences.
9. Truth round-trip through the files.
10. Property 4 of section 10.2 from the simulator side - a 2 Mb simulated genome
    gives `mtato_core::genome::scan_data` an anchor count within 10% of
    `length / densities[2]`, and its reverse complement gives the same anchor
    hashes.
