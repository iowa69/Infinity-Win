#!/usr/bin/env python3
"""Assemble the real benchmark datasets of PROJECT.md section 11.

For each dataset this script **copies** (never moves, never edits) the source
assemblies into ``bench-work/<dataset>/asm/`` and writes a ``dataset.json`` next
to them that records where every byte came from, its size and its SHA-256.

Three datasets are available on this machine:

``salmonella_trio``
    The three-isolate *Salmonella* Heidelberg outbreak set (isolate-A, isolate-B,
    isolate-C, ~5.13 Mb each). Published pairwise SNP distances exist for all
    three pairs, so this is the acceptance-check dataset
    (``scripts/expected/salmonella_trio.json``).

``ncbi_pd_cluster``
    The 32 genomes of NCBI Pathogen Detection cluster PDS000037076. A SKA2
    distance table and an IQ-TREE maximum-likelihood tree already exist for
    exactly these 32 genomes, so this dataset has an outside reference to be
    compared against.

``s_aureus_cc``
    *Staphylococcus aureus* reference genomes from the MORI collection. **These
    are different clonal complexes, i.e. far apart.** That makes them useful for
    checking that clustering separates distinct lineages and useless for checking
    SNP-level accuracy: at several percent divergence M-Tato's block shortcut
    stops helping (PROJECT.md section 15, "Divergent clusters"). The dataset.json
    and the report both say so. The collection is large, so a deterministic
    sample is taken rather than all of it.

Disk on the build machine is 99% full, which is why everything lands in
``bench-work/`` (gitignored) and why the *S. aureus* set is sampled.

Usage::

    scripts/real_datasets.py --list
    scripts/real_datasets.py --datasets salmonella_trio ncbi_pd_cluster
    scripts/real_datasets.py --datasets all
    scripts/real_datasets.py --datasets s_aureus_cc --max-per-cc 2 --max-extra-genomes 6
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

import common  # noqa: E402

SCHEMA_VERSION = 1

SALMONELLA_ROOT = Path("<salmonella-data>")
TRIO_ASM = SALMONELLA_ROOT / "asm"
NCBI_PD = SALMONELLA_ROOT / "ncbi_pd"
NCBI_PD_ASM = NCBI_PD / "asm"
SKA2_DISTANCE = NCBI_PD / "ska" / "cluster_reads_distance.tsv"
IQTREE_TREE = NCBI_PD / "tree" / "cluster.treefile"
IQTREE_BIONJ = NCBI_PD / "tree" / "cluster.bionj"
SALMONELLA_RESULTS = SALMONELLA_ROOT / "results.json"

MORI_CLOSED = Path("<mori-data>/refstrains/closed")
MORI_GENOMES = Path("<mori-data>/refstrains/genomes")

TRIO_NAMES = ("isolate-A", "isolate-B", "isolate-C")

FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")


# ---------------------------------------------------------------------------
# Source descriptions
# ---------------------------------------------------------------------------


@dataclass
class DatasetSpec:
    """What a dataset is and where it comes from."""

    name: str
    organism: str
    description: str
    licence: str
    sources: List[Path]
    notes: List[str] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)
    cluster_id: Optional[str] = None


SPECS: Dict[str, DatasetSpec] = {
    "salmonella_trio": DatasetSpec(
        name="salmonella_trio",
        organism="Salmonella enterica subsp. enterica serovar Heidelberg",
        description=(
            "Three-isolate outbreak set used as the small exact-accuracy "
            "acceptance check: published pairwise SNP distances exist for all "
            "three pairs."
        ),
        licence=(
            "Assemblies downloaded from NCBI (GenBank/Pathogen Detection); NCBI "
            "data are public domain. Copies are made for benchmarking only."
        ),
        sources=[TRIO_ASM],
        notes=[
            "Only the three isolates isolate-A, isolate-B and isolate-C are copied; "
            "the source folder holds nothing else that belongs to this set.",
            "Expected pairwise SNP distances are recorded in "
            "scripts/expected/salmonella_trio.json and are an acceptance check "
            "in REPORT.md, not a fit target.",
            "With three taxa the Robinson-Foulds distance is undefined (it needs "
            "at least four), so accuracy for this dataset is measured as SNP "
            "distance error only.",
        ],
        caveats=[
            "Three genomes is far below the scale M-Tato is designed for; this "
            "dataset tests exactness, not performance.",
        ],
    ),
    "ncbi_pd_cluster": DatasetSpec(
        name="ncbi_pd_cluster",
        organism="Salmonella enterica subsp. enterica serovar Heidelberg",
        description=(
            "All 32 genomes of NCBI Pathogen Detection cluster PDS000037076: a "
            "real, nearly identical outbreak-scale cluster."
        ),
        licence=(
            "Assemblies downloaded from NCBI Pathogen Detection; NCBI data are "
            "public domain. The SKA2 distance table and the IQ-TREE tree were "
            "produced locally from the corresponding public reads."
        ),
        sources=[NCBI_PD_ASM],
        cluster_id="PDS000037076",
        notes=[
            "The outside reference for this dataset already exists: a SKA2 "
            "distance table (496 rows for 32 genomes) and an IQ-TREE "
            "maximum-likelihood tree built from the SKA2 core alignment. Both "
            "are copied into truth/ so the dataset folder is self-contained.",
            "The IQ-TREE tree is a reference, not ground truth: it is a "
            "maximum-likelihood tree from a split-k-mer core alignment, and "
            "M-Tato's neighbour-joining tree is not expected to match it "
            "exactly.",
        ],
        caveats=[
            "SKA2 distances were computed from reads, M-Tato distances from "
            "assemblies, so small disagreements are expected and are reported "
            "rather than hidden.",
        ],
    ),
    "s_aureus_cc": DatasetSpec(
        name="s_aureus_cc",
        organism="Staphylococcus aureus",
        description=(
            "A sample of the MORI reference-genome collection: complete "
            "genomes labelled by clonal complex, plus a set of classic "
            "reference strains."
        ),
        licence=(
            "Genomes downloaded from NCBI (RefSeq, GCF_* accessions); NCBI data "
            "are public domain. Copies are made for benchmarking only."
        ),
        sources=[MORI_CLOSED, MORI_GENOMES],
        notes=[
            "The closed-genome folders are named CC<n>_<GCF accession>, so the "
            "clonal complex is known for those genomes and is recorded as a weak "
            "clustering expectation in dataset.json.",
            "The classic reference strains under genomes/ have no clonal-complex "
            "label in their file name, so they are recorded as 'unassigned'.",
            "The collection is much larger than a benchmark needs and the root "
            "partition is 99% full, so a deterministic sample is taken "
            "(--max-per-cc / --max-extra-genomes). Nothing is ever moved or "
            "modified in the source tree.",
        ],
        caveats=[
            "DIFFERENT CLONAL COMPLEXES MEAN THESE GENOMES ARE FAR APART. This "
            "dataset can show that clustering separates distinct lineages; it "
            "cannot validate SNP-level accuracy, and M-Tato's block shortcut is "
            "expected to stop helping at this divergence (PROJECT.md section 15).",
            "No ground-truth tree and no ground-truth SNP distances exist for "
            "this set, so `mtato bench` cannot score accuracy on it.",
        ],
    ),
}

DATASET_ORDER = ("salmonella_trio", "ncbi_pd_cluster", "s_aureus_cc")


# ---------------------------------------------------------------------------
# Copying
# ---------------------------------------------------------------------------


def fasta_files(folder: Path) -> List[Path]:
    out: List[Path] = []
    for suffix in FASTA_SUFFIXES:
        out.extend(sorted(folder.glob(f"*{suffix}")))
    # de-duplicate while keeping the sorted order
    seen = set()
    unique = []
    for p in out:
        if p not in seen:
            seen.add(p)
            unique.append(p)
    return unique


def genome_name(path: Path) -> str:
    name = path.name
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


@dataclass
class CopyRecord:
    name: str
    file: str
    source: str
    bytes: int
    sha256: str
    source_bytes: int
    source_sha256: str
    verified: bool
    clonal_complex: Optional[str] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        out = {
            "name": self.name,
            "file": self.file,
            "source": self.source,
            "bytes": self.bytes,
            "sha256": self.sha256,
            "source_bytes": self.source_bytes,
            "source_sha256": self.source_sha256,
            "copy_verified": self.verified,
        }
        if self.clonal_complex is not None:
            out["clonal_complex"] = self.clonal_complex
        out.update(self.extra)
        return out


def copy_one(src: Path, dst_dir: Path, name: str) -> Tuple[Path, bool, str]:
    """Copy ``src`` into ``dst_dir`` as ``<name><suffix>`` and verify the copy.

    Returns ``(destination, verified, problem)``. Verification compares size and
    SHA-256 of source and destination, which is the point of the exercise: the
    script proves it copied rather than asserting that it did.
    """
    suffix = "".join(src.suffixes)[-6:] if src.name.endswith(".gz") else src.suffix
    dst = dst_dir / f"{name}{suffix}"
    problem = ""
    if dst.exists() and dst.stat().st_size == src.stat().st_size:
        src_sum = common.sha256_of(src)
        dst_sum = common.sha256_of(dst)
        if src_sum == dst_sum:
            return dst, True, ""
        problem = f"{dst} existed but did not match {src}; re-copying"
        dst.unlink()
    dst_dir.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copy2(src, dst)
    except OSError as exc:
        return dst, False, f"copy failed {src} -> {dst}: {exc}"
    src_sum = common.sha256_of(src)
    dst_sum = common.sha256_of(dst)
    verified = src_sum == dst_sum and src.stat().st_size == dst.stat().st_size
    if not verified:
        problem = f"checksum mismatch after copy: {src} ({src_sum[:12]}) vs {dst} ({dst_sum[:12]})"
    return dst, verified, problem


CC_RE = re.compile(r"^(CC\d+)_")


def clonal_complex_of(name: str) -> Optional[str]:
    m = CC_RE.match(name)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------
# Per-dataset builders
# ---------------------------------------------------------------------------


def build_salmonella_trio(dest: Path, args: argparse.Namespace) -> Dict[str, Any]:
    spec = SPECS["salmonella_trio"]
    missing = [p for p in spec.sources if not p.is_dir()]
    if missing:
        raise common.PrereqError(
            f"source folder for salmonella_trio is missing: {missing[0]}\n"
            "       expected the three isolates "
            + ", ".join(TRIO_NAMES)
            + " as .fasta files there"
        )
    asm = common.ensure_dir(dest / "asm")
    records: List[CopyRecord] = []
    problems: List[str] = []
    for name in TRIO_NAMES:
        src = TRIO_ASM / f"{name}.fasta"
        if not src.exists():
            problems.append(f"source assembly missing: {src}")
            continue
        dst, verified, problem = copy_one(src, asm, name)
        if problem:
            problems.append(problem)
        records.append(
            CopyRecord(
                name=name,
                file=dst.name,
                source=str(src),
                bytes=dst.stat().st_size,
                sha256=common.sha256_of(dst),
                source_bytes=src.stat().st_size,
                source_sha256=common.sha256_of(src),
                verified=verified,
            )
        )
    if len(records) != 3:
        raise common.PrereqError(
            f"salmonella_trio: only {len(records)} of 3 isolates could be copied "
            f"from {TRIO_ASM}: {problems}"
        )

    truth: Dict[str, Any] = {
        "kind": "published_expected_values",
        "dir": None,
        "tree": None,
        "distances": None,
        "expected_values": str(
            common.repo_root() / "scripts" / "expected" / "salmonella_trio.json"
        ),
        "notes": [
            "Expected pairwise SNP distances come from two independent published "
            "analyses of these three isolates; see the expected-values file for "
            "per-pair provenance.",
            "Fewer than four taxa: Robinson-Foulds accuracy is undefined here, so "
            "`mtato bench` should compare pairwise distances only.",
        ],
    }
    # A self-contained 3-row SKA2 extract, when the full table is available.
    truth_dir = common.ensure_dir(dest / "truth")
    if SKA2_DISTANCE.exists():
        try:
            ska = common.parse_ska2_distance_tsv(SKA2_DISTANCE)
            wanted = {common.pair_key(a, b) for a in TRIO_NAMES for b in TRIO_NAMES if a != b}
            rows = [p for p in ska.pairs if p.key() in wanted]
            if rows:
                out_tsv = truth_dir / "ska2_trio_distance.tsv"
                common.write_tsv(
                    out_tsv,
                    ["Sample1", "Sample2", "Distance", "Mismatches (proportion)",
                     "Match count", "Mismatch count"],
                    [
                        [r.sample1, r.sample2, r.distance, r.mismatch_proportion,
                         r.match_count, r.mismatch_count]
                        for r in rows
                    ],
                )
                truth["distances"] = {
                    "format": "ska2",
                    "path": str(out_tsv),
                    "rows": len(rows),
                    "extracted_from": str(SKA2_DISTANCE),
                }
            else:
                problems.append(
                    f"no trio rows found in {SKA2_DISTANCE}; the SKA2 extract was not written"
                )
        except common.PrereqError as exc:
            problems.append(f"could not read {SKA2_DISTANCE}: {exc}")
    else:
        problems.append(f"SKA2 table not present at {SKA2_DISTANCE}; no SKA2 extract written")

    expected_path = Path(truth["expected_values"])
    if not expected_path.exists():
        problems.append(f"expected-values file missing: {expected_path}")

    return _dataset_doc(
        spec=spec,
        dest=dest,
        records=records,
        truth=truth,
        problems=problems,
        args=args,
    )


def build_ncbi_pd_cluster(dest: Path, args: argparse.Namespace) -> Dict[str, Any]:
    spec = SPECS["ncbi_pd_cluster"]
    if not NCBI_PD_ASM.is_dir():
        raise common.PrereqError(
            f"source folder for ncbi_pd_cluster is missing: {NCBI_PD_ASM}\n"
            "       expected the 32 assemblies of cluster PDS000037076 there"
        )
    sources = fasta_files(NCBI_PD_ASM)
    if not sources:
        raise common.PrereqError(f"no FASTA assemblies found in {NCBI_PD_ASM}")
    if args.limit and len(sources) > args.limit:
        sources = sources[: args.limit]
    asm = common.ensure_dir(dest / "asm")
    records: List[CopyRecord] = []
    problems: List[str] = []
    for src in sources:
        name = genome_name(src)
        dst, verified, problem = copy_one(src, asm, name)
        if problem:
            problems.append(problem)
        records.append(
            CopyRecord(
                name=name,
                file=dst.name,
                source=str(src),
                bytes=dst.stat().st_size,
                sha256=common.sha256_of(dst),
                source_bytes=src.stat().st_size,
                source_sha256=common.sha256_of(src),
                verified=verified,
            )
        )

    truth_dir = common.ensure_dir(dest / "truth")
    truth: Dict[str, Any] = {
        "kind": "external_reference",
        "dir": str(truth_dir),
        "tree": None,
        "distances": None,
        "notes": [
            "Both reference files were produced outside M-Tato and are copied "
            "here unchanged so that this dataset folder is self-contained.",
            "`mtato bench --truth` only accepts a `mtato simulate` truth folder "
            "(tree_*.nwk, pairwise_snps.tsv, clusters.tsv), so it cannot score "
            "this real dataset. The comparisons that do apply are made by the "
            "scripts: `baselines.py compare` joins these SKA2 distances against "
            "the run's cluster_*.dist.tsv and computes RF against "
            "tree/cluster.treefile, and `accuracy_over_time.py` builds the "
            "accuracy curve from events.jsonl against the same reference tree.",
        ],
    }
    if IQTREE_TREE.exists():
        dst = truth_dir / "cluster.treefile"
        shutil.copy2(IQTREE_TREE, dst)
        truth["tree"] = str(dst)
        truth["tree_provenance"] = {
            "source": str(IQTREE_TREE),
            "method": "IQ-TREE maximum likelihood on the SKA2 core alignment",
            "sha256": common.sha256_of(dst),
        }
        # The reference tree may cover more taxa than the assemblies copied here
        # (outgroups, reference genomes). Comparisons are restricted to the shared
        # leaves; record that instead of letting it show up as a mysterious gap.
        try:
            ref_tree = common.read_newick_file(dst)
            ref_leaves = set(ref_tree.leaf_names())
            copied = {r.name for r in records}
            extra = sorted(ref_leaves - copied)
            missing = sorted(copied - ref_leaves)
            truth["tree_leaf_check"] = {
                "reference_tree_leaves": len(ref_leaves),
                "assemblies_copied": len(copied),
                "in_both": len(ref_leaves & copied),
                "only_in_reference_tree": extra,
                "only_in_assemblies": missing,
            }
            if extra:
                truth.setdefault("notes", []).append(
                    f"The reference tree has {len(ref_leaves)} leaves but only "
                    f"{len(copied)} assemblies are in this dataset: "
                    f"{len(extra)} extra taxon/taxa ({', '.join(extra[:4])}"
                    f"{'...' if len(extra) > 4 else ''}) are reference or outgroup "
                    "genomes that were not part of the cluster download. Tree "
                    "comparisons are restricted to the shared leaves."
                )
            if missing:
                problems.append(
                    f"{len(missing)} copied genome(s) do not appear in the reference "
                    f"tree: {missing[:5]}"
                )
        except common.PrereqError as exc:
            problems.append(f"could not parse the copied reference tree: {exc}")
    else:
        problems.append(f"reference tree not present at {IQTREE_TREE}")
    if IQTREE_BIONJ.exists():
        dst = truth_dir / "cluster.bionj"
        shutil.copy2(IQTREE_BIONJ, dst)
        truth["tree_bionj"] = str(dst)
    if SKA2_DISTANCE.exists():
        dst = truth_dir / "ska2_cluster_reads_distance.tsv"
        shutil.copy2(SKA2_DISTANCE, dst)
        ska_samples: List[str] = []
        try:
            ska = common.parse_ska2_distance_tsv(dst)
            stats = ska.stats()
            ska_samples = ska.samples()
        except common.PrereqError as exc:
            stats = {"error": str(exc)}
        truth["distances"] = {
            "format": "ska2",
            "path": str(dst),
            "source": str(SKA2_DISTANCE),
            "sha256": common.sha256_of(dst),
            "stats": stats,
        }
    else:
        ska_samples = []
        problems.append(f"SKA2 distance table not present at {SKA2_DISTANCE}")

    names = {r.name for r in records}
    if truth.get("distances"):
        ref_samples = set(ska_samples)
        if ref_samples:
            missing_ref = sorted(names - ref_samples)
            extra_ref = sorted(ref_samples - names)
            if missing_ref:
                problems.append(
                    f"{len(missing_ref)} copied genome(s) are absent from the SKA2 "
                    f"reference table: {missing_ref[:5]}"
                )
            if extra_ref:
                problems.append(
                    f"{len(extra_ref)} sample(s) in the SKA2 reference table were "
                    f"not copied: {extra_ref[:5]}"
                )
            truth["sample_name_overlap"] = {
                "in_both": len(names & ref_samples),
                "only_in_asm": missing_ref,
                "only_in_reference": extra_ref,
            }

    return _dataset_doc(
        spec=spec, dest=dest, records=records, truth=truth, problems=problems, args=args
    )


def build_s_aureus_cc(dest: Path, args: argparse.Namespace) -> Dict[str, Any]:
    spec = SPECS["s_aureus_cc"]
    if not MORI_CLOSED.is_dir() and not MORI_GENOMES.is_dir():
        raise common.PrereqError(
            f"neither S. aureus source folder exists: {MORI_CLOSED} or {MORI_GENOMES}"
        )
    problems: List[str] = []
    if not MORI_CLOSED.is_dir():
        problems.append(f"closed-genome folder missing: {MORI_CLOSED}")
    if not MORI_GENOMES.is_dir():
        problems.append(f"genomes folder missing: {MORI_GENOMES}")

    closed: List[Tuple[Path, str, Optional[str]]] = []
    by_cc: Dict[str, List[Path]] = {}
    if MORI_CLOSED.is_dir():
        for folder in sorted(p for p in MORI_CLOSED.iterdir() if p.is_dir()):
            found = sorted(f for f in folder.glob("*.fna")) + sorted(
                f for f in folder.glob("*.fasta")
            )
            if not found:
                problems.append(f"no .fna/.fasta in {folder}")
                continue
            if len(found) > 1:
                problems.append(f"{folder} holds {len(found)} assemblies; using {found[0].name}")
            cc = clonal_complex_of(folder.name)
            by_cc.setdefault(cc or "unknown", []).append(found[0])
    for cc in sorted(by_cc):
        limit = None if args.all else args.max_per_cc
        chosen = by_cc[cc] if limit is None else by_cc[cc][:limit]
        if limit is not None and len(by_cc[cc]) > limit:
            problems.append(
                f"sampled {limit} of {len(by_cc[cc])} genomes in {cc} "
                f"(use --all to take every one)"
            )
        for path in chosen:
            label = path.stem
            closed.append((path, label, cc if cc != "unknown" else None))

    extra: List[Tuple[Path, str, Optional[str]]] = []
    if MORI_GENOMES.is_dir():
        all_extra = sorted(MORI_GENOMES.glob("*.fna")) + sorted(MORI_GENOMES.glob("*.fasta"))
        limit = None if args.all else args.max_extra_genomes
        chosen = all_extra if limit is None else all_extra[:limit]
        if limit is not None and len(all_extra) > limit:
            problems.append(
                f"sampled {limit} of {len(all_extra)} reference genomes from "
                f"{MORI_GENOMES} (use --all to take every one)"
            )
        for path in chosen:
            extra.append((path, genome_name(path), None))

    selected = closed + extra
    if not selected:
        raise common.PrereqError(
            "no S. aureus assemblies could be selected from "
            f"{MORI_CLOSED} or {MORI_GENOMES}"
        )

    asm = common.ensure_dir(dest / "asm")
    records: List[CopyRecord] = []
    used: Dict[str, int] = {}
    for src, label, cc in selected:
        name = label.replace(" ", "_")
        if name in used:
            used[name] += 1
            name = f"{name}_{used[name]}"
        else:
            used[name] = 1
        dst, verified, problem = copy_one(src, asm, name)
        if problem:
            problems.append(problem)
        records.append(
            CopyRecord(
                name=name,
                file=dst.name,
                source=str(src),
                bytes=dst.stat().st_size,
                sha256=common.sha256_of(dst),
                source_bytes=src.stat().st_size,
                source_sha256=common.sha256_of(src),
                verified=verified,
                clonal_complex=cc,
            )
        )

    expected_clusters: Dict[str, List[str]] = {}
    for rec in records:
        if rec.clonal_complex:
            expected_clusters.setdefault(rec.clonal_complex, []).append(rec.name)
    truth: Dict[str, Any] = {
        "kind": "weak_clustering_expectation",
        "dir": None,
        "tree": None,
        "distances": None,
        "expected_clusters": {
            "note": (
                "Only the closed genomes carry a clonal-complex label in their "
                "folder name. Genomes labelled 'unassigned' have no expected "
                "cluster. This is a WEAK expectation: distinct clonal complexes "
                "are far enough apart that any reasonable method separates them, "
                "so passing it says almost nothing about M-Tato's accuracy."
            ),
            "clusters": expected_clusters,
            "unassigned": [r.name for r in records if not r.clonal_complex],
        },
        "notes": [
            "No reference tree and no reference SNP distances exist for this set, "
            "so `mtato bench` cannot score accuracy on it. It is included to test "
            "that clustering separates lineages that are far apart, and to show "
            "where the block shortcut stops paying.",
        ],
    }

    return _dataset_doc(
        spec=spec, dest=dest, records=records, truth=truth, problems=problems, args=args
    )


BUILDERS = {
    "salmonella_trio": build_salmonella_trio,
    "ncbi_pd_cluster": build_ncbi_pd_cluster,
    "s_aureus_cc": build_s_aureus_cc,
}


# ---------------------------------------------------------------------------
# dataset.json
# ---------------------------------------------------------------------------


def _dataset_doc(
    *,
    spec: DatasetSpec,
    dest: Path,
    records: List[CopyRecord],
    truth: Dict[str, Any],
    problems: List[str],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    total_bytes = sum(r.bytes for r in records)
    unverified = [r.name for r in records if not r.verified]
    if unverified:
        problems.append(
            f"{len(unverified)} file(s) failed checksum verification after copying: "
            f"{unverified[:5]}"
        )
    return {
        "schema_version": SCHEMA_VERSION,
        "name": spec.name,
        "kind": "real",
        "organism": spec.organism,
        "description": spec.description,
        "cluster_id": spec.cluster_id,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "created_by": "scripts/real_datasets.py",
        "host": {
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "request": {
            "datasets": list(args.datasets),
            "all": bool(args.all),
            "max_per_cc": args.max_per_cc,
            "max_extra_genomes": args.max_extra_genomes,
            "limit": args.limit,
            "checksums": not args.no_checksums,
        },
        "asm_dir": str(dest / "asm"),
        "genome_count": len(records),
        "total_bytes": total_bytes,
        "copy_policy": (
            "sources were COPIED, never moved or modified; every copy was "
            "verified by size and SHA-256 against its source"
        ),
        "copy_verified": not unverified,
        "genomes": [r.as_dict() for r in records],
        "source": {
            "description": spec.description,
            "paths": [str(p) for p in spec.sources],
            "exists": {str(p): p.exists() for p in spec.sources},
        },
        "licence": spec.licence,
        "truth": truth,
        "notes": list(spec.notes),
        "caveats": list(spec.caveats),
        "problems": problems,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def describe_sources() -> List[Tuple[str, str, bool]]:
    rows = []
    for name in DATASET_ORDER:
        spec = SPECS[name]
        ok = all(p.exists() for p in spec.sources)
        rows.append((name, ", ".join(str(p) for p in spec.sources), ok))
    return rows


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="real_datasets.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--datasets",
        nargs="+",
        default=list(DATASET_ORDER),
        choices=[*DATASET_ORDER, "all"],
        help="which datasets to assemble (default: all three)",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help=f"benchmark work directory (default: {common.default_work_dir()})",
    )
    p.add_argument(
        "--all",
        action="store_true",
        help="for s_aureus_cc, take every genome instead of a sample",
    )
    p.add_argument(
        "--max-per-cc",
        type=int,
        default=2,
        help="s_aureus_cc: how many closed genomes per clonal complex (default 2)",
    )
    p.add_argument(
        "--max-extra-genomes",
        type=int,
        default=6,
        help="s_aureus_cc: how many classic reference genomes (default 6)",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="ncbi_pd_cluster: cap the number of genomes copied (default: all 32)",
    )
    p.add_argument(
        "--no-checksums",
        action="store_true",
        help="skip SHA-256 verification of the copies (not recommended)",
    )
    p.add_argument(
        "--clean",
        action="store_true",
        help="delete an existing bench-work/<dataset> before rebuilding it",
    )
    p.add_argument("--list", action="store_true", help="list datasets and sources, then exit")
    p.add_argument("--json", type=Path, default=None, help="also write a summary JSON here")
    p.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.list:
        print("dataset            sources present  paths")
        for name, paths, ok in describe_sources():
            print(f"{name:<18} {'yes' if ok else 'NO ':<16} {paths}")
        return 0

    work = Path(args.work_dir) if args.work_dir else common.default_work_dir()
    names = list(DATASET_ORDER) if "all" in args.datasets else list(args.datasets)
    common.ensure_dir(work)

    results: List[Dict[str, Any]] = []
    failures: List[str] = []
    print(f"work directory: {work}")
    print(f"requested:      {', '.join(names)}\n")

    for name in names:
        builder = BUILDERS.get(name)
        if builder is None:
            failures.append(f"{name}: no builder (internal error)")
            continue
        dest = work / name
        if args.clean and dest.exists():
            shutil.rmtree(dest)
            print(f"[{name}] cleaned {dest}")
        started = time.monotonic()
        try:
            doc = builder(dest, args)
        except common.PrereqError as exc:
            failures.append(f"{name}: {exc}")
            print(f"[{name}] FAILED: {exc}\n")
            results.append({"name": name, "status": "failed", "error": str(exc)})
            continue
        except OSError as exc:
            failures.append(f"{name}: {exc}")
            print(f"[{name}] FAILED: {exc}\n")
            results.append({"name": name, "status": "failed", "error": str(exc)})
            continue
        path = common.dump_json(dest / "dataset.json", doc)
        elapsed = time.monotonic() - started
        print(
            f"[{name}] {doc['genome_count']} genomes, "
            f"{common.human_bytes(doc['total_bytes'])}, "
            f"copy_verified={doc['copy_verified']}, {elapsed:.1f}s"
        )
        print(f"           asm:         {doc['asm_dir']}")
        print(f"           dataset.json: {path}")
        truth = doc.get("truth") or {}
        if truth.get("tree"):
            print(f"           truth tree:  {truth['tree']}")
        if truth.get("distances"):
            d = truth["distances"]
            print(f"           truth dist:  {d.get('path')} ({d.get('format')})")
        if truth.get("expected_values"):
            print(f"           expected:    {truth['expected_values']}")
        for problem in doc.get("problems", []):
            print(f"           NOTE: {problem}")
        for caveat in doc.get("caveats", []):
            print(f"           CAVEAT: {caveat}")
        print()
        results.append(
            {
                "name": name,
                "status": "ok",
                "dataset_json": str(path),
                "genome_count": doc["genome_count"],
                "total_bytes": doc["total_bytes"],
                "copy_verified": doc["copy_verified"],
                "problems": doc.get("problems", []),
            }
        )

    if args.json:
        common.dump_json(
            args.json,
            {
                "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "work_dir": str(work),
                "results": results,
                "failures": failures,
            },
        )
        print(f"summary written to {args.json}")

    if failures:
        print(f"{len(failures)} dataset(s) could not be assembled:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    common.guard(main)()
