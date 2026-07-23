"""Native, deterministic outbreak genome and read simulation."""

from __future__ import annotations

import hashlib
import json
import os
import random
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from genomepuzzle.contract import json_dump
from genomepuzzle.outbreak import build_outbreak_release, read_metadata
from genomepuzzle.read_generation import _simulate_short_reads
from genomepuzzle.release import ReleaseSpec, resolve_release_samples
from genomepuzzle.typing import read_fasta


BASES = "ACGT"


def _stable_seed(*parts: object) -> int:
    payload = "\0".join(str(part) for part in parts)
    value = int.from_bytes(hashlib.sha256(payload.encode("utf-8")).digest()[:4], "big")
    return value & 0x7FFFFFFF or 1


def _reference_sequence(path: Path) -> str:
    records = read_fasta(path)
    sequence = "".join(item[1] for item in records)
    if len(sequence) < 1000:
        raise ValueError("outbreak base genome must contain at least 1000 bases")
    return sequence


def _mutations(sequence: str, count: int, seed: int, excluded: set[int]) -> dict[int, str]:
    if count < 0:
        raise ValueError("mutation counts must be non-negative")
    available = [index for index, base in enumerate(sequence) if base in BASES and index not in excluded]
    if count > len(available):
        raise ValueError("requested more mutations than available reference positions")
    rng = random.Random(seed)
    selected = rng.sample(available, count)
    return {
        index: rng.choice(BASES.replace(sequence[index], ""))
        for index in selected
    }


def _apply(sequence: str, mutations: dict[int, str]) -> str:
    result = list(sequence)
    for position, base in mutations.items():
        result[position] = base
    return "".join(result)


def generate_outbreak_release(
    spec: ReleaseSpec,
    base_genome: str | os.PathLike[str],
    metadata_csv: str | os.PathLike[str],
    release_dir: str | os.PathLike[str],
    *,
    id_salt: str | None = None,
    shared_cluster_snps: int = 20,
    private_snps: int = 5,
    coverage: float = 30,
) -> dict[str, str]:
    """Simulate an outbreak cohort with shared cluster and private mutations."""

    if spec.exercise != "outbreak":
        raise ValueError("native outbreak generation requires exercise = 'outbreak'")
    destination = Path(release_dir).resolve()
    work = destination / "build" / "work" / "outbreak"
    work.mkdir(parents=True, exist_ok=True)
    metadata = read_metadata(metadata_csv)
    resolved_samples = resolve_release_samples(spec, id_salt=id_salt)
    required_sources = {sample.source_id for sample in resolved_samples}
    missing = required_sources - set(metadata)
    if missing:
        raise ValueError(
            "outbreak metadata is missing spec sources: {0}".format(
                ", ".join(sorted(missing))
            )
        )

    reference_path = Path(base_genome).resolve()
    reference = _reference_sequence(reference_path)
    cluster_mutations: dict[str, dict[int, str]] = {}
    for sample in resolved_samples:
        row = metadata[sample.source_id]
        cluster = row.get("Cluster", "").strip()
        if cluster and cluster not in cluster_mutations:
            cluster_mutations[cluster] = _mutations(
                reference,
                shared_cluster_snps,
                _stable_seed(spec.release_id, "cluster", cluster, spec.master_seed),
                set(),
            )

    def simulate_sample(sample):
        source_id = sample.source_id
        row = metadata[source_id]
        cluster = row.get("Cluster", "").strip()
        shared = cluster_mutations.get(cluster, {})
        private = _mutations(
            reference,
            private_snps,
            _stable_seed(spec.release_id, source_id, spec.master_seed),
            set(shared),
        )
        mutations = {**shared, **private}
        sample_fasta = work / "{0}.fasta".format(source_id)
        sample_fasta.write_text(
            ">{0}\n{1}\n".format(sample.sample_id, _apply(reference, mutations)),
            encoding="utf-8",
        )
        r1 = work / "{0}_R1.fastq.gz".format(source_id)
        r2 = work / "{0}_R2.fastq.gz".format(source_id)
        _simulate_short_reads(
            sample_fasta,
            r1,
            r2,
            seed=_stable_seed(spec.release_id, source_id, "reads", spec.master_seed),
            coverage=coverage,
            read_length=150,
            fragment_length=300,
            fragment_sd=50,
        )
        return {
            "source_id": source_id,
            "cluster": cluster,
            "shared_mutations": len(shared),
            "private_mutations": len(private),
            "total_mutations": len(mutations),
        }

    allocated_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    workers = min(max(1, allocated_cpus), len(resolved_samples))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        simulation_rows = list(executor.map(simulate_sample, resolved_samples))
    json_dump(
        destination / "build" / "outbreak_simulation.json",
        {
            "release_id": spec.release_id,
            "base_genome": str(reference_path),
            "base_genome_sha256": hashlib.sha256(reference.encode("ascii")).hexdigest(),
            "shared_cluster_snps": shared_cluster_snps,
            "private_snps": private_snps,
            "coverage": coverage,
            "samples": simulation_rows,
        },
    )
    return build_outbreak_release(
        spec,
        work,
        metadata_csv,
        destination,
        id_salt=id_salt,
        preanonymized=True,
    )
