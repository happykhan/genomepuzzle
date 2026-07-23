"""Assessment-side calibration for generated pilot releases.

This module is intentionally separate from release sealing. It analyses the
participant-facing files as a participant might, while writing results under
the digest-excluded ``build/calibration`` directory.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable

from genomepuzzle.contract import json_dump, sha256_file, utc_now
from genomepuzzle.runtime import require_tool


def fastq_metrics(path: str | os.PathLike[str]) -> dict[str, int | float]:
    """Stream one FASTQ and return deterministic basic read metrics."""

    records = bases = quality_sum = gc_bases = 0
    minimum_length: int | None = None
    maximum_length = 0
    with gzip.open(path, "rt", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            separator = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not sequence or not separator or not quality:
                raise ValueError("truncated FASTQ record in {0}".format(path))
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError("malformed FASTQ record in {0}".format(path))
            if len(sequence) != len(quality):
                raise ValueError("sequence/quality length mismatch in {0}".format(path))
            length = len(sequence)
            records += 1
            bases += length
            gc_bases += sum(base in "GCgc" for base in sequence)
            quality_sum += sum(ord(value) - 33 for value in quality)
            minimum_length = (
                length if minimum_length is None else min(minimum_length, length)
            )
            maximum_length = max(maximum_length, length)
    if not records or not bases:
        raise ValueError("FASTQ contains no reads: {0}".format(path))
    return {
        "records": records,
        "bases": bases,
        "minimum_length": int(minimum_length or 0),
        "maximum_length": maximum_length,
        "mean_length": round(bases / records, 3),
        "mean_quality": round(quality_sum / bases, 3),
        "gc_fraction": round(gc_bases / bases, 6),
    }


def fasta_metrics(path: str | os.PathLike[str]) -> dict[str, int | float]:
    """Return assembly size, contiguity and ambiguous-base measurements."""

    lengths: list[int] = []
    ambiguous = 0
    current = 0
    with open(path, encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            sequence = line.strip()
            current += len(sequence)
            ambiguous += sum(base not in "ACGTacgt" for base in sequence)
    if current:
        lengths.append(current)
    if not lengths:
        raise ValueError("FASTA contains no sequences: {0}".format(path))
    total = sum(lengths)
    threshold = total / 2
    cumulative = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            n50 = length
            break
    return {
        "contigs": len(lengths),
        "total_bases": total,
        "n50": n50,
        "maximum_contig": max(lengths),
        "ambiguous_bases": ambiguous,
        "ambiguous_fraction": round(ambiguous / total, 8),
    }


def _run(command: Iterable[str], commands: list[str], *, stdout: Path | None = None) -> None:
    values = [str(value) for value in command]
    commands.append(shlex.join(values))
    if stdout is None:
        subprocess.run(values, check=True)
        return
    with open(stdout, "w", encoding="utf-8") as handle:
        subprocess.run(values, check=True, stdout=handle)


def _load_release(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    if not (root / "COMPLETE.json").is_file():
        raise ValueError("pilot review requires a sealed release")
    manifest = json.loads((root / "public" / "manifest.json").read_text())
    provenance = json.loads((root / "private" / "provenance.json").read_text())
    return manifest, provenance


def calibrate_release(
    release_dir: str | os.PathLike[str],
    *,
    reference_dir: str | os.PathLike[str],
    sample_ids: Iterable[str] = (),
    threads: int = 8,
    memory_gb: int = 32,
) -> Path:
    """Run read QC, participant-like assemblies, typing and relatedness checks."""

    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("pilot calibration is heavy and must run inside SLURM")
    root = Path(release_dir).resolve()
    manifest, provenance = _load_release(root)
    exercise = manifest["exercise"]
    if exercise not in {"assembly", "hybrid", "outbreak"}:
        raise ValueError("pilot calibration requires an assembly, hybrid or outbreak release")

    output = root / "build" / "calibration"
    if output.exists():
        shutil.rmtree(output)
    assemblies = output / "assemblies"
    mash_output = output / "mash"
    assemblies.mkdir(parents=True)
    mash_output.mkdir()
    commands: list[str] = []
    sample_truth = {item["sample_id"]: item for item in provenance["samples"]}
    requested = set(sample_ids)
    available = {item["sample_id"] for item in manifest["samples"]}
    unknown = requested - available
    if unknown:
        raise ValueError(
            "unknown calibration sample IDs: {0}".format(", ".join(sorted(unknown)))
        )
    selected_samples = [
        item
        for item in manifest["samples"]
        if not requested or item["sample_id"] in requested
    ]

    reference_path = Path(reference_dir).resolve()
    source_paths = []
    for item in provenance["samples"]:
        source_id = item["source_id"]
        for suffix in ("", ".fasta", ".fna", ".fa"):
            candidate = reference_path / "{0}{1}".format(source_id, suffix)
            if candidate.is_file():
                source_paths.append(str(candidate))
                break
    source_paths = sorted(set(source_paths))

    mash = require_tool("mash")
    reference_sketch = mash_output / "sources"
    if source_paths:
        _run(
            [mash, "sketch", "-o", str(reference_sketch), *source_paths],
            commands,
        )

    sample_reports: list[dict[str, Any]] = []
    contigs: list[str] = []
    for sample in selected_samples:
        sample_id = sample["sample_id"]
        roles = sample["files"]
        paths = {
            role: root / "public" / "files" / details["filename"]
            for role, details in roles.items()
        }
        report: dict[str, Any] = {
            "sample_id": sample_id,
            "implant": sample_truth[sample_id]["implant"]["type"],
            "read_metrics": {
                role: fastq_metrics(path)
                for role, path in paths.items()
                if role in {"read_1", "read_2", "long_reads"}
            },
        }
        if {"read_1", "read_2"} <= paths.keys():
            if (
                report["read_metrics"]["read_1"]["records"]
                != report["read_metrics"]["read_2"]["records"]
            ):
                raise ValueError("paired read counts differ for {0}".format(sample_id))
            if source_paths:
                screen_path = mash_output / "{0}.screen.tsv".format(sample_id)
                _run(
                    [
                        mash,
                        "screen",
                        str(reference_sketch) + ".msh",
                        str(paths["read_1"]),
                        str(paths["read_2"]),
                    ],
                    commands,
                    stdout=screen_path,
                )
                report["mash_screen"] = str(screen_path.relative_to(root))

            assembly_dir = assemblies / sample_id
            command = [
                require_tool("spades.py"),
                "--careful",
                "-t",
                str(threads),
                "-m",
                str(memory_gb),
                "-1",
                str(paths["read_1"]),
                "-2",
                str(paths["read_2"]),
                "-o",
                str(assembly_dir),
            ]
            if "long_reads" in paths:
                command.extend(["--nanopore", str(paths["long_reads"])])
            _run(command, commands)
            contig_path = assembly_dir / "contigs.fasta"
            report["assembly_metrics"] = fasta_metrics(contig_path)
            report["assembly"] = str(contig_path.relative_to(root))
            contigs.append(str(contig_path))
        sample_reports.append(report)

    kleborate_output = output / "kleborate"
    _run(
        [require_tool("kleborate"), "-a", *contigs, "-o", str(kleborate_output)],
        commands,
    )
    tree_path: Path | None = None
    if exercise == "outbreak":
        tree_path = output / "mashtree.dnd"
        _run(
            [
                require_tool("mashtree"),
                "--numcpus",
                str(threads),
                "--outtree",
                str(tree_path),
                *contigs,
            ],
            commands,
        )

    report_path = output / "report.json"
    repository = Path(__file__).resolve().parents[1]
    git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    json_dump(
        report_path,
        {
            "schema_version": "1.0",
            "release_id": manifest["release_id"],
            "exercise": exercise,
            "slurm_job_id": os.environ["SLURM_JOB_ID"],
            "git_commit": git_commit,
            "pixi_lock_sha256": sha256_file(repository / "pixi.lock"),
            "completed_at": utc_now(),
            "commands": commands,
            "samples": sample_reports,
            "kleborate_output": str(kleborate_output.relative_to(root)),
            "tree": str(tree_path.relative_to(root)) if tree_path else None,
        },
    )
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument(
        "--sample-id",
        action="append",
        default=[],
        help="Calibrate only this public sample ID; repeat as needed.",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-gb", type=int, default=32)
    args = parser.parse_args()
    print(
        calibrate_release(
            args.release_dir,
            reference_dir=args.reference_dir,
            sample_ids=args.sample_id,
            threads=args.threads,
            memory_gb=args.memory_gb,
        )
    )


if __name__ == "__main__":
    main()
