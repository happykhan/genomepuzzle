"""Fetch and stage NCBI assembly sources named by a release specification."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path
from typing import Callable, Iterable

from genomepuzzle.release import ReleaseSpec, sha256_file


def required_assembly_accessions(spec: ReleaseSpec) -> tuple[str, ...]:
    """Return every target and contaminant assembly required by ``spec``."""

    accessions: set[str] = set()
    for sample in spec.samples:
        accessions.add(sample.source_id)
        for parameter in ("contaminant_source_id", "replacement_source_id"):
            additional = sample.implant_parameters.get(parameter)
            if additional is None:
                continue
            if not isinstance(additional, str) or not additional.strip():
                raise ValueError(
                    "implant_parameters.{0} must be a non-empty string".format(
                        parameter
                    )
                )
            accessions.add(additional.strip())
    return tuple(sorted(accessions))


def _safe_extract(package: Path, destination: Path) -> None:
    """Extract an NCBI package without permitting paths outside destination."""

    destination = destination.resolve()
    with zipfile.ZipFile(package) as archive:
        for member in archive.infolist():
            member_path = (destination / member.filename).resolve()
            if member_path != destination and destination not in member_path.parents:
                raise ValueError(
                    "NCBI package contains an unsafe path: {0}".format(member.filename)
                )
        archive.extractall(destination)


def _find_genome(package_root: Path, accession: str) -> Path:
    data_dir = package_root / "ncbi_dataset" / "data"
    candidates = sorted(data_dir.glob("{0}*/*.fna".format(accession)))
    candidates = [
        path for path in candidates if path.name != "rna.fna" and path.is_file()
    ]
    if len(candidates) != 1:
        raise ValueError(
            "expected one genomic FASTA for {0}, found {1}".format(
                accession, len(candidates)
            )
        )
    return candidates[0]


def fetch_assembly_sources(
    spec: ReleaseSpec,
    output_dir: str | Path,
    datasets_executable: str = "datasets",
    refresh: bool = False,
    runner: Callable[..., object] = subprocess.run,
) -> Path:
    """Download and stage source FASTAs with stable filenames and checksums.

    Existing staged FASTAs are reused by default. ``refresh`` downloads every
    source again and replaces files only after the complete package validates.
    """

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    accessions = required_assembly_accessions(spec)
    missing = [
        accession
        for accession in accessions
        if refresh or not (destination / "{0}.fasta".format(accession)).is_file()
    ]

    if missing:
        with tempfile.TemporaryDirectory(prefix="genomepuzzle-sources-") as temporary:
            temporary_path = Path(temporary)
            package = temporary_path / "ncbi_dataset.zip"
            extracted = temporary_path / "extracted"
            runner(
                [
                    datasets_executable,
                    "download",
                    "genome",
                    "accession",
                    *missing,
                    "--include",
                    "genome",
                    "--filename",
                    str(package),
                    "--no-progressbar",
                ],
                check=True,
            )
            if not package.is_file():
                raise ValueError(
                    "NCBI datasets did not create package: {0}".format(package)
                )
            _safe_extract(package, extracted)
            staged: list[tuple[Path, Path]] = []
            for accession in missing:
                source = _find_genome(extracted, accession)
                temporary_fasta = destination / ".{0}.fasta.tmp".format(accession)
                shutil.copyfile(source, temporary_fasta)
                if temporary_fasta.stat().st_size == 0:
                    raise ValueError(
                        "downloaded FASTA is empty for {0}".format(accession)
                    )
                staged.append(
                    (temporary_fasta, destination / "{0}.fasta".format(accession))
                )
            for temporary_fasta, final_fasta in staged:
                temporary_fasta.replace(final_fasta)

    manifest = {
        "schema_version": "1.0",
        "release_id": spec.release_id,
        "sources": [
            {
                "source_id": accession,
                "file": "{0}.fasta".format(accession),
                "sha256": sha256_file(
                    destination / "{0}.fasta".format(accession)
                ),
            }
            for accession in accessions
        ],
    }
    manifest_path = destination / "sources.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest_path
