"""Package TreeToReads output as reproducible outbreak assessment releases."""

from __future__ import annotations

import csv
import gzip
import random
import shutil
from pathlib import Path
from typing import Iterable

from genomepuzzle.release import (
    ReleaseArtifactSample,
    ReleaseSpec,
    ResolvedReleaseSample,
    resolve_release_samples,
    write_release_manifests,
)


OUTBREAK_IMPLANTS = {"NORMAL", "NONE", "LOW_COVERAGE", "CONTAMINATED"}
PRIVATE_METADATA_FIELDS = {"Cluster", "SPECIES"}


def _fastq_pair(source_dir: Path, source_id: str) -> tuple[Path, Path]:
    patterns = [
        ("{0}_R1.fastq.gz", "{0}_R2.fastq.gz"),
        ("{0}_1.fq.gz", "{0}_2.fq.gz"),
        ("{0}_1.fastq.gz", "{0}_2.fastq.gz"),
    ]
    for r1_pattern, r2_pattern in patterns:
        r1 = source_dir / r1_pattern.format(source_id)
        r2 = source_dir / r2_pattern.format(source_id)
        if r1.is_file() and r2.is_file():
            return r1, r2
    raise ValueError("paired FASTQ source not found for {0}".format(source_id))


def read_metadata(path: str | Path) -> dict[str, dict[str, str]]:
    with open(path, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "Sample" not in rows[0]:
        raise ValueError("outbreak metadata must contain a Sample column")
    result = {}
    for row in rows:
        sample = row.get("Sample", "").strip()
        if not sample or sample in result:
            raise ValueError("outbreak metadata has blank or duplicate Sample")
        result[sample] = row
    return result


def _copy_selected_pairs(
    r1: Path,
    r2: Path,
    output_r1: Path,
    output_r2: Path,
    fraction: float,
    seed: int,
) -> int:
    if not 0 < fraction <= 1:
        raise ValueError("read fraction must be greater than 0 and at most 1")
    rng = random.Random(seed)
    kept = 0
    with gzip.open(r1, "rt") as input_r1, gzip.open(
        r2, "rt"
    ) as input_r2, gzip.open(output_r1, "wt") as target_r1, gzip.open(
        output_r2, "wt"
    ) as target_r2:
        while True:
            record_r1 = [input_r1.readline() for _ in range(4)]
            record_r2 = [input_r2.readline() for _ in range(4)]
            if not record_r1[0] and not record_r2[0]:
                break
            if not record_r1[0] or not record_r2[0]:
                raise ValueError("paired FASTQ files have different record counts")
            if any(not line for line in record_r1 + record_r2):
                raise ValueError("truncated FASTQ record")
            if rng.random() <= fraction:
                target_r1.writelines(record_r1)
                target_r2.writelines(record_r2)
                kept += 1
    if kept == 0:
        raise ValueError("implant removed every read; increase fraction or source size")
    return kept


def _append_gzip_members(inputs: Iterable[Path], output: Path) -> None:
    with open(output, "wb") as target:
        for path in inputs:
            with open(path, "rb") as source:
                shutil.copyfileobj(source, target)


def build_outbreak_release(
    spec: ReleaseSpec,
    source_dir: str | Path,
    metadata_csv: str | Path,
    release_dir: str | Path,
    id_salt: str | None = None,
) -> dict[str, str]:
    """Build an outbreak release from frozen TreeToReads FASTQ output."""

    if spec.exercise != "outbreak":
        raise ValueError("outbreak builder requires exercise = 'outbreak'")
    destination = Path(release_dir)
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("release directory is not empty: {0}".format(destination))
    files_dir = destination / "public" / "files"
    files_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_metadata(metadata_csv)
    artifacts = [
        _build_outbreak_sample(sample, Path(source_dir), files_dir, metadata)
        for sample in resolve_release_samples(spec, id_salt=id_salt)
    ]
    _write_sample_sheet(destination / "public" / "sample_sheet.csv", artifacts)
    manifests = write_release_manifests(
        destination,
        spec,
        artifacts,
        generator={"module": "genomepuzzle.outbreak", "simulator": "TreeToReads"},
    )
    (destination / "COMPLETE").write_text("release complete\n", encoding="utf-8")
    return manifests


def _build_outbreak_sample(
    sample: ResolvedReleaseSample,
    source_dir: Path,
    files_dir: Path,
    metadata: dict[str, dict[str, str]],
) -> ReleaseArtifactSample:
    if sample.implant not in OUTBREAK_IMPLANTS:
        raise ValueError("unsupported outbreak implant: {0}".format(sample.implant))
    if sample.source_id not in metadata:
        raise ValueError("metadata missing source sample {0}".format(sample.source_id))
    source_r1, source_r2 = _fastq_pair(source_dir, sample.source_id)
    output_r1 = files_dir / "{0}_R1.fastq.gz".format(sample.sample_id)
    output_r2 = files_dir / "{0}_R2.fastq.gz".format(sample.sample_id)
    provenance: dict[str, object] = {
        "source_r1": str(source_r1),
        "source_r2": str(source_r2),
    }

    if sample.implant in {"NORMAL", "NONE"}:
        shutil.copyfile(source_r1, output_r1)
        shutil.copyfile(source_r2, output_r2)
    elif sample.implant == "LOW_COVERAGE":
        fraction = float(sample.implant_parameters.get("read_fraction", 0.1))
        provenance["retained_pairs"] = _copy_selected_pairs(
            source_r1,
            source_r2,
            output_r1,
            output_r2,
            fraction,
            sample.random_seed,
        )
    else:
        contaminant_id = sample.implant_parameters.get("contaminant_source_id")
        if not isinstance(contaminant_id, str) or not contaminant_id:
            raise ValueError(
                "CONTAMINATED requires implant_parameters.contaminant_source_id"
            )
        contaminant_r1, contaminant_r2 = _fastq_pair(source_dir, contaminant_id)
        fraction = float(
            sample.implant_parameters.get("contamination_fraction", 0.1)
        )
        temp_r1 = files_dir / ".{0}_contaminant_R1.gz".format(sample.sample_id)
        temp_r2 = files_dir / ".{0}_contaminant_R2.gz".format(sample.sample_id)
        try:
            provenance["contaminant_pairs"] = _copy_selected_pairs(
                contaminant_r1,
                contaminant_r2,
                temp_r1,
                temp_r2,
                fraction,
                sample.random_seed,
            )
            _append_gzip_members([source_r1, temp_r1], output_r1)
            _append_gzip_members([source_r2, temp_r2], output_r2)
        finally:
            temp_r1.unlink(missing_ok=True)
            temp_r2.unlink(missing_ok=True)
        provenance["contaminant_source_id"] = contaminant_id

    row = metadata[sample.source_id]
    expected = {
        "cluster": row.get("Cluster", ""),
        "species": row.get("SPECIES", ""),
        "qc_decision": "exclude"
        if sample.implant not in {"NORMAL", "NONE"}
        else "include",
    }
    public_metadata = {
        key: value
        for key, value in row.items()
        if key not in PRIVATE_METADATA_FIELDS and key != "Sample"
    }
    return ReleaseArtifactSample(
        sample_id=sample.sample_id,
        source_id=sample.source_id,
        random_seed=sample.random_seed,
        files={"read_1": str(output_r1), "read_2": str(output_r2)},
        expected_answers=expected,
        implant=sample.implant,
        implant_parameters=sample.implant_parameters,
        public_metadata=public_metadata,
        private_provenance=provenance,
    )


def _write_sample_sheet(path: Path, samples: list[ReleaseArtifactSample]) -> None:
    metadata_fields = sorted(
        {key for sample in samples for key in sample.public_metadata.keys()}
    )
    fieldnames = ["Sample"] + metadata_fields + ["R1", "R2"]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample in samples:
            row = {"Sample": sample.sample_id}
            row.update(sample.public_metadata)
            row["R1"] = Path(sample.files["read_1"]).name
            row["R2"] = Path(sample.files["read_2"]).name
            writer.writerow(row)
