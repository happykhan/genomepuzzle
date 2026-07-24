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
    require_available_release_directory,
    sha256_file,
    write_release_manifests,
)
from genomepuzzle.contract import (
    complete_release,
    failure_reason_for_implant,
    retained_read_pairs_for_fault,
)
from genomepuzzle.sequence_io import (
    anonymize_paired_fastq_in_place,
    anonymize_single_fastq_in_place,
    validate_paired_fastq,
    validate_single_fastq,
)


OUTBREAK_IMPLANTS = {
    "NORMAL",
    "NONE",
    "LOW_COVERAGE",
    "CONTAMINATED",
    "ZERO_BYTE_R1",
    "ZERO_BYTE_R2",
    "MISSING_R1",
    "MISSING_R2",
    "TEN_READ_PAIRS",
    "TRUNCATE_TO_READ_PAIRS",
    "WRONG_ORGANISM",
}
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


def _copy_first_pairs(
    r1: Path,
    r2: Path,
    output_r1: Path,
    output_r2: Path,
    count: int,
) -> int:
    if count < 1:
        raise ValueError("pair count must be positive")
    copied = 0
    with gzip.open(r1, "rt") as input_r1, gzip.open(
        r2, "rt"
    ) as input_r2, gzip.open(output_r1, "wt") as target_r1, gzip.open(
        output_r2, "wt"
    ) as target_r2:
        while copied < count:
            record_r1 = [input_r1.readline() for _ in range(4)]
            record_r2 = [input_r2.readline() for _ in range(4)]
            if not record_r1[0] and not record_r2[0]:
                break
            if any(not line for line in record_r1 + record_r2):
                raise ValueError("truncated or unpaired FASTQ record")
            target_r1.writelines(record_r1)
            target_r2.writelines(record_r2)
            copied += 1
    if copied != count:
        raise ValueError(
            "source contains only {0} pairs; cannot retain {1}".format(copied, count)
        )
    return copied


def _count_pairs(r1: Path, r2: Path) -> int:
    return validate_paired_fastq(r1, r2)


def build_outbreak_release(
    spec: ReleaseSpec,
    source_dir: str | Path,
    metadata_csv: str | Path,
    release_dir: str | Path,
    id_salt: str | None = None,
    preanonymized: bool = False,
) -> dict[str, str]:
    """Build an outbreak release from frozen TreeToReads FASTQ output."""

    if spec.exercise != "outbreak":
        raise ValueError("outbreak builder requires exercise = 'outbreak'")
    destination = require_available_release_directory(release_dir)
    files_dir = destination / "public" / "files"
    files_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_metadata(metadata_csv)
    artifacts = [
        _build_outbreak_sample(
            sample,
            Path(source_dir),
            files_dir,
            metadata,
            preanonymized=preanonymized,
        )
        for sample in resolve_release_samples(spec, id_salt=id_salt)
    ]
    _write_sample_sheet(destination / "public" / "sample_sheet.csv", artifacts)
    manifests = write_release_manifests(
        destination,
        spec,
        artifacts,
        generator={"module": "genomepuzzle.outbreak", "simulator": "TreeToReads"},
    )
    complete_release(destination)
    return manifests


def _build_outbreak_sample(
    sample: ResolvedReleaseSample,
    source_dir: Path,
    files_dir: Path,
    metadata: dict[str, dict[str, str]],
    *,
    preanonymized: bool = False,
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
        "source_r1_sha256": sha256_file(source_r1),
        "source_r2_sha256": sha256_file(source_r2),
    }

    if sample.implant in {"NORMAL", "NONE"}:
        shutil.copyfile(source_r1, output_r1)
        shutil.copyfile(source_r2, output_r2)
    elif sample.implant == "LOW_COVERAGE":
        fraction = float(sample.implant_parameters.get("read_fraction", 0.02))
        source_coverage = float(
            sample.implant_parameters.get("source_coverage", 30)
        )
        expected_coverage = source_coverage * fraction
        if expected_coverage > 1.0:
            raise ValueError(
                "LOW_COVERAGE must produce no more than 1x expected coverage"
            )
        provenance["retained_pairs"] = _copy_selected_pairs(
            source_r1,
            source_r2,
            output_r1,
            output_r2,
            fraction,
            sample.random_seed,
        )
        provenance["expected_coverage"] = round(expected_coverage, 4)
    elif sample.implant in {"TEN_READ_PAIRS", "TRUNCATE_TO_READ_PAIRS"}:
        retained_pairs = retained_read_pairs_for_fault(
            sample.implant, sample.implant_parameters
        )
        provenance["retained_pairs"] = _copy_first_pairs(
            source_r1, source_r2, output_r1, output_r2, retained_pairs
        )
    elif sample.implant in {"ZERO_BYTE_R1", "ZERO_BYTE_R2"}:
        shutil.copyfile(source_r1, output_r1)
        shutil.copyfile(source_r2, output_r2)
        empty_path = output_r1 if sample.implant == "ZERO_BYTE_R1" else output_r2
        empty_path.write_bytes(b"")
        provenance["zero_byte_role"] = (
            "read_1" if sample.implant == "ZERO_BYTE_R1" else "read_2"
        )
    elif sample.implant in {"MISSING_R1", "MISSING_R2"}:
        if sample.implant != "MISSING_R1":
            shutil.copyfile(source_r1, output_r1)
        if sample.implant != "MISSING_R2":
            shutil.copyfile(source_r2, output_r2)
        provenance["missing_role"] = (
            "read_1" if sample.implant == "MISSING_R1" else "read_2"
        )
    elif sample.implant == "WRONG_ORGANISM":
        replacement_id = sample.implant_parameters.get("replacement_source_id")
        if not isinstance(replacement_id, str) or not replacement_id:
            raise ValueError(
                "WRONG_ORGANISM requires implant_parameters.replacement_source_id"
            )
        if replacement_id == sample.source_id:
            raise ValueError("WRONG_ORGANISM replacement must differ from target")
        replacement_r1, replacement_r2 = _fastq_pair(source_dir, replacement_id)
        shutil.copyfile(replacement_r1, output_r1)
        shutil.copyfile(replacement_r2, output_r2)
        provenance["replacement_source_id"] = replacement_id
        provenance["replacement_r1_sha256"] = sha256_file(replacement_r1)
        provenance["replacement_r2_sha256"] = sha256_file(replacement_r2)
        provenance["replacement_fraction"] = 1.0
    else:
        contaminant_id = sample.implant_parameters.get("contaminant_source_id")
        if not isinstance(contaminant_id, str) or not contaminant_id:
            raise ValueError(
                "CONTAMINATED requires implant_parameters.contaminant_source_id"
            )
        if contaminant_id == sample.source_id:
            raise ValueError("CONTAMINATED source must differ from target")
        contaminant_r1, contaminant_r2 = _fastq_pair(source_dir, contaminant_id)
        fraction = float(
            sample.implant_parameters.get("contamination_fraction", 0.5)
        )
        if not 0.3 <= fraction <= 0.9:
            raise ValueError(
                "contamination_fraction must be between 0.30 and 0.90"
            )
        clean_r1 = files_dir / ".{0}_clean_R1.gz".format(sample.sample_id)
        clean_r2 = files_dir / ".{0}_clean_R2.gz".format(sample.sample_id)
        temp_r1 = files_dir / ".{0}_contaminant_R1.gz".format(sample.sample_id)
        temp_r2 = files_dir / ".{0}_contaminant_R2.gz".format(sample.sample_id)
        try:
            total_pairs = min(
                _count_pairs(source_r1, source_r2),
                _count_pairs(contaminant_r1, contaminant_r2),
            )
            contaminant_pairs = round(total_pairs * fraction)
            clean_pairs = total_pairs - contaminant_pairs
            provenance["clean_pairs"] = _copy_first_pairs(
                source_r1,
                source_r2,
                clean_r1,
                clean_r2,
                clean_pairs,
            )
            provenance["contaminant_pairs"] = _copy_first_pairs(
                contaminant_r1,
                contaminant_r2,
                temp_r1,
                temp_r2,
                contaminant_pairs,
            )
            _append_gzip_members([clean_r1, temp_r1], output_r1)
            _append_gzip_members([clean_r2, temp_r2], output_r2)
        finally:
            clean_r1.unlink(missing_ok=True)
            clean_r2.unlink(missing_ok=True)
            temp_r1.unlink(missing_ok=True)
            temp_r2.unlink(missing_ok=True)
        achieved_fraction = contaminant_pairs / total_pairs
        provenance["requested_contamination_fraction"] = fraction
        provenance["achieved_contamination_fraction"] = round(
            achieved_fraction, 6
        )
        provenance["contaminant_source_id"] = contaminant_id
        provenance["contaminant_r1_sha256"] = sha256_file(contaminant_r1)
        provenance["contaminant_r2_sha256"] = sha256_file(contaminant_r2)

    files = {}
    faulted_pair = sample.implant in {
        "ZERO_BYTE_R1",
        "ZERO_BYTE_R2",
        "MISSING_R1",
        "MISSING_R2",
    }
    if not faulted_pair and preanonymized and sample.implant != "CONTAMINATED":
        provenance["participant_read_pairs"] = validate_paired_fastq(
            output_r1, output_r2, sample.sample_id
        )
        files = {"read_1": str(output_r1), "read_2": str(output_r2)}
    elif not faulted_pair:
        # A contaminated sample combines reads generated under two public IDs,
        # so its final headers must still be normalised to the target ID.
        provenance["participant_read_pairs"] = anonymize_paired_fastq_in_place(
            output_r1, output_r2, sample.sample_id
        )
        files = {"read_1": str(output_r1), "read_2": str(output_r2)}
    else:
        counts = {}
        for role, path, zero_fault in (
            ("read_1", output_r1, sample.implant == "ZERO_BYTE_R1"),
            ("read_2", output_r2, sample.implant == "ZERO_BYTE_R2"),
        ):
            if not path.is_file():
                continue
            if zero_fault:
                if path.stat().st_size != 0:
                    raise ValueError("zero-byte outbreak fault was not materialized")
                counts[role] = 0
            elif preanonymized:
                counts[role] = validate_single_fastq(path, sample.sample_id)
            else:
                counts[role] = anonymize_single_fastq_in_place(
                    path, sample.sample_id, role=role
                )
            files[role] = str(path)
        provenance["participant_read_counts"] = counts

    row = metadata[sample.source_id]
    expected = {
        "cluster": row.get("Cluster", ""),
        "species": row.get("SPECIES", ""),
        "qc_status": "FAIL"
        if sample.implant not in {"NORMAL", "NONE"}
        else "PASS",
        "failure_reason": failure_reason_for_implant("outbreak", sample.implant),
    }
    public_metadata = {
        key: value
        for key, value in row.items()
        if key not in PRIVATE_METADATA_FIELDS and key != "Sample"
    }
    provenance["validation"] = {
        "status": "passed",
        "checks": [
            "paired_fastq_structure",
            "anonymous_headers",
            "implant_materialized",
        ],
        "fault_type": sample.implant,
    }
    return ReleaseArtifactSample(
        sample_id=sample.sample_id,
        source_id=sample.source_id,
        random_seed=sample.random_seed,
        files=files,
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
            read_1 = sample.files.get("read_1")
            read_2 = sample.files.get("read_2")
            row["R1"] = Path(read_1).name if read_1 else ""
            row["R2"] = Path(read_2).name if read_2 else ""
            writer.writerow(row)
