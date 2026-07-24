"""Package final short-read or hybrid assets under the shared release contract."""

from __future__ import annotations

import csv
import json
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Mapping

from genomepuzzle.release import (
    ReleaseArtifactSample,
    ReleaseSpec,
    resolve_release_samples,
    require_available_release_directory,
    sha256_file,
    write_release_manifests,
)
from genomepuzzle.contract import complete_release
from genomepuzzle.sequence_io import (
    anonymize_paired_fastq,
    anonymize_single_fastq_in_place,
    validate_paired_fastq,
    validate_single_fastq,
)


def load_expected_answers(path: str | Path) -> dict[str, Mapping[str, object]]:
    """Load a private JSON object keyed by source ID."""

    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict) or not payload:
        raise ValueError("expected answers must be a non-empty JSON object")
    for source_id, answers in payload.items():
        if not isinstance(source_id, str) or not isinstance(answers, dict):
            raise ValueError("expected answers must map source IDs to objects")
    return payload


def _find_asset(source_dir: Path, source_id: str, suffixes: tuple[str, ...]) -> Path:
    for suffix in suffixes:
        candidate = source_dir / "{0}{1}".format(source_id, suffix)
        if candidate.is_file():
            return candidate
    raise ValueError(
        "missing source asset for {0}: expected one of {1}".format(
            source_id, ", ".join(source_id + suffix for suffix in suffixes)
        )
    )


def _find_asset_optional(
    source_dir: Path, source_id: str, suffixes: tuple[str, ...]
) -> Path | None:
    for suffix in suffixes:
        candidate = source_dir / "{0}{1}".format(source_id, suffix)
        if candidate.is_file():
            return candidate
    return None


def package_read_release(
    spec: ReleaseSpec,
    source_dir: str | Path,
    expected_answers: Mapping[str, Mapping[str, object]] | None,
    release_dir: str | Path,
    id_salt: str | None = None,
    implant_validations: Mapping[str, Mapping[str, object]] | None = None,
    preanonymized: bool = False,
) -> dict[str, str]:
    """Package final, already-implanted reads for assembly or hybrid exercises."""

    if spec.exercise not in {"assembly", "hybrid"}:
        raise ValueError("read packager requires exercise = 'assembly' or 'hybrid'")
    destination = require_available_release_directory(release_dir)
    files_dir = destination / "public" / "files"
    files_dir.mkdir(parents=True, exist_ok=True)
    source_path = Path(source_dir)
    answers_by_source = dict(expected_answers or {})
    for sample_spec in spec.samples:
        if sample_spec.expected_answers:
            answers_by_source.setdefault(
                sample_spec.source_id, sample_spec.expected_answers
            )
    samples = list(resolve_release_samples(spec, id_salt=id_salt))

    def package_sample(sample) -> ReleaseArtifactSample:
        if sample.source_id not in answers_by_source:
            raise ValueError("expected answers missing source {0}".format(sample.source_id))
        r1_suffixes = ("_R1.fastq.gz", "_1.fastq.gz", "_1.fq.gz")
        r2_suffixes = ("_R2.fastq.gz", "_2.fastq.gz", "_2.fq.gz")
        source_r1 = _find_asset_optional(source_path, sample.source_id, r1_suffixes)
        source_r2 = _find_asset_optional(source_path, sample.source_id, r2_suffixes)
        missing_r1 = sample.implant == "MISSING_R1"
        missing_r2 = sample.implant == "MISSING_R2"
        if missing_r1 and source_r1 is not None:
            raise ValueError("MISSING_R1 fault materialized an R1 asset")
        if missing_r2 and source_r2 is not None:
            raise ValueError("MISSING_R2 fault materialized an R2 asset")
        if source_r1 is None and not missing_r1:
            source_r1 = _find_asset(source_path, sample.source_id, r1_suffixes)
        if source_r2 is None and not missing_r2:
            source_r2 = _find_asset(source_path, sample.source_id, r2_suffixes)
        output_r1 = files_dir / "{0}_R1.fastq.gz".format(sample.sample_id)
        output_r2 = files_dir / "{0}_R2.fastq.gz".format(sample.sample_id)
        zero_r1 = sample.implant == "ZERO_BYTE_R1"
        zero_r2 = sample.implant == "ZERO_BYTE_R2"
        faulted_pair = (
            missing_r1
            or missing_r2
            or zero_r1
            or zero_r2
        )
        files: dict[str, str] = {}
        provenance: dict[str, object] = {}
        if not faulted_pair:
            assert source_r1 is not None and source_r2 is not None
            if preanonymized:
                participant_pairs = validate_paired_fastq(
                    source_r1, source_r2, sample.sample_id
                )
                shutil.copyfile(source_r1, output_r1)
                shutil.copyfile(source_r2, output_r2)
            else:
                participant_pairs = anonymize_paired_fastq(
                    source_r1,
                    source_r2,
                    output_r1,
                    output_r2,
                    sample.sample_id,
                )
            files.update({"read_1": str(output_r1), "read_2": str(output_r2)})
            provenance["participant_read_pairs"] = participant_pairs
        else:
            participant_counts: dict[str, int] = {}
            for role, source, output, is_zero in (
                ("read_1", source_r1, output_r1, zero_r1),
                ("read_2", source_r2, output_r2, zero_r2),
            ):
                if source is None:
                    continue
                if is_zero:
                    if source.stat().st_size != 0:
                        raise ValueError(
                            "{0} fault did not produce a zero-byte file".format(
                                sample.implant
                            )
                        )
                    shutil.copyfile(source, output)
                    participant_counts[role] = 0
                else:
                    if source.stat().st_size == 0:
                        raise ValueError(
                            "undeclared zero-byte participant file for {0}".format(role)
                        )
                    shutil.copyfile(source, output)
                    if preanonymized:
                        participant_counts[role] = validate_single_fastq(
                            output, sample.sample_id
                        )
                    else:
                        participant_counts[role] = anonymize_single_fastq_in_place(
                            output, sample.sample_id, role=role
                        )
                files[role] = str(output)
            provenance["participant_read_counts"] = participant_counts
        for role, source in (("r1", source_r1), ("r2", source_r2)):
            if source is not None:
                provenance["source_{0}".format(role)] = str(source)
                provenance["source_{0}_sha256".format(role)] = sha256_file(source)
        validation = dict((implant_validations or {}).get(sample.source_id, {}))
        if sample.implant in {"NORMAL", "NONE"} and not validation:
            validation = {
                "status": "passed",
                "checks": ["paired_fastq_structure", "anonymous_headers"],
                "fault_type": sample.implant,
            }
        validation.setdefault("fault_type", sample.implant)
        if validation.get("status") != "passed":
            raise ValueError(
                "troublesome sample {0} requires a passing implant validation".format(
                    sample.source_id
                )
            )
        provenance["validation"] = validation
        if spec.exercise == "hybrid":
            long_suffixes = ("_long.fastq.gz", "_ONT.fastq.gz", "_long.fq.gz")
            source_long = _find_asset_optional(
                source_path, sample.source_id, long_suffixes
            )
            missing_long = sample.implant == "MISSING_LONG_READS"
            zero_long = sample.implant == "ZERO_BYTE_LONG_READS"
            if missing_long and source_long is not None:
                raise ValueError(
                    "MISSING_LONG_READS fault materialized a long-read asset"
                )
            if source_long is None and not missing_long:
                source_long = _find_asset(
                    source_path, sample.source_id, long_suffixes
                )
            if source_long is not None:
                output_long = files_dir / "{0}_long.fastq.gz".format(sample.sample_id)
                if zero_long:
                    if source_long.stat().st_size != 0:
                        raise ValueError(
                            "ZERO_BYTE_LONG_READS fault did not produce a zero-byte file"
                        )
                    shutil.copyfile(source_long, output_long)
                    provenance["participant_long_reads"] = 0
                else:
                    if source_long.stat().st_size == 0:
                        raise ValueError("undeclared zero-byte long-read file")
                    shutil.copyfile(source_long, output_long)
                    if preanonymized:
                        provenance["participant_long_reads"] = validate_single_fastq(
                            output_long, sample.sample_id
                        )
                    else:
                        provenance[
                            "participant_long_reads"
                        ] = anonymize_single_fastq_in_place(
                            output_long, sample.sample_id
                        )
                files["long_reads"] = str(output_long)
                provenance["source_long_reads"] = str(source_long)
                provenance["source_long_reads_sha256"] = sha256_file(source_long)
        return ReleaseArtifactSample(
            sample_id=sample.sample_id,
            source_id=sample.source_id,
            random_seed=sample.random_seed,
            files=files,
            expected_answers=answers_by_source[sample.source_id],
            implant=sample.implant,
            implant_parameters=sample.implant_parameters,
            public_metadata={"format": "paired FASTQ"}
            if spec.exercise == "assembly"
            else {"format": "paired short-read and long-read FASTQ"},
            private_provenance=provenance,
        )

    allocated_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    workers = min(4, max(1, allocated_cpus), len(samples))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        artifacts = list(executor.map(package_sample, samples))
    _write_sample_sheet(destination / "public" / "sample_sheet.csv", artifacts)
    manifests = write_release_manifests(
        destination,
        spec,
        artifacts,
        generator={
            "module": "genomepuzzle.reads_release",
            "input_stage": "final_implanted_reads",
        },
    )
    complete_release(destination)
    return manifests


def _write_sample_sheet(path: Path, samples: list[ReleaseArtifactSample]) -> None:
    hybrid = any("long_reads" in sample.files for sample in samples)
    fields = ["sample_id", "read_1", "read_2"]
    if hybrid:
        fields.append("long_reads")
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sample in samples:
            writer.writerow(
                {
                    role: Path(filename).name if role != "sample_id" else filename
                    for role, filename in {
                        "sample_id": sample.sample_id,
                        **sample.files,
                    }.items()
                }
            )
