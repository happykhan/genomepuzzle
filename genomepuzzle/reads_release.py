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
        source_r1 = _find_asset(
            source_path, sample.source_id, ("_R1.fastq.gz", "_1.fastq.gz", "_1.fq.gz")
        )
        source_r2 = _find_asset(
            source_path, sample.source_id, ("_R2.fastq.gz", "_2.fastq.gz", "_2.fq.gz")
        )
        output_r1 = files_dir / "{0}_R1.fastq.gz".format(sample.sample_id)
        output_r2 = files_dir / "{0}_R2.fastq.gz".format(sample.sample_id)
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
        files = {"read_1": str(output_r1), "read_2": str(output_r2)}
        provenance: dict[str, object] = {
            "source_r1": str(source_r1),
            "source_r2": str(source_r2),
            "source_r1_sha256": sha256_file(source_r1),
            "source_r2_sha256": sha256_file(source_r2),
            "participant_read_pairs": participant_pairs,
        }
        validation = dict((implant_validations or {}).get(sample.source_id, {}))
        if sample.implant in {"NORMAL", "NONE"} and not validation:
            validation = {
                "status": "passed",
                "checks": ["paired_fastq_structure", "anonymous_headers"],
                "implant": sample.implant,
            }
        if validation.get("status") != "passed":
            raise ValueError(
                "troublesome sample {0} requires a passing implant validation".format(
                    sample.source_id
                )
            )
        provenance["validation"] = validation
        if spec.exercise == "hybrid":
            source_long = _find_asset(
                source_path,
                sample.source_id,
                ("_long.fastq.gz", "_ONT.fastq.gz", "_long.fq.gz"),
            )
            output_long = files_dir / "{0}_long.fastq.gz".format(sample.sample_id)
            shutil.copyfile(source_long, output_long)
            if preanonymized:
                provenance["participant_long_reads"] = validate_single_fastq(
                    output_long, sample.sample_id
                )
            else:
                provenance["participant_long_reads"] = anonymize_single_fastq_in_place(
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
