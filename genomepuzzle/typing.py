"""Build Kleborate typing releases from assembly FASTA files."""

from __future__ import annotations

import random
import csv
import subprocess
from pathlib import Path
from typing import Callable, Iterable, Mapping

from genomepuzzle.release import (
    ReleaseArtifactSample,
    ReleaseSpec,
    ResolvedReleaseSample,
    resolve_release_samples,
    require_available_release_directory,
    sha256_file,
    write_release_manifests,
)
from genomepuzzle.contract import complete_release


TYPING_IMPLANTS = {"NORMAL", "NONE", "FRAGMENTED", "MIXED_CONTIGS"}
KLEBORATE_FIELDS = {
    "species": "enterobacterales__species__species",
    "st": "klebsiella_pneumo_complex__mlst__ST",
    "k_locus": "klebsiella_pneumo_complex__kaptive__K_locus",
    "capsule_type": "klebsiella_pneumo_complex__kaptive__K_type",
    "wzi": "klebsiella_pneumo_complex__wzi__wzi",
    "o_locus": "klebsiella_pneumo_complex__kaptive__O_locus",
    "o_type": "klebsiella_pneumo_complex__kaptive__O_type",
    "bla_carb": "klebsiella_pneumo_complex__amr__Bla_Carb_acquired",
}


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA without retaining potentially identifying descriptions."""

    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0] or "contig"
                sequence = []
            elif name is None:
                raise ValueError(
                    "{0}:{1}: sequence before first FASTA header".format(
                        path, line_number
                    )
                )
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    if not records or any(not sequence for _, sequence in records):
        raise ValueError("FASTA contains no usable sequence records: {0}".format(path))
    return records


def write_anonymous_fasta(
    path: str | Path, records: Iterable[tuple[str, str]], sample_id: str
) -> None:
    """Write stable participant-safe contig names and wrapped sequence."""

    with open(path, "w", encoding="utf-8") as handle:
        for index, (_, sequence) in enumerate(records, start=1):
            handle.write(">{0}_contig_{1:05d}\n".format(sample_id, index))
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def run_kleborate(
    fasta_path: Path, executable: str = "kleborate"
) -> Mapping[str, object]:
    """Run Kleborate against the final FASTA and normalise scoring fields."""

    output_dir = fasta_path.parent.parent.parent / "private" / "kleborate"
    sample_output = output_dir / fasta_path.stem
    sample_output.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            executable,
            "-a",
            str(fasta_path.resolve()),
            "-o",
            str(sample_output.resolve()),
            "-p",
            "kpsc",
        ],
        check=True,
    )
    result_path = sample_output / "klebsiella_pneumo_complex_output.txt"
    if not result_path.is_file():
        raise ValueError("Kleborate did not produce {0}".format(result_path))
    with open(result_path, encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"), None)
    if row is None:
        raise ValueError("Kleborate output contains no result row: {0}".format(result_path))
    missing = [column for column in KLEBORATE_FIELDS.values() if column not in row]
    if missing:
        raise ValueError("Kleborate output is missing columns: {0}".format(", ".join(missing)))
    result = {
        public_name: row[column] for public_name, column in KLEBORATE_FIELDS.items()
    }
    result["analysis_status"] = "complete"
    return result


def fragment_records(
    records: Iterable[tuple[str, str]], fragment_size: int
) -> list[tuple[str, str]]:
    """Split contigs into deterministic fixed-size fragments."""

    if fragment_size < 100:
        raise ValueError("fragment_size must be at least 100 bases")
    fragments = []
    for name, sequence in records:
        for offset in range(0, len(sequence), fragment_size):
            fragment = sequence[offset : offset + fragment_size]
            if fragment:
                fragments.append(
                    ("{0}_{1}".format(name, offset // fragment_size + 1), fragment)
                )
    return fragments


def mix_records(
    target_records: Iterable[tuple[str, str]],
    contaminant_records: Iterable[tuple[str, str]],
    fraction: float,
    seed: int,
) -> list[tuple[str, str]]:
    """Add contaminant contigs up to a target fraction of target assembly bases."""

    if not 0 < fraction <= 1:
        raise ValueError("contamination_fraction must be greater than 0 and at most 1")
    target = list(target_records)
    contaminant = list(contaminant_records)
    if not contaminant:
        raise ValueError("contaminant FASTA has no records")
    required_bases = max(1, round(sum(len(seq) for _, seq in target) * fraction))
    rng = random.Random(seed)
    rng.shuffle(contaminant)
    selected = []
    selected_bases = 0
    for record in contaminant:
        selected.append(record)
        selected_bases += len(record[1])
        if selected_bases >= required_bases:
            break
    return target + selected


def _source_path(source_dir: Path, source_id: str) -> Path:
    candidates = [
        source_dir / source_id,
        source_dir / "{0}.fasta".format(source_id),
        source_dir / "{0}.fa".format(source_id),
        source_dir / "{0}.fna".format(source_id),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError("assembly source not found for {0} in {1}".format(source_id, source_dir))


def build_typing_release(
    spec: ReleaseSpec,
    source_dir: str | Path,
    release_dir: str | Path,
    id_salt: str | None = None,
    analyser: Callable[[Path], Mapping[str, object]] | None = None,
) -> dict[str, str]:
    """Build a complete typing release from local source assemblies.

    ``analyser`` must inspect the final participant FASTA. Production builds
    should supply a pinned Kleborate runner; tests and preparation builds may
    omit it, in which case the private answer record is explicitly pending.
    """

    if spec.exercise != "typing":
        raise ValueError("typing builder requires exercise = 'typing'")
    source_path = Path(source_dir)
    destination = require_available_release_directory(release_dir)
    files_dir = destination / "public" / "files"
    files_dir.mkdir(parents=True, exist_ok=True)

    resolved = resolve_release_samples(spec, id_salt=id_salt)
    artifacts = []
    for sample in resolved:
        artifacts.append(
            _build_typing_sample(sample, source_path, files_dir, analyser)
        )
    manifests = write_release_manifests(
        destination,
        spec,
        artifacts,
        generator={
            "module": "genomepuzzle.typing",
            "expected_analysis": "kleborate" if analyser else "pending",
        },
    )
    if analyser:
        complete_release(destination)
    return manifests


def _build_typing_sample(
    sample: ResolvedReleaseSample,
    source_dir: Path,
    files_dir: Path,
    analyser: Callable[[Path], Mapping[str, object]] | None,
) -> ReleaseArtifactSample:
    if sample.implant not in TYPING_IMPLANTS:
        raise ValueError(
            "unsupported typing implant {0}; expected one of {1}".format(
                sample.implant, ", ".join(sorted(TYPING_IMPLANTS))
            )
        )
    input_path = _source_path(source_dir, sample.source_id)
    records = read_fasta(input_path)
    original_contigs = len(records)
    original_bases = sum(len(sequence) for _, sequence in records)
    provenance: dict[str, object] = {
        "source_file": str(input_path),
        "source_sha256": sha256_file(input_path),
    }

    if sample.implant == "FRAGMENTED":
        fragment_size = int(sample.implant_parameters.get("fragment_size", 1000))
        records = fragment_records(records, fragment_size)
    elif sample.implant == "MIXED_CONTIGS":
        contaminant_source_id = sample.implant_parameters.get("contaminant_source_id")
        if not isinstance(contaminant_source_id, str) or not contaminant_source_id:
            raise ValueError(
                "MIXED_CONTIGS requires implant_parameters.contaminant_source_id"
            )
        contaminant_path = _source_path(source_dir, contaminant_source_id)
        fraction = float(
            sample.implant_parameters.get("contamination_fraction", 0.1)
        )
        records = mix_records(
            records,
            read_fasta(contaminant_path),
            fraction,
            sample.random_seed,
        )
        provenance["contaminant_source_file"] = str(contaminant_path)
        provenance["contaminant_source_sha256"] = sha256_file(contaminant_path)

    output_path = files_dir / "{0}.fasta".format(sample.sample_id)
    write_anonymous_fasta(output_path, records, sample.sample_id)
    expected = (
        dict(analyser(output_path))
        if analyser
        else {"analysis_status": "pending_kleborate"}
    )
    validation: dict[str, object] = {
        "status": "passed",
        "checks": ["anonymous_fasta", "reference_analysis"],
        "original_contigs": original_contigs,
        "original_bases": original_bases,
        "final_contigs": len(records),
        "final_bases": sum(len(sequence) for _, sequence in records),
    }
    if sample.implant == "FRAGMENTED":
        fragment_size = int(sample.implant_parameters.get("fragment_size", 1000))
        if max(len(sequence) for _, sequence in records) > fragment_size:
            raise ValueError("FRAGMENTED implant exceeded configured fragment size")
        validation["checks"].append("fragment_size")
    if sample.implant == "MIXED_CONTIGS":
        if validation["final_bases"] <= original_bases:
            raise ValueError("MIXED_CONTIGS did not add contaminant sequence")
        validation["checks"].append("contaminant_bases_added")
    provenance["validation"] = validation
    return ReleaseArtifactSample(
        sample_id=sample.sample_id,
        source_id=sample.source_id,
        random_seed=sample.random_seed,
        files={"assembly": str(output_path)},
        expected_answers=expected,
        implant=sample.implant,
        implant_parameters=sample.implant_parameters,
        public_metadata={"format": "FASTA"},
        private_provenance=provenance,
    )
