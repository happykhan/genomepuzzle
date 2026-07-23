"""Spec-driven short-read and hybrid dataset generation."""

from __future__ import annotations

import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from genomepuzzle.create_error import (
    concatenate_fastqs,
    degrade_quality,
    subsample_paired_fastq,
    subsample_single_fastq,
    truncate_fastq,
)
from genomepuzzle.reads_release import package_read_release
from genomepuzzle.release import ReleaseSpec, resolve_release_samples
from genomepuzzle.runtime import require_tool


ASSEMBLY_IMPLANTS = {
    "NORMAL",
    "NONE",
    "LOW_COVERAGE",
    "POOR_QUALITY",
    "TRUNCATED",
    "CONTAMINATED",
}
HYBRID_IMPLANTS = {
    "NORMAL",
    "NONE",
    "LOW_SHORT_COVERAGE",
    "LOW_LONG_COVERAGE",
    "LONG_READ_QUALITY",
    "CONTAMINATED",
}


def _compression_threads() -> str:
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    return str(max(1, allocated // 4))


def _source_fasta(source_dir: Path, source_id: str) -> Path:
    for suffix in ("", ".fasta", ".fna", ".fa"):
        candidate = source_dir / "{0}{1}".format(source_id, suffix)
        if candidate.is_file():
            return candidate
    raise ValueError("source assembly not found for {0}".format(source_id))


def _simulate_short_reads(
    reference: Path,
    output_r1: Path,
    output_r2: Path,
    *,
    seed: int,
    coverage: float,
    read_length: int,
    fragment_length: int,
    fragment_sd: int,
) -> None:
    art = require_tool("art_illumina")
    prefix = output_r1.parent / ".art-{0}-".format(output_r1.stem)
    command = [
        art,
        "-ss",
        "HS25",
        "-i",
        str(reference),
        "-l",
        str(read_length),
        "-f",
        str(coverage),
        "-o",
        str(prefix),
        "-p",
        "-m",
        str(fragment_length),
        "-s",
        str(fragment_sd),
        "--rndSeed",
        str(seed),
        "-na",
    ]
    subprocess.run(command, check=True)
    generated_r1 = Path(str(prefix) + "1.fq")
    generated_r2 = Path(str(prefix) + "2.fq")
    pigz = require_tool("pigz")
    subprocess.run(
        [
            pigz,
            "-n",
            "-p",
            _compression_threads(),
            "-f",
            str(generated_r1),
            str(generated_r2),
        ],
        check=True,
    )
    os.replace(str(generated_r1) + ".gz", output_r1)
    os.replace(str(generated_r2) + ".gz", output_r2)


def _simulate_long_reads(
    reference: Path, output: Path, *, seed: int, quantity: str
) -> None:
    badread = require_tool("badread")
    command = [
        badread,
        "simulate",
        "--reference",
        str(reference),
        "--quantity",
        quantity,
        "--seed",
        str(seed),
    ]
    pigz = require_tool("pigz")
    progress_log = output.with_name(".{0}.badread.log".format(output.name))
    with open(output, "wb") as handle, open(
        progress_log, "w", encoding="utf-8"
    ) as progress:
        simulator = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=progress,
        )
        compressor = subprocess.Popen(
            [pigz, "-n", "-p", _compression_threads(), "-c"],
            stdin=simulator.stdout,
            stdout=handle,
        )
        if simulator.stdout:
            simulator.stdout.close()
        compressor_code = compressor.wait()
        simulator_code = simulator.wait()
    if simulator_code or compressor_code:
        raise subprocess.CalledProcessError(
            simulator_code or compressor_code,
            command,
            stderr="Badread progress retained in {0}".format(progress_log),
        )
    progress_log.unlink()


def _parameter(sample, name: str, default):
    return sample.implant_parameters.get(name, default)


def _copy_pair(source_r1: Path, source_r2: Path, target_r1: Path, target_r2: Path) -> None:
    shutil.copyfile(source_r1, target_r1)
    shutil.copyfile(source_r2, target_r2)


def _generate_contaminant(
    source_dir: Path,
    work_dir: Path,
    source_id: str,
    *,
    seed: int,
    short_parameters: dict,
    include_long: bool,
    long_quantity: str,
) -> tuple[Path, Path, Path | None]:
    reference = _source_fasta(source_dir, source_id)
    prefix = work_dir / "contaminant-{0}".format(seed)
    r1 = Path(str(prefix) + "_R1.fastq.gz")
    r2 = Path(str(prefix) + "_R2.fastq.gz")
    _simulate_short_reads(reference, r1, r2, seed=seed, **short_parameters)
    long_reads = None
    if include_long:
        long_reads = Path(str(prefix) + "_long.fastq.gz")
        _simulate_long_reads(
            reference, long_reads, seed=seed + 1, quantity=long_quantity
        )
    return r1, r2, long_reads


def _generate_sample(
    spec: ReleaseSpec,
    sample,
    source_path: Path,
    work_dir: Path,
) -> tuple[str, dict[str, object], dict[str, object]]:
    allowed = ASSEMBLY_IMPLANTS if spec.exercise == "assembly" else HYBRID_IMPLANTS
    if sample.implant not in allowed:
        raise ValueError(
            "unsupported {0} implant: {1}".format(spec.exercise, sample.implant)
        )
    reference = _source_fasta(source_path, sample.source_id)
    short_parameters = {
        "coverage": float(_parameter(sample, "short_coverage", 40)),
        "read_length": int(_parameter(sample, "read_length", 150)),
        "fragment_length": int(_parameter(sample, "fragment_length", 300)),
        "fragment_sd": int(_parameter(sample, "fragment_sd", 50)),
    }
    base_r1 = work_dir / ".{0}_base_R1.fastq.gz".format(sample.source_id)
    base_r2 = work_dir / ".{0}_base_R2.fastq.gz".format(sample.source_id)
    final_r1 = work_dir / "{0}_R1.fastq.gz".format(sample.source_id)
    final_r2 = work_dir / "{0}_R2.fastq.gz".format(sample.source_id)
    _simulate_short_reads(
        reference,
        base_r1,
        base_r2,
        seed=sample.random_seed,
        **short_parameters,
    )
    contaminant_long_for_hybrid: Path | None = None
    achieved: dict[str, object] = {
        "status": "passed",
        "implant": sample.implant,
        "checks": ["read_simulation", "implant_materialized"],
    }

    if sample.implant in {"LOW_COVERAGE", "LOW_SHORT_COVERAGE"}:
        fraction = float(_parameter(sample, "read_fraction", 0.15))
        subsample_paired_fastq(
            str(base_r1),
            str(base_r2),
            str(final_r1),
            str(final_r2),
            fraction,
            sample.random_seed,
        )
        achieved["read_fraction"] = fraction
    elif sample.implant == "POOR_QUALITY":
        minimum = int(_parameter(sample, "min_quality", 5))
        maximum = int(_parameter(sample, "max_quality", 14))
        degrade_quality(
            str(base_r1), str(final_r1), minimum, maximum, sample.random_seed
        )
        degrade_quality(
            str(base_r2), str(final_r2), minimum, maximum, sample.random_seed + 1
        )
        achieved["quality_range"] = [minimum, maximum]
    elif sample.implant == "TRUNCATED":
        length = int(_parameter(sample, "read_length", 35))
        truncate_fastq(str(base_r1), str(final_r1), length)
        truncate_fastq(str(base_r2), str(final_r2), length)
        achieved["read_length"] = length
    elif sample.implant == "CONTAMINATED":
        contaminant_id = sample.implant_parameters.get("contaminant_source_id")
        if not isinstance(contaminant_id, str) or not contaminant_id:
            raise ValueError("CONTAMINATED requires contaminant_source_id")
        (
            contaminant_r1,
            contaminant_r2,
            contaminant_long_for_hybrid,
        ) = _generate_contaminant(
            source_path,
            work_dir,
            contaminant_id,
            seed=sample.random_seed + 101,
            short_parameters=short_parameters,
            include_long=spec.exercise == "hybrid",
            long_quantity=str(_parameter(sample, "long_quantity", "30x")),
        )
        fraction = float(_parameter(sample, "contamination_fraction", 0.2))
        dirty_r1 = work_dir / ".{0}_dirty_R1.fastq.gz".format(sample.source_id)
        dirty_r2 = work_dir / ".{0}_dirty_R2.fastq.gz".format(sample.source_id)
        subsample_paired_fastq(
            str(contaminant_r1),
            str(contaminant_r2),
            str(dirty_r1),
            str(dirty_r2),
            fraction,
            sample.random_seed,
        )
        concatenate_fastqs([str(base_r1), str(dirty_r1)], str(final_r1))
        concatenate_fastqs([str(base_r2), str(dirty_r2)], str(final_r2))
        achieved["contamination_fraction"] = fraction
        achieved["contaminant_source_id"] = contaminant_id
    else:
        _copy_pair(base_r1, base_r2, final_r1, final_r2)

    if spec.exercise == "hybrid":
        base_long = work_dir / ".{0}_base_long.fastq.gz".format(sample.source_id)
        final_long = work_dir / "{0}_long.fastq.gz".format(sample.source_id)
        quantity = str(_parameter(sample, "long_quantity", "30x"))
        _simulate_long_reads(
            reference, base_long, seed=sample.random_seed + 1, quantity=quantity
        )
        if sample.implant == "LOW_LONG_COVERAGE":
            fraction = float(_parameter(sample, "read_fraction", 0.15))
            subsample_single_fastq(
                str(base_long),
                str(final_long),
                fraction,
                sample.random_seed,
            )
            achieved["long_read_fraction"] = fraction
        elif sample.implant == "LONG_READ_QUALITY":
            minimum = int(_parameter(sample, "min_quality", 5))
            maximum = int(_parameter(sample, "max_quality", 14))
            degrade_quality(
                str(base_long),
                str(final_long),
                minimum,
                maximum,
                sample.random_seed,
            )
            achieved["long_quality_range"] = [minimum, maximum]
        elif sample.implant == "CONTAMINATED":
            if contaminant_long_for_hybrid is None:
                raise RuntimeError("hybrid contaminant long reads were not generated")
            dirty_long = work_dir / ".{0}_dirty_long.fastq.gz".format(
                sample.source_id
            )
            fraction = float(_parameter(sample, "contamination_fraction", 0.2))
            subsample_single_fastq(
                str(contaminant_long_for_hybrid),
                str(dirty_long),
                fraction,
                sample.random_seed,
            )
            concatenate_fastqs([str(base_long), str(dirty_long)], str(final_long))
        else:
            shutil.copyfile(base_long, final_long)

    answers = dict(sample.expected_answers)
    answers.setdefault(
        "qc", "pass" if sample.implant in {"NORMAL", "NONE"} else "fail"
    )
    answers.setdefault(
        "error", "none" if sample.implant in {"NORMAL", "NONE"} else sample.implant
    )
    if "species" not in answers:
        species = sample.implant_parameters.get("species")
        if not species:
            raise ValueError(
                "{0} requires expected_answers.species or "
                "implant_parameters.species".format(sample.source_id)
            )
        answers["species"] = species
    return sample.source_id, answers, achieved


def generate_read_release(
    spec: ReleaseSpec,
    source_dir: str | os.PathLike[str],
    release_dir: str | os.PathLike[str],
    *,
    id_salt: str | None = None,
) -> dict[str, str]:
    """Generate final participant reads and immediately package the release."""

    if spec.exercise not in {"assembly", "hybrid"}:
        raise ValueError("read generation requires assembly or hybrid exercise")
    source_path = Path(source_dir).resolve()
    destination = Path(release_dir).resolve()
    work_dir = destination / "build" / "work" / "generated-reads"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)
    samples = list(resolve_release_samples(spec, id_salt=id_salt))
    source_ids = [sample.source_id for sample in samples]
    if len(source_ids) != len(set(source_ids)):
        raise ValueError(
            "assembly and hybrid releases require unique sample source_id values"
        )
    allocated_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    workers = min(4, max(1, allocated_cpus), len(samples))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(
            executor.map(
                lambda sample: _generate_sample(
                    spec, sample, source_path, work_dir
                ),
                samples,
            )
        )
    expected = {source_id: answers for source_id, answers, _ in results}
    validations = {
        source_id: validation for source_id, _, validation in results
    }

    return package_read_release(
        spec,
        work_dir,
        expected,
        destination,
        id_salt=id_salt,
        implant_validations=validations,
    )
