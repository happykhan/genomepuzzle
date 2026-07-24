"""Spec-driven short-read and hybrid dataset generation."""

from __future__ import annotations

import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from genomepuzzle.contract import failure_reason_for_implant
from genomepuzzle.create_error import (
    concatenate_fastqs,
    count_reads,
    subsample_paired_fastq,
    subsample_paired_read_by_count,
    subsample_single_fastq,
    subsample_single_fastq_by_count,
)
from genomepuzzle.reads_release import package_read_release
from genomepuzzle.release import ReleaseSpec, resolve_release_samples
from genomepuzzle.runtime import require_tool
from genomepuzzle.typing import read_fasta, write_anonymous_fasta


ASSEMBLY_IMPLANTS = {
    "NORMAL",
    "NONE",
    "LOW_COVERAGE",
    "CONTAMINATED",
    "ZERO_BYTE_R1",
    "ZERO_BYTE_R2",
    "MISSING_R1",
    "MISSING_R2",
    "TEN_READ_PAIRS",
    "TRUNCATED_R1_TO_10_READS",
    "TRUNCATED_R2_TO_10_READS",
    "WRONG_ORGANISM",
}
HYBRID_IMPLANTS = {
    "NORMAL",
    "NONE",
    "LOW_SHORT_COVERAGE",
    "CONTAMINATED",
    "ZERO_BYTE_R1",
    "ZERO_BYTE_R2",
    "MISSING_R1",
    "MISSING_R2",
    "TEN_READ_PAIRS",
    "TRUNCATED_R1_TO_10_READS",
    "TRUNCATED_R2_TO_10_READS",
    "MISSING_LONG_READS",
    "ZERO_BYTE_LONG_READS",
    "TEN_LONG_READS",
    "WRONG_ORGANISM",
    "DISCORDANT_READ_SETS",
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


def _simulation_reference(
    source: Path, work_dir: Path, filename: str, sample_id: str
) -> Path:
    """Create a reference whose contig names are safe for simulator headers."""

    target = work_dir / filename
    write_anonymous_fasta(target, read_fasta(source), sample_id)
    return target


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
        "--junk_reads",
        "0",
        "--random_reads",
        "0",
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
    sample_id: str,
) -> tuple[Path, Path, Path | None]:
    source_reference = _source_fasta(source_dir, source_id)
    prefix = work_dir / "contaminant-{0}".format(seed)
    reference = _simulation_reference(
        source_reference,
        work_dir,
        ".contaminant-{0}.reference.fasta".format(seed),
        sample_id,
    )
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
    source_reference = _source_fasta(source_path, sample.source_id)
    reference = _simulation_reference(
        source_reference,
        work_dir,
        ".{0}.reference.fasta".format(sample.source_id),
        sample.sample_id,
    )
    short_parameters = {
        "coverage": float(_parameter(sample, "short_coverage", 30)),
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
    contaminant_r1: Path | None = None
    contaminant_r2: Path | None = None
    contaminant_long_for_hybrid: Path | None = None
    achieved: dict[str, object] = {
        "status": "passed",
        "fault_type": sample.implant,
        "checks": ["read_simulation", "implant_materialized"],
    }
    if sample.implant in {
        "CONTAMINATED",
        "WRONG_ORGANISM",
        "DISCORDANT_READ_SETS",
    }:
        source_parameter = (
            "replacement_source_id"
            if sample.implant == "WRONG_ORGANISM"
            else "contaminant_source_id"
        )
        contaminant_id = sample.implant_parameters.get(source_parameter)
        if not isinstance(contaminant_id, str) or not contaminant_id:
            raise ValueError(
                "{0} requires {1}".format(sample.implant, source_parameter)
            )
        if contaminant_id == sample.source_id:
            raise ValueError(
                "{0} fault source must differ from target".format(sample.implant)
            )
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
            long_quantity=str(_parameter(sample, "long_quantity", "10x")),
            sample_id=sample.sample_id,
        )
        achieved[source_parameter] = contaminant_id

    if sample.implant in {"LOW_COVERAGE", "LOW_SHORT_COVERAGE"}:
        fraction = float(_parameter(sample, "read_fraction", 0.02))
        achieved_coverage = short_parameters["coverage"] * fraction
        if achieved_coverage > 1.0:
            raise ValueError(
                "LOW_COVERAGE must produce no more than 1x expected coverage"
            )
        subsample_paired_fastq(
            str(base_r1),
            str(base_r2),
            str(final_r1),
            str(final_r2),
            fraction,
            sample.random_seed,
        )
        achieved["read_fraction"] = fraction
        achieved["expected_short_coverage"] = round(achieved_coverage, 4)
    elif sample.implant == "TEN_READ_PAIRS":
        subsample_paired_read_by_count(
            str(base_r1),
            str(base_r2),
            str(final_r1),
            str(final_r2),
            num_reads=10,
            random_seed=sample.random_seed,
        )
        achieved["short_read_pairs"] = 10
    elif sample.implant in {
        "TRUNCATED_R1_TO_10_READS",
        "TRUNCATED_R2_TO_10_READS",
    }:
        if sample.implant == "TRUNCATED_R1_TO_10_READS":
            subsample_single_fastq_by_count(
                str(base_r1),
                str(final_r1),
                num_reads=10,
                random_seed=sample.random_seed,
            )
            shutil.copyfile(base_r2, final_r2)
            truncated_role = "read_1"
        else:
            shutil.copyfile(base_r1, final_r1)
            subsample_single_fastq_by_count(
                str(base_r2),
                str(final_r2),
                num_reads=10,
                random_seed=sample.random_seed,
            )
            truncated_role = "read_2"
        achieved["truncated_role"] = truncated_role
        achieved["truncated_role_reads"] = 10
    elif sample.implant in {"ZERO_BYTE_R1", "ZERO_BYTE_R2"}:
        _copy_pair(base_r1, base_r2, final_r1, final_r2)
        empty_path = final_r1 if sample.implant == "ZERO_BYTE_R1" else final_r2
        empty_path.write_bytes(b"")
        achieved["zero_byte_role"] = (
            "read_1" if sample.implant == "ZERO_BYTE_R1" else "read_2"
        )
        achieved["zero_byte_size"] = empty_path.stat().st_size
    elif sample.implant in {"MISSING_R1", "MISSING_R2"}:
        _copy_pair(base_r1, base_r2, final_r1, final_r2)
        missing_path = final_r1 if sample.implant == "MISSING_R1" else final_r2
        missing_path.unlink()
        achieved["missing_role"] = (
            "read_1" if sample.implant == "MISSING_R1" else "read_2"
        )
    elif sample.implant == "CONTAMINATED":
        assert contaminant_r1 is not None and contaminant_r2 is not None
        fraction = float(_parameter(sample, "contamination_fraction", 0.5))
        if not 0.3 <= fraction <= 0.9:
            raise ValueError(
                "contamination_fraction must be between 0.30 and 0.90"
            )
        clean_r1 = work_dir / ".{0}_clean_R1.fastq.gz".format(sample.source_id)
        clean_r2 = work_dir / ".{0}_clean_R2.fastq.gz".format(sample.source_id)
        dirty_r1 = work_dir / ".{0}_dirty_R1.fastq.gz".format(sample.source_id)
        dirty_r2 = work_dir / ".{0}_dirty_R2.fastq.gz".format(sample.source_id)
        base_pairs = count_reads(str(base_r1))
        available_contaminant_pairs = count_reads(str(contaminant_r1))
        total_pairs = int(
            min(
                base_pairs / (1 - fraction),
                available_contaminant_pairs / fraction,
            )
        )
        requested_contaminant_pairs = round(total_pairs * fraction)
        requested_clean_pairs = total_pairs - requested_contaminant_pairs
        subsample_paired_read_by_count(
            str(base_r1),
            str(base_r2),
            str(clean_r1),
            str(clean_r2),
            num_reads=requested_clean_pairs,
            random_seed=sample.random_seed,
        )
        subsample_paired_read_by_count(
            str(contaminant_r1),
            str(contaminant_r2),
            str(dirty_r1),
            str(dirty_r2),
            num_reads=requested_contaminant_pairs,
            random_seed=sample.random_seed + 1,
        )
        concatenate_fastqs([str(clean_r1), str(dirty_r1)], str(final_r1))
        concatenate_fastqs([str(clean_r2), str(dirty_r2)], str(final_r2))
        clean_pairs = count_reads(str(clean_r1))
        contaminant_pairs = count_reads(str(dirty_r1))
        achieved_fraction = contaminant_pairs / (clean_pairs + contaminant_pairs)
        if abs(achieved_fraction - fraction) > 1 / total_pairs:
            raise ValueError(
                "achieved contamination fraction differs from requested fraction"
            )
        achieved["requested_contamination_fraction"] = fraction
        achieved["achieved_contamination_fraction"] = round(
            achieved_fraction, 6
        )
        achieved["clean_pairs"] = clean_pairs
        achieved["contaminant_pairs"] = contaminant_pairs
    elif sample.implant == "WRONG_ORGANISM":
        assert contaminant_r1 is not None and contaminant_r2 is not None
        _copy_pair(contaminant_r1, contaminant_r2, final_r1, final_r2)
        achieved["replacement_fraction"] = 1.0
    else:
        _copy_pair(base_r1, base_r2, final_r1, final_r2)

    if spec.exercise == "hybrid":
        base_long = work_dir / ".{0}_base_long.fastq.gz".format(sample.source_id)
        final_long = work_dir / "{0}_long.fastq.gz".format(sample.source_id)
        quantity = str(_parameter(sample, "long_quantity", "10x"))
        _simulate_long_reads(
            reference, base_long, seed=sample.random_seed + 1, quantity=quantity
        )
        if sample.implant == "TEN_LONG_READS":
            subsample_single_fastq_by_count(
                str(base_long),
                str(final_long),
                num_reads=10,
                random_seed=sample.random_seed,
            )
            achieved["long_read_count"] = 10
        elif sample.implant == "ZERO_BYTE_LONG_READS":
            final_long.write_bytes(b"")
            achieved["zero_byte_role"] = "long_reads"
            achieved["zero_byte_size"] = final_long.stat().st_size
        elif sample.implant == "MISSING_LONG_READS":
            achieved["missing_role"] = "long_reads"
        elif sample.implant == "CONTAMINATED":
            if contaminant_long_for_hybrid is None:
                raise RuntimeError("hybrid contaminant long reads were not generated")
            clean_long = work_dir / ".{0}_clean_long.fastq.gz".format(
                sample.source_id
            )
            dirty_long = work_dir / ".{0}_dirty_long.fastq.gz".format(
                sample.source_id
            )
            fraction = float(_parameter(sample, "contamination_fraction", 0.5))
            base_long_count = count_reads(str(base_long))
            available_contaminant_long = count_reads(
                str(contaminant_long_for_hybrid)
            )
            total_long = int(
                min(
                    base_long_count / (1 - fraction),
                    available_contaminant_long / fraction,
                )
            )
            requested_contaminant_long = round(total_long * fraction)
            requested_clean_long = total_long - requested_contaminant_long
            subsample_single_fastq_by_count(
                str(base_long),
                str(clean_long),
                num_reads=requested_clean_long,
                random_seed=sample.random_seed,
            )
            subsample_single_fastq_by_count(
                str(contaminant_long_for_hybrid),
                str(dirty_long),
                num_reads=requested_contaminant_long,
                random_seed=sample.random_seed + 1,
            )
            concatenate_fastqs([str(clean_long), str(dirty_long)], str(final_long))
            clean_long_count = count_reads(str(clean_long))
            contaminant_long_count = count_reads(str(dirty_long))
            long_fraction = contaminant_long_count / (
                clean_long_count + contaminant_long_count
            )
            if abs(long_fraction - fraction) > 1 / total_long:
                raise ValueError(
                    "achieved long-read contamination differs from requested fraction"
                )
            achieved["achieved_long_contamination_fraction"] = round(
                long_fraction, 6
            )
        elif sample.implant in {"WRONG_ORGANISM", "DISCORDANT_READ_SETS"}:
            if contaminant_long_for_hybrid is None:
                raise RuntimeError("hybrid contaminant long reads were not generated")
            shutil.copyfile(contaminant_long_for_hybrid, final_long)
            achieved["long_read_source"] = "contaminant"
        else:
            shutil.copyfile(base_long, final_long)

    answers = dict(sample.expected_answers)
    answers.setdefault(
        "qc_status", "PASS" if sample.implant in {"NORMAL", "NONE"} else "FAIL"
    )
    answers.setdefault(
        "failure_reason",
        failure_reason_for_implant(spec.exercise, sample.implant),
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
        preanonymized=True,
    )
