"""SLURM helpers for cluster-backed dataset generation."""

import os
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class SlurmResources:
    partition: str = "short"
    cpus: int = 4
    memory_gb: int = 16
    time_limit: str = "06:00:00"


def build_stage_sbatch_script(
    *,
    plan_path: Path,
    stage_name: str,
    repo_dir: Path,
    log_dir: Path,
    resources: SlurmResources,
    job_name: str,
) -> str:
    command = shlex.join(
        [
            "pixi",
            "run",
            "genomepuzzle",
            "release",
            "run-stage",
            "--plan",
            str(plan_path),
            "--stage",
            stage_name,
        ]
    )
    return """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -p {partition}
#SBATCH -c {cpus}
#SBATCH --mem={memory_gb}G
#SBATCH -t {time_limit}
#SBATCH -o {stdout}
#SBATCH -e {stderr}

set -euo pipefail
cd {repo_dir}
{command}
""".format(
        job_name=job_name,
        partition=resources.partition,
        cpus=resources.cpus,
        memory_gb=resources.memory_gb,
        time_limit=resources.time_limit,
        stdout=shlex.quote(str(log_dir / "%x-%j.out")),
        stderr=shlex.quote(str(log_dir / "%x-%j.err")),
        repo_dir=shlex.quote(str(repo_dir)),
        command=command,
    )


def build_hybrid_sbatch_script(
    samplelist,
    output_dir,
    mode="challenge",
    contamination_list=None,
    random_seed=42,
    job_name="genomepuzzle-hybrid",
    partition="short",
    cpus_per_task=8,
    mem_gb=32,
    time_limit="12:00:00",
    repo_dir=None,
    log_dir=None,
):
    repo_dir = repo_dir or os.getcwd()
    log_dir = log_dir or os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    command = [
        "pixi",
        "run",
        "genomepuzzle",
        "long",
        "hybrid",
        "--samplelist",
        samplelist,
        "--mode",
        mode,
        "--output-dir",
        output_dir,
        "--random-seed",
        str(random_seed),
    ]
    if contamination_list:
        command.extend(["--contamination-list", contamination_list])
    wrapped_command = " ".join(subprocess.list2cmdline([part]) for part in command)
    return """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -p {partition}
#SBATCH -c {cpus}
#SBATCH --mem={mem}G
#SBATCH -t {time_limit}
#SBATCH -o {log_dir}/%x-%j.out
#SBATCH -e {log_dir}/%x-%j.err

set -euo pipefail
cd {repo_dir}
mkdir -p {output_dir}
{wrapped_command}
""".format(
        job_name=job_name,
        partition=partition,
        cpus=cpus_per_task,
        mem=mem_gb,
        time_limit=time_limit,
        log_dir=log_dir,
        repo_dir=repo_dir,
        output_dir=output_dir,
        wrapped_command=wrapped_command,
    )


def submit_sbatch_script(
    script_path: str | os.PathLike[str],
    dependency_job_ids: Sequence[str] = (),
):
    command = ["sbatch", "--parsable"]
    if dependency_job_ids:
        command.append(
            "--dependency=afterok:{0}".format(":".join(dependency_job_ids))
        )
    command.append(str(script_path))
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()
