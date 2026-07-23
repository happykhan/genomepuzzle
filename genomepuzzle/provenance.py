"""Reproducibility metadata captured for every generated release."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

from genomepuzzle.contract import sha256_file, utc_now


TRACKED_PACKAGES = {
    "art",
    "badread",
    "flye",
    "iqtree",
    "kleborate",
    "mashtree",
    "minimap2",
    "ncbi-datasets-cli",
    "pigz",
    "python",
    "seqtk",
    "spades",
    "sra-tools",
}


def _git_commit(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _pixi_packages() -> dict[str, dict[str, str]]:
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        return {}
    metadata_dir = Path(prefix) / "conda-meta"
    packages: dict[str, dict[str, str]] = {}
    for path in metadata_dir.glob("*.json"):
        try:
            with open(path, encoding="utf-8") as handle:
                payload = json.load(handle)
        except (OSError, json.JSONDecodeError):
            continue
        name = payload.get("name")
        if name in TRACKED_PACKAGES:
            packages[name] = {
                "version": str(payload.get("version", "")),
                "build": str(payload.get("build", "")),
                "channel": str(payload.get("channel", "")),
            }
    return dict(sorted(packages.items()))


def runtime_provenance() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    lock = root / "pixi.lock"
    return {
        "generated_at": utc_now(),
        "genomepuzzle_git_commit": _git_commit(root),
        "pixi_lock_sha256": sha256_file(lock) if lock.is_file() else None,
        "pixi_packages": _pixi_packages(),
        "slurm": {
            "job_id": os.environ.get("SLURM_JOB_ID"),
            "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
            "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
            "job_partition": os.environ.get("SLURM_JOB_PARTITION"),
            "node": os.environ.get("SLURMD_NODENAME"),
        },
    }
