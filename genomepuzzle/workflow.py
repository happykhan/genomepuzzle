"""Resumable release planning and SLURM execution."""

from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from genomepuzzle.contract import json_dump, sha256_file, utc_now
from genomepuzzle.provenance import PLANNED_GIT_COMMIT_ENV
from genomepuzzle.release import ReleaseSpec, load_release_spec, resolve_release_samples
from genomepuzzle.slurm import SlurmResources, build_stage_sbatch_script, submit_sbatch_script


PLAN_SCHEMA_VERSION = "1.0"
TERMINAL_SUCCESS = {"completed"}
RETRYABLE = {"planned", "failed", "cancelled"}
SCHEDULER_STATUS = {
    "PENDING": "submitted",
    "CONFIGURING": "submitted",
    "RUNNING": "running",
    "COMPLETING": "running",
    "COMPLETED": "completed",
    "CANCELLED": "cancelled",
    "FAILED": "failed",
    "TIMEOUT": "failed",
    "OUT_OF_MEMORY": "failed",
    "NODE_FAIL": "failed",
    "BOOT_FAIL": "failed",
    "PREEMPTED": "failed",
}


@dataclass(frozen=True)
class WorkflowStage:
    name: str
    command: tuple[str, ...]
    resources: SlurmResources
    dependencies: tuple[str, ...] = ()
    description: str = ""


@dataclass(frozen=True)
class WorkflowPlan:
    release_id: str
    spec_path: str
    output_dir: str
    repo_dir: str
    created_at: str
    git_commit: str
    pixi_lock_sha256: str
    samples: tuple[Mapping[str, Any], ...]
    stages: tuple[WorkflowStage, ...]
    schema_version: str = PLAN_SCHEMA_VERSION


def _git_commit(repo_dir: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_dir,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _resolve_input(spec_path: Path, raw: str) -> str:
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = spec_path.parent / path
    return str(path.resolve())


def _required_input(spec: ReleaseSpec, spec_path: Path, name: str) -> str:
    raw = spec.inputs.get(name)
    if not raw:
        raise ValueError(
            "release specification requires [inputs].{0} for {1}".format(
                name, spec.exercise
            )
        )
    return _resolve_input(spec_path, raw)


def build_release_plan(
    spec_file: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    *,
    repo_dir: str | os.PathLike[str] | None = None,
    id_salt: str | None = None,
    partition: str = "short",
) -> WorkflowPlan:
    """Resolve identities and create the deterministic release workflow."""

    spec_path = Path(spec_file).resolve()
    spec = load_release_spec(spec_path)
    output = Path(output_dir).resolve()
    repository = Path(repo_dir or Path(__file__).resolve().parents[1]).resolve()
    resolved = resolve_release_samples(spec, id_salt=id_salt)
    sample_plan = tuple(
        {
            "sample_id": sample.sample_id,
            "source_id": sample.source_id,
            "random_seed": sample.random_seed,
            "implant": sample.implant,
            "implant_parameters": dict(sample.implant_parameters),
        }
        for sample in resolved
    )

    base = [
        "pixi",
        "run",
        "genomepuzzle",
        "release",
    ]
    source_stage = None
    if spec.exercise == "typing":
        source_dir = _required_input(spec, spec_path, "source_dir")
        source_stage = WorkflowStage(
            name="sources",
            command=tuple(
                base
                + [
                    "fetch-assemblies",
                    "--spec",
                    str(spec_path),
                    "--output-dir",
                    source_dir,
                ]
            ),
            resources=SlurmResources(
                partition=partition,
                cpus=2,
                memory_gb=8,
                time_limit="02:00:00",
            ),
            description="Fetch and checksum every required target and fault source.",
        )
        build_command = base + [
            "build-typing",
            "--spec",
            str(spec_path),
            "--source-dir",
            source_dir,
            "--output-dir",
            str(output),
        ]
        resources = SlurmResources(partition=partition, cpus=4, memory_gb=16, time_limit="06:00:00")
    elif spec.exercise == "outbreak":
        build_command = base + ["generate-outbreak", "--spec", str(spec_path)]
        if spec.inputs.get("base_genome"):
            build_command.extend(
                [
                    "--base-genome",
                    _required_input(spec, spec_path, "base_genome"),
                ]
            )
            if spec.inputs.get("fault_genome"):
                build_command.extend(
                    [
                        "--fault-genome",
                        _required_input(spec, spec_path, "fault_genome"),
                    ]
                )
        else:
            build_command.extend(
                [
                    "--source-dir",
                    _required_input(spec, spec_path, "source_dir"),
                ]
            )
        build_command.extend(
            [
                "--metadata",
                _required_input(spec, spec_path, "metadata_csv"),
                "--output-dir",
                str(output),
            ]
        )
        resources = SlurmResources(partition=partition, cpus=4, memory_gb=16, time_limit="06:00:00")
    else:
        source_dir = _required_input(spec, spec_path, "source_dir")
        source_stage = WorkflowStage(
            name="sources",
            command=tuple(
                base
                + [
                    "fetch-assemblies",
                    "--spec",
                    str(spec_path),
                    "--output-dir",
                    source_dir,
                ]
            ),
            resources=SlurmResources(
                partition=partition,
                cpus=2,
                memory_gb=8,
                time_limit="02:00:00",
            ),
            description="Fetch and checksum every required target and fault source.",
        )
        build_command = base + [
            "generate-reads",
            "--spec",
            str(spec_path),
            "--source-dir",
            source_dir,
            "--output-dir",
            str(output),
        ]
        resources = SlurmResources(partition=partition, cpus=8, memory_gb=32, time_limit="12:00:00")

    validate_command = base + [
        "validate",
        "--release-dir",
        str(output),
        "--require-complete",
    ]
    generate_stage = WorkflowStage(
            name="generate",
            command=tuple(build_command),
            resources=resources,
            dependencies=("sources",) if source_stage else (),
            description="Generate, analyse, package and seal the release.",
        )
    stages = tuple(
        stage
        for stage in (
            source_stage,
            generate_stage,
            WorkflowStage(
                name="validate",
                command=tuple(validate_command),
                resources=SlurmResources(
                    partition=partition, cpus=1, memory_gb=4, time_limit="01:00:00"
                ),
                dependencies=("generate",),
                description="Independently validate the completed release contract.",
            ),
        )
        if stage is not None
    )
    lock_path = repository / "pixi.lock"
    return WorkflowPlan(
        release_id=spec.release_id,
        spec_path=str(spec_path),
        output_dir=str(output),
        repo_dir=str(repository),
        created_at=utc_now(),
        git_commit=_git_commit(repository),
        pixi_lock_sha256=sha256_file(lock_path) if lock_path.is_file() else "",
        samples=sample_plan,
        stages=stages,
    )


def _plan_payload(plan: WorkflowPlan) -> dict[str, Any]:
    payload = asdict(plan)
    payload["stages"] = [
        {
            **asdict(stage),
            "command": list(stage.command),
            "dependencies": list(stage.dependencies),
            "resources": asdict(stage.resources),
        }
        for stage in plan.stages
    ]
    return payload


def write_release_plan(plan: WorkflowPlan) -> Path:
    root = Path(plan.output_dir)
    build_dir = root / "build"
    stages_dir = build_dir / "stages"
    scripts_dir = build_dir / "scripts"
    logs_dir = build_dir / "logs"
    for path in (stages_dir, scripts_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)
    plan_path = build_dir / "plan.json"
    json_dump(plan_path, _plan_payload(plan))
    for stage in plan.stages:
        stage_path = stages_dir / "{0}.json".format(stage.name)
        if not stage_path.exists():
            json_dump(
                stage_path,
                {
                    "schema_version": PLAN_SCHEMA_VERSION,
                    "release_id": plan.release_id,
                    "stage": stage.name,
                    "description": stage.description,
                    "status": "planned",
                    "command": list(stage.command),
                    "command_display": shlex.join(stage.command),
                    "dependencies": list(stage.dependencies),
                    "resources": asdict(stage.resources),
                    "job_id": None,
                    "attempts": [],
                },
            )
        script = build_stage_sbatch_script(
            plan_path=plan_path,
            stage_name=stage.name,
            repo_dir=Path(plan.repo_dir),
            log_dir=logs_dir,
            resources=stage.resources,
            job_name="gp-{0}-{1}".format(plan.release_id[:28], stage.name),
        )
        script_path = scripts_dir / "{0}.sbatch".format(stage.name)
        script_path.write_text(script, encoding="utf-8")
    return plan_path


def load_plan(path: str | os.PathLike[str]) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("schema_version") != PLAN_SCHEMA_VERSION:
        raise ValueError("unsupported workflow plan schema")
    return payload


def _stage_path(plan: Mapping[str, Any], stage_name: str) -> Path:
    return Path(plan["output_dir"]) / "build" / "stages" / "{0}.json".format(stage_name)


def _read_stage(plan: Mapping[str, Any], stage_name: str) -> dict[str, Any]:
    path = _stage_path(plan, stage_name)
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def _write_stage(plan: Mapping[str, Any], stage: Mapping[str, Any]) -> None:
    json_dump(_stage_path(plan, str(stage["stage"])), dict(stage))


def _scheduler_state(job_id: str) -> str | None:
    """Return a normalised stage status from SLURM accounting when available."""

    try:
        result = subprocess.run(
            [
                "sacct",
                "-j",
                str(job_id),
                "--format=State",
                "-n",
                "-X",
                "-P",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    for line in result.stdout.splitlines():
        raw = line.strip().split()[0] if line.strip() else ""
        raw = raw.split("+", 1)[0]
        if raw in SCHEDULER_STATUS:
            return SCHEDULER_STATUS[raw]
    return None


def _reconcile_stage(plan: Mapping[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    job_id = stage.get("job_id")
    if not job_id or stage["status"] in TERMINAL_SUCCESS:
        return stage
    scheduler_status = _scheduler_state(str(job_id))
    if scheduler_status and scheduler_status != stage["status"]:
        stage["status"] = scheduler_status
        stage["scheduler_checked_at"] = utc_now()
        _write_stage(plan, stage)
    return stage


def submit_plan(
    plan_path: str | os.PathLike[str], *, retry: bool = False
) -> dict[str, str]:
    """Submit eligible stages with afterok dependencies and record job IDs."""

    plan = load_plan(plan_path)
    submitted: dict[str, str] = {}
    known_jobs: dict[str, str] = {}
    for stage_definition in plan["stages"]:
        name = stage_definition["name"]
        state = _reconcile_stage(plan, _read_stage(plan, name))
        if state["status"] in TERMINAL_SUCCESS:
            continue
        if state["status"] not in RETRYABLE:
            known = state.get("job_id")
            if known:
                known_jobs[name] = known
            continue
        dependency_jobs = []
        blocked = False
        for dependency in stage_definition["dependencies"]:
            dependency_state = _reconcile_stage(
                plan, _read_stage(plan, dependency)
            )
            if dependency_state["status"] in TERMINAL_SUCCESS:
                continue
            job_id = submitted.get(dependency) or known_jobs.get(dependency) or dependency_state.get("job_id")
            if not job_id:
                blocked = True
                break
            dependency_jobs.append(job_id)
        if blocked:
            continue
        script = (
            Path(plan["output_dir"]) / "build" / "scripts" / "{0}.sbatch".format(name)
        )
        job_id = submit_sbatch_script(script, dependency_job_ids=dependency_jobs)
        state["status"] = "submitted"
        state["job_id"] = job_id
        state["submitted_at"] = utc_now()
        state["attempts"].append(
            {"job_id": job_id, "submitted_at": state["submitted_at"]}
        )
        _write_stage(plan, state)
        submitted[name] = job_id
        known_jobs[name] = job_id
    return submitted


def run_stage(plan_path: str | os.PathLike[str], stage_name: str) -> None:
    """Execute one planned stage inside its SLURM allocation."""

    plan = load_plan(plan_path)
    definitions = {stage["name"]: stage for stage in plan["stages"]}
    if stage_name not in definitions:
        raise ValueError("unknown workflow stage: {0}".format(stage_name))
    definition = definitions[stage_name]
    state = _read_stage(plan, stage_name)
    prior_status = state["status"]
    if (
        stage_name == "generate"
        and len(state.get("attempts", [])) > 1
        and prior_status != "completed"
    ):
        root = Path(plan["output_dir"]).resolve()
        for name in ("public", "private"):
            partial_dir = root / name
            if partial_dir.is_dir():
                shutil.rmtree(partial_dir)
        for name in ("release.json", "COMPLETE", "COMPLETE.json"):
            partial_file = root / name
            if partial_file.is_file():
                partial_file.unlink()
    state["status"] = "running"
    state["started_at"] = utc_now()
    state["slurm"] = {
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "cpus": os.environ.get("SLURM_CPUS_PER_TASK"),
        "node": os.environ.get("SLURMD_NODENAME"),
    }
    _write_stage(plan, state)
    try:
        actual_commit = _git_commit(Path(plan["repo_dir"]))
        if actual_commit != plan["git_commit"]:
            raise RuntimeError(
                "workflow plan was created at Git commit {0}, but the repository "
                "is now at {1}; create a new plan before running this stage".format(
                    plan["git_commit"], actual_commit
                )
            )
        stage_environment = os.environ.copy()
        stage_environment[PLANNED_GIT_COMMIT_ENV] = plan["git_commit"]
        subprocess.run(
            definition["command"],
            cwd=plan["repo_dir"],
            env=stage_environment,
            check=True,
        )
    except BaseException as exc:
        state["status"] = "failed"
        state["finished_at"] = utc_now()
        state["error"] = "{0}: {1}".format(type(exc).__name__, exc)
        _write_stage(plan, state)
        raise
    state["status"] = "completed"
    state["finished_at"] = utc_now()
    state.pop("error", None)
    _write_stage(plan, state)


def workflow_status(plan_path: str | os.PathLike[str]) -> list[dict[str, Any]]:
    plan = load_plan(plan_path)
    return [
        _reconcile_stage(plan, _read_stage(plan, stage["name"]))
        for stage in plan["stages"]
    ]
