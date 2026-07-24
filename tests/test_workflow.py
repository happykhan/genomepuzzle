import json
from pathlib import Path

from genomepuzzle.workflow import (
    build_release_plan,
    run_stage,
    submit_plan,
    workflow_status,
    write_release_plan,
)


def _spec(path: Path, source_dir: Path):
    path.write_text(
        """
release_id = "workflow-test"
exercise = "assembly"
mode = "practice"

[inputs]
source_dir = "{source_dir}"

[[samples]]
source_id = "source-a"
public_id = "Sample_a"
[samples.expected_answers]
species = "Klebsiella pneumoniae"
qc = "pass"
error = "none"
""".format(source_dir=source_dir).strip()
        + "\n",
        encoding="utf-8",
    )


def test_workflow_plan_writes_resumable_slurm_state(tmp_path):
    source = tmp_path / "sources"
    source.mkdir()
    spec = tmp_path / "release.toml"
    _spec(spec, source)
    output = tmp_path / "output"

    plan = build_release_plan(spec, output, repo_dir=Path(__file__).parents[1])
    plan_path = write_release_plan(plan)
    payload = json.loads(plan_path.read_text())
    script = (output / "build/scripts/generate.sbatch").read_text()

    assert [stage["name"] for stage in payload["stages"]] == [
        "sources",
        "generate",
        "validate",
    ]
    assert payload["stages"][1]["dependencies"] == ["sources"]
    assert "pixi run genomepuzzle release run-stage" in script
    assert workflow_status(plan_path)[0]["status"] == "planned"


def test_submit_plan_records_afterok_dependencies(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    spec = tmp_path / "release.toml"
    _spec(spec, source)
    output = tmp_path / "output"
    plan_path = write_release_plan(
        build_release_plan(spec, output, repo_dir=Path(__file__).parents[1])
    )
    calls = []

    def fake_submit(path, dependency_job_ids=()):
        calls.append((Path(path).name, tuple(dependency_job_ids)))
        return str(100 + len(calls))

    monkeypatch.setattr("genomepuzzle.workflow.submit_sbatch_script", fake_submit)
    monkeypatch.setattr(
        "genomepuzzle.workflow._scheduler_state", lambda job_id: "submitted"
    )
    jobs = submit_plan(plan_path)

    assert jobs == {"sources": "101", "generate": "102", "validate": "103"}
    assert calls == [
        ("sources.sbatch", ()),
        ("generate.sbatch", ("101",)),
        ("validate.sbatch", ("102",)),
    ]


def test_resume_does_not_duplicate_a_running_job(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    spec = tmp_path / "release.toml"
    _spec(spec, source)
    output = tmp_path / "output"
    plan_path = write_release_plan(
        build_release_plan(spec, output, repo_dir=Path(__file__).parents[1])
    )
    sources_path = output / "build/stages/sources.json"
    sources = json.loads(sources_path.read_text())
    sources.update({"status": "running", "job_id": "122"})
    sources_path.write_text(json.dumps(sources))
    stage_path = output / "build/stages/generate.json"
    stage = json.loads(stage_path.read_text())
    stage.update({"status": "running", "job_id": "123"})
    stage_path.write_text(json.dumps(stage))
    validation_path = output / "build/stages/validate.json"
    validation = json.loads(validation_path.read_text())
    validation.update({"status": "submitted", "job_id": "124"})
    validation_path.write_text(json.dumps(validation))
    monkeypatch.setattr(
        "genomepuzzle.workflow._scheduler_state", lambda job_id: "running"
    )
    monkeypatch.setattr(
        "genomepuzzle.workflow.submit_sbatch_script",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("duplicate")),
    )

    assert submit_plan(plan_path, retry=True) == {}


def test_status_reconciles_cancelled_slurm_job(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    spec = tmp_path / "release.toml"
    _spec(spec, source)
    output = tmp_path / "output"
    plan_path = write_release_plan(
        build_release_plan(spec, output, repo_dir=Path(__file__).parents[1])
    )
    stage_path = output / "build/stages/generate.json"
    stage = json.loads(stage_path.read_text())
    stage.update({"status": "running", "job_id": "123"})
    stage_path.write_text(json.dumps(stage))
    monkeypatch.setattr(
        "genomepuzzle.workflow._scheduler_state", lambda job_id: "cancelled"
    )

    statuses = {row["stage"]: row["status"] for row in workflow_status(plan_path)}
    assert statuses["generate"] == "cancelled"


def test_outbreak_plan_uses_native_generator_when_base_genome_is_declared(tmp_path):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">r\n" + "A" * 1000 + "\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "Sample,Cluster,SPECIES\none,A,Klebsiella pneumoniae\n",
        encoding="utf-8",
    )
    spec = tmp_path / "outbreak.toml"
    spec.write_text(
        """
release_id = "outbreak-plan"
exercise = "outbreak"
mode = "practice"
[inputs]
base_genome = "{reference}"
metadata_csv = "{metadata}"
[[samples]]
source_id = "one"
public_id = "Sample_one"
""".format(reference=reference, metadata=metadata).strip()
        + "\n",
        encoding="utf-8",
    )
    plan = build_release_plan(
        spec, tmp_path / "output", repo_dir=Path(__file__).parents[1]
    )
    assert "generate-outbreak" in plan.stages[0].command
    assert "--base-genome" in plan.stages[0].command


def test_run_stage_freezes_planned_git_commit_in_environment(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    spec = tmp_path / "release.toml"
    _spec(spec, source)
    output = tmp_path / "output"
    plan_path = write_release_plan(
        build_release_plan(spec, output, repo_dir=Path(__file__).parents[1])
    )
    plan = json.loads(plan_path.read_text())
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))

    monkeypatch.setattr("genomepuzzle.workflow.subprocess.run", fake_run)
    run_stage(plan_path, "generate")

    assert calls[0][1]["env"]["GENOMEPUZZLE_PLANNED_GIT_COMMIT"] == plan["git_commit"]
