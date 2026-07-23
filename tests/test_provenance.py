from genomepuzzle.provenance import runtime_provenance


def test_runtime_provenance_prefers_frozen_workflow_commit(monkeypatch):
    monkeypatch.setenv("GENOMEPUZZLE_PLANNED_GIT_COMMIT", "frozen-commit")
    monkeypatch.setattr(
        "genomepuzzle.provenance._git_commit",
        lambda root: "moving-head",
    )

    assert runtime_provenance()["genomepuzzle_git_commit"] == "frozen-commit"
