import pytest

typer = pytest.importorskip("typer")
from typer.testing import CliRunner

from genomepuzzle.main import app


runner = CliRunner()


def test_cli_shows_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "long" in result.stdout
    assert "legacy" in result.stdout
    assert "release" in result.stdout


def test_hybrid_help_mentions_mode():
    result = runner.invoke(app, ["legacy", "hybrid", "--help"])
    assert result.exit_code == 0
    assert "mode" in result.stdout
    assert "challenge" in result.stdout


def test_long_qc_help_mentions_sample_sheet():
    result = runner.invoke(app, ["long", "qc", "--help"])
    assert result.exit_code == 0
    assert "sample-sheet" in result.stdout
    assert "output-csv" in result.stdout


def test_long_report_help_mentions_manifest():
    result = runner.invoke(app, ["long", "report", "--help"])
    assert result.exit_code == 0
    assert "manifest" in result.stdout


def test_long_hybrid_slurm_help_mentions_partition():
    result = runner.invoke(app, ["legacy", "hybrid-slurm", "--help"])
    assert result.exit_code == 0
    assert "partition" in result.stdout


def test_release_validate_spec_writes_private_mapping(tmp_path):
    spec = tmp_path / "release.toml"
    spec.write_text(
        """
schema_version = "1.0"
release_id = "test-typing-practice"
exercise = "typing"
mode = "practice"
master_seed = 42

[[samples]]
source_id = "GCA_000001.1"
public_id = "Sample_fixed123"
""",
        encoding="utf-8",
    )
    output = tmp_path / "resolved.json"
    result = runner.invoke(
        app,
        [
            "release",
            "validate-spec",
            "--spec",
            str(spec),
            "--output-json",
            str(output),
        ],
    )
    assert result.exit_code == 0
    assert output.exists()
    assert "Sample_fixed123" in output.read_text(encoding="utf-8")
