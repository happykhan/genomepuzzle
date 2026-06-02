import pytest

typer = pytest.importorskip("typer")
from typer.testing import CliRunner

from genomepuzzle.main import app


runner = CliRunner()


def test_cli_shows_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "short" in result.stdout
    assert "long" in result.stdout
    assert "rapid" in result.stdout


def test_hybrid_help_mentions_mode():
    result = runner.invoke(app, ["long", "hybrid", "--help"])
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
    result = runner.invoke(app, ["long", "hybrid-slurm", "--help"])
    assert result.exit_code == 0
    assert "partition" in result.stdout
