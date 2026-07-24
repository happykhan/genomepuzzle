import gzip
from pathlib import Path

import pytest

from genomepuzzle.pilot_review import (
    fasta_metrics,
    fastq_metrics,
    stage_named_assemblies,
)


def test_fastq_metrics_streams_read_statistics(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt", encoding="ascii") as handle:
        handle.write("@one\nACGT\n+\nIIII\n@two\nGG\n+\n!!\n")

    assert fastq_metrics(path) == {
        "records": 2,
        "bases": 6,
        "minimum_length": 2,
        "maximum_length": 4,
        "mean_length": 3.0,
        "mean_quality": 26.667,
        "gc_fraction": 0.666667,
    }


def test_fastq_metrics_rejects_truncated_records(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt", encoding="ascii") as handle:
        handle.write("@one\nACGT\n+\n")
    with pytest.raises(ValueError, match="truncated"):
        fastq_metrics(path)


def test_fasta_metrics_reports_n50(tmp_path):
    path = tmp_path / "assembly.fasta"
    path.write_text(">one\nAAAAAA\n>two\nGCNN\n>three\nA\n", encoding="ascii")

    assert fasta_metrics(path) == {
        "contigs": 3,
        "total_bases": 11,
        "n50": 6,
        "maximum_contig": 6,
        "ambiguous_bases": 2,
        "ambiguous_fraction": 0.18181818,
    }


def test_stage_named_assemblies_preserves_sample_identity(tmp_path):
    first = tmp_path / "one" / "contigs.fasta"
    second = tmp_path / "two" / "contigs.fasta"
    first.parent.mkdir()
    second.parent.mkdir()
    first.write_text(">one\nAAAA\n", encoding="ascii")
    second.write_text(">two\nCCCC\n", encoding="ascii")

    staged = stage_named_assemblies(
        tmp_path / "analysis-inputs",
        [("Sample_A", str(first)), ("Sample_B", str(second))],
    )

    assert [Path(path) for path in staged] == [
        tmp_path / "analysis-inputs" / "Sample_A.fasta",
        tmp_path / "analysis-inputs" / "Sample_B.fasta",
    ]
    assert (tmp_path / "analysis-inputs" / "Sample_A.fasta").read_text() == (
        ">one\nAAAA\n"
    )
