import gzip

import pytest

from genomepuzzle.pilot_review import fasta_metrics, fastq_metrics


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
