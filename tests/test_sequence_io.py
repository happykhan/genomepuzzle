import gzip

import pytest

from genomepuzzle.sequence_io import (
    anonymize_paired_fastq,
    fastq_records,
    is_anonymous_fastq_header,
)


def _write(path, header):
    with gzip.open(path, "wt") as handle:
        handle.write("@{0}\nACGT\n+\nIIII\n".format(header))


def test_anonymize_paired_fastq_reheaders_and_validates_pair(tmp_path):
    source_r1 = tmp_path / "source_R1.fastq.gz"
    source_r2 = tmp_path / "source_R2.fastq.gz"
    output_r1 = tmp_path / "Sample_x_R1.fastq.gz"
    output_r2 = tmp_path / "Sample_x_R2.fastq.gz"
    _write(source_r1, "secret/1")
    _write(source_r2, "secret/2")

    assert anonymize_paired_fastq(
        source_r1, source_r2, output_r1, output_r2, "Sample_x"
    ) == 1
    assert next(fastq_records(output_r1))[0] == "@Sample_x_read_000000001/1\n"
    assert "secret" not in gzip.open(output_r1, "rt").read()


def test_fastq_validation_rejects_sequence_quality_mismatch(tmp_path):
    path = tmp_path / "bad.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@read\nACGT\n+\nIII\n")
    with pytest.raises(ValueError, match="length mismatch"):
        list(fastq_records(path))


def test_badread_generic_headers_are_safe_but_source_headers_are_not():
    assert is_anonymous_fastq_header(
        "@uuid junk_seq length=1000\n", "Sample_x"
    )
    assert is_anonymous_fastq_header(
        "@uuid random_seq length=1000\n", "Sample_x"
    )
    assert is_anonymous_fastq_header(
        "@uuid Sample_x_contig_00001,+strand,1-1000\n", "Sample_x"
    )
    assert not is_anonymous_fastq_header(
        "@uuid CABFXZ010000001.1,+strand,1-1000\n", "Sample_x"
    )
