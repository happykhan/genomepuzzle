import csv
import gzip
import json
import os

from genomepuzzle.long_qc import (
    summarize_fastq,
    summarize_hybrid_dataset,
    write_qc_outputs,
)


def write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for header, sequence, quality in records:
            handle.write("@{0}\n{1}\n+\n{2}\n".format(header, sequence, quality))


def test_summarize_fastq(tmp_path):
    fastq_path = tmp_path / "reads.fastq.gz"
    write_fastq(
        fastq_path,
        [
            ("r1", "ACTG", "IIII"),
            ("r2", "ACTGACTG", "JJJJJJJJ"),
        ],
    )
    stats = summarize_fastq(str(fastq_path))
    assert stats["read_count"] == 2
    assert stats["total_bases"] == 12
    assert stats["mean_length"] == 6.0
    assert stats["n50_length"] == 8


def test_summarize_hybrid_dataset_flags_low_long_read_count(tmp_path):
    sample_sheet = tmp_path / "sample_sheet.csv"
    write_fastq(tmp_path / "sample1_R1.fastq.gz", [("a", "ACTG", "IIII")] * 10)
    write_fastq(tmp_path / "sample1_R2.fastq.gz", [("a", "ACTG", "IIII")] * 10)
    write_fastq(tmp_path / "sample1_long.fastq.gz", [("a", "ACTGACTG", "IIIIIIII")] * 12)
    write_fastq(tmp_path / "sample2_R1.fastq.gz", [("b", "ACTG", "IIII")] * 10)
    write_fastq(tmp_path / "sample2_R2.fastq.gz", [("b", "ACTG", "IIII")] * 10)
    write_fastq(tmp_path / "sample2_long.fastq.gz", [("b", "ACTGACTG", "IIIIIIII")] * 1)

    with open(sample_sheet, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_name",
                "reference_accession",
                "species",
                "tax_classification",
                "r1",
                "r2",
                "long_reads",
                "assembler",
                "qc",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_name": "sample1",
                "reference_accession": "GCF_1",
                "species": "Kp",
                "tax_classification": "",
                "r1": "sample1_R1.fastq.gz",
                "r2": "sample1_R2.fastq.gz",
                "long_reads": "sample1_long.fastq.gz",
                "assembler": "",
                "qc": "",
                "notes": "",
            }
        )
        writer.writerow(
            {
                "sample_name": "sample2",
                "reference_accession": "GCF_2",
                "species": "Kp",
                "tax_classification": "",
                "r1": "sample2_R1.fastq.gz",
                "r2": "sample2_R2.fastq.gz",
                "long_reads": "sample2_long.fastq.gz",
                "assembler": "",
                "qc": "",
                "notes": "",
            }
        )

    rows = summarize_hybrid_dataset(str(sample_sheet), dataset_dir=str(tmp_path))
    flags = {row["sample_name"]: row["flags"] for row in rows}
    assert flags["sample1"] == "OK"
    assert "LOW_LONG_READ_COUNT" in flags["sample2"]


def test_write_qc_outputs(tmp_path):
    rows = [
        {
            "sample_name": "sample1",
            "species": "Kp",
            "reference_accession": "GCF_1",
            "r1": "a.fastq.gz",
            "r2": "b.fastq.gz",
            "long_reads": "c.fastq.gz",
            "short_read_count": 10,
            "short_total_bases": 3000,
            "short_mean_read_length": 150.0,
            "short_mean_quality": 30.0,
            "long_read_count": 50,
            "long_total_bases": 100000,
            "long_mean_read_length": 2000.0,
            "long_mean_quality": 14.0,
            "long_n50_length": 2500,
            "flags": "OK",
        }
    ]
    output_csv = tmp_path / "qc.csv"
    output_json = tmp_path / "qc.json"
    write_qc_outputs(rows, str(output_csv), str(output_json))
    assert output_csv.exists()
    assert json.loads(output_json.read_text(encoding="utf-8"))[0]["sample_name"] == "sample1"
