"""
Long-read and hybrid dataset QC summaries.
"""

import csv
import json
import os
from statistics import median

from genomepuzzle.runtime import open_maybe_gzip


def summarize_fastq(fastq_path):
    read_count = 0
    total_bases = 0
    total_quality = 0
    quality_bases = 0
    lengths = []

    with open_maybe_gzip(fastq_path, "rt") as handle:
        while True:
            header = handle.readline().strip()
            sequence = handle.readline().strip()
            plus = handle.readline().strip()
            quality = handle.readline().strip()
            if not header:
                break
            if not sequence or not plus:
                raise ValueError("Malformed FASTQ record in {0}".format(fastq_path))
            read_count += 1
            read_length = len(sequence)
            lengths.append(read_length)
            total_bases += read_length
            total_quality += sum(ord(char) - 33 for char in quality)
            quality_bases += len(quality)

    mean_length = round(total_bases / float(read_count), 2) if read_count else 0.0
    mean_quality = round(total_quality / float(quality_bases), 2) if quality_bases else 0.0
    n50 = calculate_n50(lengths)
    return {
        "read_count": read_count,
        "total_bases": total_bases,
        "mean_length": mean_length,
        "mean_quality": mean_quality,
        "n50_length": n50,
    }


def calculate_n50(lengths):
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    cumulative = 0
    for length in ordered:
        cumulative += length
        if cumulative >= half:
            return length
    return 0


def load_hybrid_sample_sheet(sample_sheet):
    with open(sample_sheet, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("sample sheet contains no records")
    required = {"sample_name", "r1", "r2", "long_reads"}
    missing = required.difference(rows[0].keys())
    if missing:
        raise ValueError(
            "sample sheet is missing required columns: {0}".format(
                ", ".join(sorted(missing))
            )
        )
    return rows


def resolve_dataset_path(dataset_dir, filename):
    if os.path.isabs(filename):
        return filename
    return os.path.join(dataset_dir, filename)


def classify_sample(short_read_count, long_read_count, short_baseline, long_baseline):
    flags = []
    if short_baseline and short_read_count < (0.25 * short_baseline):
        flags.append("LOW_SHORT_READ_COUNT")
    if long_baseline and long_read_count < (0.25 * long_baseline):
        flags.append("LOW_LONG_READ_COUNT")
    return flags or ["OK"]


def summarize_hybrid_dataset(sample_sheet, dataset_dir=None):
    rows = load_hybrid_sample_sheet(sample_sheet)
    dataset_dir = dataset_dir or os.path.dirname(sample_sheet) or "."

    summaries = []
    for row in rows:
        short_r1_path = resolve_dataset_path(dataset_dir, row["r1"])
        short_r2_path = resolve_dataset_path(dataset_dir, row["r2"])
        long_reads_path = resolve_dataset_path(dataset_dir, row["long_reads"])

        r1_stats = summarize_fastq(short_r1_path)
        r2_stats = summarize_fastq(short_r2_path)
        long_stats = summarize_fastq(long_reads_path)

        summaries.append(
            {
                "sample_name": row["sample_name"],
                "species": row.get("species", ""),
                "reference_accession": row.get("reference_accession", ""),
                "r1": row["r1"],
                "r2": row["r2"],
                "long_reads": row["long_reads"],
                "short_read_count": min(r1_stats["read_count"], r2_stats["read_count"]),
                "short_total_bases": r1_stats["total_bases"] + r2_stats["total_bases"],
                "short_mean_read_length": round(
                    (r1_stats["mean_length"] + r2_stats["mean_length"]) / 2.0, 2
                ),
                "short_mean_quality": round(
                    (r1_stats["mean_quality"] + r2_stats["mean_quality"]) / 2.0, 2
                ),
                "long_read_count": long_stats["read_count"],
                "long_total_bases": long_stats["total_bases"],
                "long_mean_read_length": long_stats["mean_length"],
                "long_mean_quality": long_stats["mean_quality"],
                "long_n50_length": long_stats["n50_length"],
            }
        )

    short_baseline = median([row["short_read_count"] for row in summaries])
    long_baseline = median([row["long_read_count"] for row in summaries])
    for row in summaries:
        row["flags"] = ";".join(
            classify_sample(
                row["short_read_count"],
                row["long_read_count"],
                short_baseline,
                long_baseline,
            )
        )
    return summaries


def _load_csv_rows(path):
    with open(path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def expected_flags_for_error(error_type):
    mapping = {
        "NORMAL": ["OK"],
        "LOW_SHORT_COVERAGE": ["LOW_SHORT_READ_COUNT"],
        "LOW_LONG_COVERAGE": ["LOW_LONG_READ_COUNT"],
        "LONG_READ_QUALITY": [],
        "CONTAMINATED": [],
    }
    return mapping.get(error_type, [])


def compare_qc_to_manifest(qc_rows, manifest_path):
    manifest_rows = _load_csv_rows(manifest_path)
    qc_by_sample = {row["sample_name"]: row for row in qc_rows}
    comparison_rows = []
    for manifest_row in manifest_rows:
        sample_name = manifest_row["sample_name"]
        qc_row = qc_by_sample.get(sample_name)
        if qc_row is None:
            comparison_rows.append(
                {
                    "sample_name": sample_name,
                    "error_type": manifest_row["error_type"],
                    "expected_flags": "",
                    "observed_flags": "MISSING_QC",
                    "status": "missing_qc",
                    "notes": "No QC summary row found for this sample",
                }
            )
            continue
        observed_flags = [
            flag for flag in qc_row.get("flags", "").split(";") if flag and flag != "OK"
        ]
        expected_flags = expected_flags_for_error(manifest_row["error_type"])
        if not expected_flags:
            status = "not_assessed"
            notes = "Current QC summary does not directly assess this implant type"
        elif all(flag in observed_flags for flag in expected_flags):
            status = "detected"
            notes = "Observed QC flags match the expected implant signal"
        else:
            status = "missed"
            notes = "Expected implant signal was not fully present in QC flags"
        comparison_rows.append(
            {
                "sample_name": sample_name,
                "error_type": manifest_row["error_type"],
                "expected_flags": ";".join(expected_flags) or "N/A",
                "observed_flags": qc_row.get("flags", "OK"),
                "status": status,
                "notes": notes,
            }
        )
    return comparison_rows


def write_qc_outputs(sample_rows, output_csv, output_json=None):
    fieldnames = [
        "sample_name",
        "species",
        "reference_accession",
        "r1",
        "r2",
        "long_reads",
        "short_read_count",
        "short_total_bases",
        "short_mean_read_length",
        "short_mean_quality",
        "long_read_count",
        "long_total_bases",
        "long_mean_read_length",
        "long_mean_quality",
        "long_n50_length",
        "flags",
    ]
    with open(output_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sample_rows)
    if output_json:
        with open(output_json, "w", encoding="utf-8") as handle:
            json.dump(sample_rows, handle, indent=2)


def write_report_outputs(report_rows, output_csv, output_json=None):
    fieldnames = [
        "sample_name",
        "error_type",
        "expected_flags",
        "observed_flags",
        "status",
        "notes",
    ]
    with open(output_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)
    if output_json:
        with open(output_json, "w", encoding="utf-8") as handle:
            json.dump(report_rows, handle, indent=2)
