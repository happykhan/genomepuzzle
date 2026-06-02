"""
Generate hybrid assembly datasets with explicit implanted errors.
"""

import csv
import gzip
import hashlib
import json
import logging
import os
import random
import shutil
from collections import Counter

from genomepuzzle.create_error import (
    concatenate_fastqs,
    degrade_quality,
    subsample_paired_fastq,
    subsample_single_fastq,
    subsample_single_fastq_by_count,
)
from genomepuzzle.sample_util import run_badread
from genomepuzzle.simulate_reads import cleanup_output_dir, fetch_assembly, run_art


HYBRID_ERROR_TYPES = [
    "NORMAL",
    "LOW_SHORT_COVERAGE",
    "LOW_LONG_COVERAGE",
    "LONG_READ_QUALITY",
    "CONTAMINATED",
]


def count_reads(fastq_file):
    count = 0
    with gzip.open(fastq_file, "rt") as handle:
        for _ in handle:
            count += 1
    return count // 4


def stable_public_name(accession, random_seed):
    digest = hashlib.md5(
        "{accession}:{seed}".format(accession=accession, seed=random_seed).encode("utf-8")
    ).hexdigest()[:10]
    return "Sample_{digest}".format(digest=digest)


def build_hybrid_error_plan(num_samples, mode="challenge", random_seed=42):
    """
    Build a reproducible hybrid implant plan.
    """
    if mode not in ["practice", "challenge", "none"]:
        raise ValueError("mode must be one of practice, challenge, or none")
    if mode == "none":
        return ["NORMAL"] * num_samples

    if mode == "practice":
        plan = [
            "CONTAMINATED",
            "LOW_SHORT_COVERAGE",
            "LOW_LONG_COVERAGE",
            "LONG_READ_QUALITY",
        ] + ["NORMAL"] * max(0, num_samples - 4)
    else:
        plan = [
            "CONTAMINATED",
            "CONTAMINATED",
            "LOW_SHORT_COVERAGE",
            "LOW_SHORT_COVERAGE",
            "LOW_LONG_COVERAGE",
            "LOW_LONG_COVERAGE",
            "LONG_READ_QUALITY",
            "LONG_READ_QUALITY",
        ] + ["NORMAL"] * max(0, num_samples - 8)
    plan = plan[:num_samples]
    random.seed(random_seed)
    random.shuffle(plan)
    return plan


def read_records(samplelist):
    with open(samplelist, encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle)]


def get_assembly_path(output_dir, accession):
    assembly_dir = os.path.join(output_dir, "ncbi_dataset", "data", accession)
    fasta_files = [
        os.path.join(assembly_dir, filename)
        for filename in os.listdir(assembly_dir)
        if filename.endswith(".fna") or filename.endswith(".fasta")
    ]
    return fasta_files[0]


def simulate_base_reads(record, output_dir, random_seed=42):
    accession = record["accession"]
    public_name = record["public_name"]
    assembly_path = get_assembly_path(output_dir, accession)
    short_r1 = os.path.join(output_dir, "{name}_R1.fastq.gz".format(name=public_name))
    short_r2 = os.path.join(output_dir, "{name}_R2.fastq.gz".format(name=public_name))
    sample = {
        "public_name": public_name,
        "platform": "HS25",
        "read_length": 150,
        "coverage": float(record.get("short_read_coverage", 40)),
        "fragment_length": 300,
        "standard_deviation": 50,
        "random_seed": random_seed,
    }
    run_art(sample, output_dir, assembly_path, short_r1, short_r2)
    long_output = run_badread(public_name, assembly_path, output_dir)["long_reads_path"]
    return short_r1, short_r2, long_output


def choose_contaminant(target_record, contaminant_records, random_seed):
    compatible = [
        row for row in contaminant_records
        if row.get("organism_organismname") != target_record.get("organism_organismname")
    ]
    if not compatible:
        compatible = [row for row in contaminant_records if row["accession"] != target_record["accession"]]
    random.seed(random_seed)
    return random.choice(compatible) if compatible else None


def apply_low_short_coverage(r1, r2, output_prefix, random_seed):
    fraction = random.uniform(0.08, 0.25)
    new_r1 = "{prefix}_R1.fastq.gz".format(prefix=output_prefix)
    new_r2 = "{prefix}_R2.fastq.gz".format(prefix=output_prefix)
    subsample_paired_fastq(r1, r2, new_r1, new_r2, fraction, random_seed)
    return new_r1, new_r2, {
        "error_type": "LOW_SHORT_COVERAGE",
        "severity": round(fraction, 4),
        "notes": "Paired reads subsampled to {:.1f}% of original".format(fraction * 100),
    }


def apply_low_long_coverage(long_reads, output_prefix, random_seed):
    fraction = random.uniform(0.08, 0.25)
    new_long = "{prefix}_long.fastq.gz".format(prefix=output_prefix)
    subsample_single_fastq(long_reads, new_long, fraction, random_seed)
    return new_long, {
        "error_type": "LOW_LONG_COVERAGE",
        "severity": round(fraction, 4),
        "notes": "Long reads subsampled to {:.1f}% of original".format(fraction * 100),
    }


def apply_long_read_quality(long_reads, output_prefix, random_seed):
    min_quality = random.randint(3, 10)
    max_quality = random.randint(12, 22)
    new_long = "{prefix}_long.fastq.gz".format(prefix=output_prefix)
    degrade_quality(long_reads, new_long, min_quality=min_quality, max_quality=max_quality, random_seed=random_seed)
    return new_long, {
        "error_type": "LONG_READ_QUALITY",
        "severity": "{min_q}-{max_q}".format(min_q=min_quality, max_q=max_quality),
        "notes": "Long-read qualities degraded",
    }


def apply_contamination(
    record,
    r1,
    r2,
    long_reads,
    contaminant_record,
    output_dir,
    output_prefix,
    random_seed,
):
    contamination_fraction = random.uniform(0.15, 0.45)
    contaminant_public_name = "contaminant_{name}".format(name=record["public_name"])
    contaminant_tmp = dict(contaminant_record)
    contaminant_tmp["public_name"] = contaminant_public_name
    contam_r1, contam_r2, contam_long = simulate_base_reads(contaminant_tmp, output_dir, random_seed + 101)

    original_short_reads = count_reads(r1)
    contam_short_reads = max(1, int(round(original_short_reads * contamination_fraction)))
    clean_short_reads = max(1, original_short_reads - contam_short_reads)

    original_long_reads = count_reads(long_reads)
    contam_long_reads = max(1, int(round(original_long_reads * contamination_fraction)))
    clean_long_reads = max(1, original_long_reads - contam_long_reads)

    clean_r1 = os.path.join(output_dir, "tmp_clean_r1.fastq.gz")
    clean_r2 = os.path.join(output_dir, "tmp_clean_r2.fastq.gz")
    dirty_r1 = os.path.join(output_dir, "tmp_dirty_r1.fastq.gz")
    dirty_r2 = os.path.join(output_dir, "tmp_dirty_r2.fastq.gz")
    clean_long = os.path.join(output_dir, "tmp_clean_long.fastq.gz")
    dirty_long = os.path.join(output_dir, "tmp_dirty_long.fastq.gz")

    subsample_single_fastq_by_count(long_reads, clean_long, clean_long_reads, random_seed)
    subsample_single_fastq_by_count(contam_long, dirty_long, contam_long_reads, random_seed)
    subsample_single_fastq_by_count(contam_r1, dirty_r1, contam_short_reads, random_seed)
    subsample_single_fastq_by_count(contam_r2, dirty_r2, contam_short_reads, random_seed)
    subsample_single_fastq_by_count(r1, clean_r1, clean_short_reads, random_seed)
    subsample_single_fastq_by_count(r2, clean_r2, clean_short_reads, random_seed)

    final_r1 = "{prefix}_R1.fastq.gz".format(prefix=output_prefix)
    final_r2 = "{prefix}_R2.fastq.gz".format(prefix=output_prefix)
    final_long = "{prefix}_long.fastq.gz".format(prefix=output_prefix)
    concatenate_fastqs([clean_r1, dirty_r1], final_r1)
    concatenate_fastqs([clean_r2, dirty_r2], final_r2)
    concatenate_fastqs([clean_long, dirty_long], final_long)

    for temp_file in [
        clean_r1,
        clean_r2,
        dirty_r1,
        dirty_r2,
        clean_long,
        dirty_long,
        contam_r1,
        contam_r2,
        contam_long,
    ]:
        if os.path.exists(temp_file):
            os.remove(temp_file)

    return final_r1, final_r2, final_long, {
        "error_type": "CONTAMINATED",
        "severity": round(contamination_fraction, 4),
        "notes": "Mixed with contaminant accession {accession}".format(accession=contaminant_record["accession"]),
    }


def write_csv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def create_hybrid_dataset(
    output_dir,
    samplelist,
    contamination_list=None,
    mode="challenge",
    random_seed=42,
):
    """
    Create a hybrid dataset with implanted errors and a private manifest.
    """
    records = read_records(samplelist)
    if not records:
        raise ValueError("samplelist contains no records")
    os.makedirs(output_dir, exist_ok=True)

    contaminant_records = read_records(contamination_list) if contamination_list else list(records)
    all_accessions = list({row["accession"] for row in records + contaminant_records})
    fetch_assembly(all_accessions, output_dir)

    error_plan = build_hybrid_error_plan(len(records), mode=mode, random_seed=random_seed)
    answer_rows = []
    sample_rows = []
    manifest_rows = []

    for index, (record, error_type) in enumerate(zip(records, error_plan)):
        seed = random_seed + index
        public_name = record.get("public_name") or stable_public_name(record["accession"], seed)
        record["public_name"] = public_name
        species = record.get("organism_organismname") or record.get("species") or "Unknown"
        r1, r2, long_reads = simulate_base_reads(record, output_dir, random_seed=seed)
        metadata = {
            "error_type": "NORMAL",
            "severity": "none",
            "notes": "No implant",
        }

        final_r1 = r1
        final_r2 = r2
        final_long = long_reads
        final_prefix = os.path.join(output_dir, public_name)

        if error_type == "LOW_SHORT_COVERAGE":
            final_r1, final_r2, metadata = apply_low_short_coverage(r1, r2, final_prefix, seed)
            if os.path.exists(r1):
                os.remove(r1)
            if os.path.exists(r2):
                os.remove(r2)
        elif error_type == "LOW_LONG_COVERAGE":
            final_long, metadata = apply_low_long_coverage(long_reads, final_prefix, seed)
            if os.path.exists(long_reads):
                os.remove(long_reads)
        elif error_type == "LONG_READ_QUALITY":
            final_long, metadata = apply_long_read_quality(long_reads, final_prefix, seed)
            if os.path.exists(long_reads):
                os.remove(long_reads)
        elif error_type == "CONTAMINATED":
            contaminant_record = choose_contaminant(record, contaminant_records, seed)
            if not contaminant_record:
                logging.warning("No contaminant available for %s, leaving sample as NORMAL", public_name)
            else:
                final_r1, final_r2, final_long, metadata = apply_contamination(
                    record,
                    r1,
                    r2,
                    long_reads,
                    contaminant_record,
                    output_dir,
                    final_prefix,
                    seed,
                )
                for old_path in [r1, r2, long_reads]:
                    if os.path.exists(old_path):
                        os.remove(old_path)

        expected_qc = "PASSED" if metadata["error_type"] == "NORMAL" else "FAILED"
        answer_rows.append({
            "public_name": public_name,
            "species": species,
            "reference_accession": record["accession"],
            "tax_classification": species,
            "assembler": "Unknown",
            "qc": expected_qc,
            "notes": metadata["notes"],
            "error_type": metadata["error_type"],
            "severity": metadata["severity"],
        })
        sample_rows.append({
            "sample_name": public_name,
            "reference_accession": record["accession"],
            "species": species,
            "tax_classification": "",
            "r1": os.path.basename(final_r1),
            "r2": os.path.basename(final_r2),
            "long_reads": os.path.basename(final_long),
            "assembler": "",
            "qc": "",
            "notes": "",
        })
        manifest_rows.append({
            "sample_name": public_name,
            "reference_accession": record["accession"],
            "species": species,
            "error_type": metadata["error_type"],
            "severity": metadata["severity"],
            "notes": metadata["notes"],
            "short_read_count": count_reads(final_r1),
            "long_read_count": count_reads(final_long),
        })

    write_csv(
        os.path.join(output_dir, "answer_sheet.csv"),
        answer_rows,
        ["public_name", "species", "reference_accession", "tax_classification", "assembler", "qc", "notes", "error_type", "severity"],
    )
    write_csv(
        os.path.join(output_dir, "sample_sheet.csv"),
        sample_rows,
        ["sample_name", "reference_accession", "species", "tax_classification", "r1", "r2", "long_reads", "assembler", "qc", "notes"],
    )
    write_csv(
        os.path.join(output_dir, "implant_manifest.csv"),
        manifest_rows,
        ["sample_name", "reference_accession", "species", "error_type", "severity", "notes", "short_read_count", "long_read_count"],
    )

    summary = Counter(plan for plan in error_plan)
    with open(os.path.join(output_dir, "implant_manifest.json"), "w", encoding="utf-8") as handle:
        json.dump(
            {
                "mode": mode,
                "random_seed": random_seed,
                "summary": dict(summary),
                "samples": manifest_rows,
            },
            handle,
            indent=2,
        )
    cleanup_output_dir(output_dir)

