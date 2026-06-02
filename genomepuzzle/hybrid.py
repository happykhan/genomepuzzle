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

SHORT_COVERAGE_FRACTION = 0.15
LONG_COVERAGE_FRACTION = 0.15
LONG_READ_MIN_QUALITY = 5
LONG_READ_MAX_QUALITY = 14
CONTAMINATION_FRACTION = 0.30


class HybridImplant(object):
    def __init__(self, error_type, severity, notes):
        self.error_type = error_type
        self.severity = severity
        self.notes = notes


class HybridSampleContext(object):
    def __init__(self, record, public_name, species, seed, implant):
        self.record = record
        self.public_name = public_name
        self.species = species
        self.seed = seed
        self.r1 = None
        self.r2 = None
        self.long_reads = None
        self.implant = implant
        self.implant_count = 0


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
    rng = random.Random(random_seed)
    rng.shuffle(plan)
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


def simulate_base_reads(context, output_dir):
    assembly_path = get_assembly_path(output_dir, context.record["accession"])
    short_r1 = os.path.join(output_dir, "{name}_R1.fastq.gz".format(name=context.public_name))
    short_r2 = os.path.join(output_dir, "{name}_R2.fastq.gz".format(name=context.public_name))
    sample = {
        "public_name": context.public_name,
        "platform": "HS25",
        "read_length": 150,
        "coverage": float(context.record.get("short_read_coverage", 40)),
        "fragment_length": 300,
        "standard_deviation": 50,
        "random_seed": context.seed,
    }
    run_art(sample, output_dir, assembly_path, short_r1, short_r2)
    long_output = run_badread(context.public_name, assembly_path, output_dir)["long_reads_path"]
    context.r1 = short_r1
    context.r2 = short_r2
    context.long_reads = long_output


def choose_contaminant(target_record, contaminant_records, rng):
    compatible = [
        row for row in contaminant_records
        if row.get("organism_organismname") != target_record.get("organism_organismname")
    ]
    if not compatible:
        compatible = [row for row in contaminant_records if row["accession"] != target_record["accession"]]
    return rng.choice(compatible) if compatible else None


def _make_temp_path(output_dir, public_name, suffix):
    return os.path.join(output_dir, "{name}_{suffix}".format(name=public_name, suffix=suffix))


def apply_low_short_coverage(context, output_dir):
    fraction = SHORT_COVERAGE_FRACTION
    new_r1 = os.path.join(output_dir, "{name}_R1.fastq.gz".format(name=context.public_name))
    new_r2 = os.path.join(output_dir, "{name}_R2.fastq.gz".format(name=context.public_name))
    subsample_paired_fastq(context.r1, context.r2, new_r1, new_r2, fraction, context.seed)
    os.remove(context.r1)
    os.remove(context.r2)
    context.r1 = new_r1
    context.r2 = new_r2
    context.implant = HybridImplant(
        "LOW_SHORT_COVERAGE",
        round(fraction, 4),
        "Paired reads subsampled to {:.0f}% of original".format(fraction * 100),
    )
    context.implant_count += 1


def apply_low_long_coverage(context, output_dir):
    fraction = LONG_COVERAGE_FRACTION
    new_long = os.path.join(output_dir, "{name}_long.fastq.gz".format(name=context.public_name))
    subsample_single_fastq(context.long_reads, new_long, fraction, context.seed)
    os.remove(context.long_reads)
    context.long_reads = new_long
    context.implant = HybridImplant(
        "LOW_LONG_COVERAGE",
        round(fraction, 4),
        "Long reads subsampled to {:.0f}% of original".format(fraction * 100),
    )
    context.implant_count += 1


def apply_long_read_quality(context, output_dir):
    min_quality = LONG_READ_MIN_QUALITY
    max_quality = LONG_READ_MAX_QUALITY
    new_long = os.path.join(output_dir, "{name}_long.fastq.gz".format(name=context.public_name))
    degrade_quality(
        context.long_reads,
        new_long,
        min_quality=min_quality,
        max_quality=max_quality,
        random_seed=context.seed,
    )
    os.remove(context.long_reads)
    context.long_reads = new_long
    context.implant = HybridImplant(
        "LONG_READ_QUALITY",
        "{min_q}-{max_q}".format(min_q=min_quality, max_q=max_quality),
        "Long-read qualities degraded",
    )
    context.implant_count += 1


def apply_contamination(context, contaminant_record, output_dir):
    contamination_fraction = CONTAMINATION_FRACTION
    contaminant_context = HybridSampleContext(
        record=dict(contaminant_record),
        public_name="contaminant_{name}_{accession}".format(
            name=context.public_name,
            accession=contaminant_record["accession"].replace(".", "_"),
        ),
        species=contaminant_record.get("organism_organismname")
        or contaminant_record.get("species")
        or "Unknown",
        seed=context.seed + 101,
    )
    simulate_base_reads(contaminant_context, output_dir)

    original_short_reads = count_reads(context.r1)
    contaminant_short_reads = max(1, int(round(original_short_reads * contamination_fraction)))
    clean_short_reads = max(1, original_short_reads - contaminant_short_reads)
    original_long_reads = count_reads(context.long_reads)
    contaminant_long_reads = max(1, int(round(original_long_reads * contamination_fraction)))
    clean_long_reads = max(1, original_long_reads - contaminant_long_reads)

    temp_clean_r1 = _make_temp_path(output_dir, context.public_name, "tmp_clean_r1.fastq.gz")
    temp_clean_r2 = _make_temp_path(output_dir, context.public_name, "tmp_clean_r2.fastq.gz")
    temp_dirty_r1 = _make_temp_path(output_dir, context.public_name, "tmp_dirty_r1.fastq.gz")
    temp_dirty_r2 = _make_temp_path(output_dir, context.public_name, "tmp_dirty_r2.fastq.gz")
    temp_clean_long = _make_temp_path(output_dir, context.public_name, "tmp_clean_long.fastq.gz")
    temp_dirty_long = _make_temp_path(output_dir, context.public_name, "tmp_dirty_long.fastq.gz")

    subsample_single_fastq_by_count(context.long_reads, temp_clean_long, clean_long_reads, context.seed)
    subsample_single_fastq_by_count(contaminant_context.long_reads, temp_dirty_long, contaminant_long_reads, context.seed)
    subsample_single_fastq_by_count(contaminant_context.r1, temp_dirty_r1, contaminant_short_reads, context.seed)
    subsample_single_fastq_by_count(contaminant_context.r2, temp_dirty_r2, contaminant_short_reads, context.seed)
    subsample_single_fastq_by_count(context.r1, temp_clean_r1, clean_short_reads, context.seed)
    subsample_single_fastq_by_count(context.r2, temp_clean_r2, clean_short_reads, context.seed)

    final_r1 = os.path.join(output_dir, "{name}_R1.fastq.gz".format(name=context.public_name))
    final_r2 = os.path.join(output_dir, "{name}_R2.fastq.gz".format(name=context.public_name))
    final_long = os.path.join(output_dir, "{name}_long.fastq.gz".format(name=context.public_name))
    concatenate_fastqs([temp_clean_r1, temp_dirty_r1], final_r1)
    concatenate_fastqs([temp_clean_r2, temp_dirty_r2], final_r2)
    concatenate_fastqs([temp_clean_long, temp_dirty_long], final_long)

    for path in [
        context.r1,
        context.r2,
        context.long_reads,
        contaminant_context.r1,
        contaminant_context.r2,
        contaminant_context.long_reads,
        temp_clean_r1,
        temp_clean_r2,
        temp_dirty_r1,
        temp_dirty_r2,
        temp_clean_long,
        temp_dirty_long,
    ]:
        if os.path.exists(path):
            os.remove(path)

    context.r1 = final_r1
    context.r2 = final_r2
    context.long_reads = final_long
    context.implant = HybridImplant(
        "CONTAMINATED",
        round(contamination_fraction, 4),
        "Mixed with contaminant accession {accession}".format(
            accession=contaminant_record["accession"]
        ),
    )
    context.implant_count += 1


def write_csv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_sample_context(record, index, random_seed):
    seed = random_seed + index
    public_name = record.get("public_name") or stable_public_name(record["accession"], seed)
    species = record.get("organism_organismname") or record.get("species") or "Unknown"
    return HybridSampleContext(
        record=dict(record),
        public_name=public_name,
        species=species,
        seed=seed,
        implant=HybridImplant("NORMAL", "none", "No implant"),
    )


def implant_context(context, error_type, contaminant_records, output_dir):
    if error_type == "LOW_SHORT_COVERAGE":
        apply_low_short_coverage(context, output_dir)
    elif error_type == "LOW_LONG_COVERAGE":
        apply_low_long_coverage(context, output_dir)
    elif error_type == "LONG_READ_QUALITY":
        apply_long_read_quality(context, output_dir)
    elif error_type == "CONTAMINATED":
        rng = random.Random(context.seed)
        contaminant_record = choose_contaminant(context.record, contaminant_records, rng)
        if contaminant_record:
            apply_contamination(context, contaminant_record, output_dir)
        else:
            logging.warning("No contaminant available for %s, leaving sample as NORMAL", context.public_name)
    if context.implant_count > 1:
        raise ValueError(
            "Hybrid sample {name} received multiple implants".format(
                name=context.public_name
            )
        )


def answer_row(context):
    expected_qc = "PASSED" if context.implant.error_type == "NORMAL" else "FAILED"
    return {
        "public_name": context.public_name,
        "species": context.species,
        "reference_accession": context.record["accession"],
        "tax_classification": context.species,
        "assembler": "Unknown",
        "qc": expected_qc,
        "notes": context.implant.notes,
        "error_type": context.implant.error_type,
        "severity": context.implant.severity,
    }


def sample_sheet_row(context):
    return {
        "sample_name": context.public_name,
        "reference_accession": context.record["accession"],
        "species": context.species,
        "tax_classification": "",
        "r1": os.path.basename(context.r1),
        "r2": os.path.basename(context.r2),
        "long_reads": os.path.basename(context.long_reads),
        "assembler": "",
        "qc": "",
        "notes": "",
    }


def manifest_row(context):
    return {
        "sample_name": context.public_name,
        "reference_accession": context.record["accession"],
        "species": context.species,
        "error_type": context.implant.error_type,
        "severity": context.implant.severity,
        "notes": context.implant.notes,
        "short_read_count": count_reads(context.r1),
        "long_read_count": count_reads(context.long_reads),
    }


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
    contexts = []
    for index, (record, error_type) in enumerate(zip(records, error_plan)):
        context = make_sample_context(record, index, random_seed)
        simulate_base_reads(context, output_dir)
        implant_context(context, error_type, contaminant_records, output_dir)
        contexts.append(context)

    answer_rows = [answer_row(context) for context in contexts]
    sample_rows = [sample_sheet_row(context) for context in contexts]
    manifest_rows = [manifest_row(context) for context in contexts]

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

    summary = Counter(error_plan)
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
