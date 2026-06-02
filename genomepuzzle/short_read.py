"""
Generate short-read assembly datasets with explicit implanted errors.
"""

import csv
import json
import logging
import os
import random
from collections import Counter

from genomepuzzle.create_error import (
    contamination,
    corrupt_read,
    count_reads,
    degrade_quality,
    pass_through,
    subsample_paired_fastq,
    truncate_fastq,
)


SHORT_READ_ERROR_TYPES = [
    "NORMAL",
    "CONTAMINATED",
    "LOW_COVERAGE",
    "POOR_QUALITY",
    "TRUNCATED",
    "CORRUPT",
]

LOW_COVERAGE_FRACTION = 0.15
CONTAMINATION_FRACTION = 30
POOR_QUALITY_MIN = 5
POOR_QUALITY_MAX = 14
TRUNCATED_READ_LENGTH = 35


class ShortReadImplant(object):
    def __init__(self, error_type, severity, notes):
        self.error_type = error_type
        self.severity = severity
        self.notes = notes


class ShortReadSampleContext(object):
    def __init__(self, record, public_name, species, seed, source_dir, implant):
        self.record = record
        self.public_name = public_name
        self.species = species
        self.seed = seed
        self.source_dir = source_dir
        self.source_r1 = os.path.join(source_dir, os.path.basename(record["r1"]))
        self.source_r2 = os.path.join(source_dir, os.path.basename(record["r2"]))
        self.r1 = self.source_r1
        self.r2 = self.source_r2
        self.implant = implant
        self.implant_count = 0


def _write_csv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_records(sample_sheet):
    with open(sample_sheet, "r", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle)]


def build_short_read_error_plan(num_samples, error_proportion, random_seed=42):
    if num_samples < 0:
        raise ValueError("num_samples must be non-negative")
    if not 0 <= error_proportion <= 1:
        raise ValueError("error_proportion must be between 0 and 1")

    plan = ["NORMAL"] * num_samples
    error_count = min(num_samples, int(num_samples * error_proportion))
    if error_count == 0:
        return plan

    rng = random.Random(random_seed)
    indices = list(range(num_samples))
    rng.shuffle(indices)
    implants = [
        SHORT_READ_ERROR_TYPES[1 + (i % (len(SHORT_READ_ERROR_TYPES) - 1))]
        for i in range(error_count)
    ]
    rng.shuffle(implants)
    for index, implant in zip(indices[:error_count], implants):
        plan[index] = implant
    return plan


def make_sample_context(record, index, random_seed, source_dir):
    public_name = "sample{count:02d}".format(count=index + 1)
    species = record.get("SPECIES") or record.get("species") or "Unknown"
    return ShortReadSampleContext(
        record=dict(record),
        public_name=public_name,
        species=species,
        seed=random_seed + index,
        source_dir=source_dir,
        implant=ShortReadImplant("NORMAL", "none", "No implant"),
    )


def _sample_output_paths(output_dir, public_name):
    return (
        os.path.join(output_dir, "{name}_R1.fastq.gz".format(name=public_name)),
        os.path.join(output_dir, "{name}_R2.fastq.gz".format(name=public_name)),
    )


def apply_normal(context, output_dir):
    output_r1, output_r2 = _sample_output_paths(output_dir, context.public_name)
    pass_through(context.source_r1, context.source_r2, output_r1, output_r2)
    context.r1 = output_r1
    context.r2 = output_r2


def apply_low_coverage(context, output_dir):
    output_r1, output_r2 = _sample_output_paths(output_dir, context.public_name)
    subsample_paired_fastq(
        context.source_r1,
        context.source_r2,
        output_r1,
        output_r2,
        subsample_fraction=LOW_COVERAGE_FRACTION,
        random_seed=context.seed,
    )
    context.r1 = output_r1
    context.r2 = output_r2
    context.implant = ShortReadImplant(
        "LOW_COVERAGE",
        round(LOW_COVERAGE_FRACTION, 4),
        "Paired reads subsampled to {:.0f}% of original".format(
            LOW_COVERAGE_FRACTION * 100
        ),
    )
    context.implant_count += 1


def apply_poor_quality(context, output_dir):
    output_r1, output_r2 = _sample_output_paths(output_dir, context.public_name)
    degrade_quality(
        context.source_r1,
        output_r1,
        min_quality=POOR_QUALITY_MIN,
        max_quality=POOR_QUALITY_MAX,
        random_seed=context.seed,
    )
    degrade_quality(
        context.source_r2,
        output_r2,
        min_quality=POOR_QUALITY_MIN,
        max_quality=POOR_QUALITY_MAX,
        random_seed=context.seed,
    )
    context.r1 = output_r1
    context.r2 = output_r2
    context.implant = ShortReadImplant(
        "POOR_QUALITY",
        "{min_q}-{max_q}".format(
            min_q=POOR_QUALITY_MIN, max_q=POOR_QUALITY_MAX
        ),
        "Base qualities degraded across both reads",
    )
    context.implant_count += 1


def apply_truncated(context, output_dir):
    output_r1, output_r2 = _sample_output_paths(output_dir, context.public_name)
    truncate_fastq(context.source_r1, output_r1, truncate_length=TRUNCATED_READ_LENGTH)
    truncate_fastq(context.source_r2, output_r2, truncate_length=TRUNCATED_READ_LENGTH)
    context.r1 = output_r1
    context.r2 = output_r2
    context.implant = ShortReadImplant(
        "TRUNCATED",
        TRUNCATED_READ_LENGTH,
        "Both reads truncated to {length} bases".format(
            length=TRUNCATED_READ_LENGTH
        ),
    )
    context.implant_count += 1


def apply_corrupt(context, output_dir):
    output_r1, output_r2 = _sample_output_paths(output_dir, context.public_name)
    pass_through(context.source_r1, context.source_r2, output_r1, output_r2)
    corrupt_read(output_r1, random_seed=context.seed)
    corrupt_read(output_r2, random_seed=context.seed + 1)
    context.r1 = output_r1
    context.r2 = output_r2
    context.implant = ShortReadImplant(
        "CORRUPT",
        "null-byte",
        "Injected one random null byte into each gzipped FASTQ",
    )
    context.implant_count += 1


def _choose_contaminant(species, contaminant_records, rng):
    compatible = [
        record
        for record in contaminant_records
        if record.get("SPECIES") != species and record.get("ASSEMBLY")
    ]
    return rng.choice(compatible) if compatible else None


def apply_contaminated(context, output_dir, contaminant_records):
    output_r1, output_r2 = _sample_output_paths(output_dir, context.public_name)
    rng = random.Random(context.seed)
    contaminant_record = _choose_contaminant(context.species, contaminant_records, rng)
    if contaminant_record is None:
        logging.warning(
            "No contaminant available for %s, falling back to NORMAL",
            context.public_name,
        )
        apply_normal(context, output_dir)
        return
    contamination(
        context.source_r1,
        context.source_r2,
        output_r1,
        output_r2,
        contaminant_record["ASSEMBLY"],
        output_dir,
        percentage=CONTAMINATION_FRACTION,
        random_seed=context.seed,
    )
    context.r1 = output_r1
    context.r2 = output_r2
    context.implant = ShortReadImplant(
        "CONTAMINATED",
        round(CONTAMINATION_FRACTION / 100.0, 4),
        "Mixed with contaminant assembly {assembly}".format(
            assembly=contaminant_record["ASSEMBLY"]
        ),
    )
    context.implant_count += 1


def implant_context(context, error_type, output_dir, contaminant_records):
    if error_type == "NORMAL":
        apply_normal(context, output_dir)
    elif error_type == "LOW_COVERAGE":
        apply_low_coverage(context, output_dir)
    elif error_type == "POOR_QUALITY":
        apply_poor_quality(context, output_dir)
    elif error_type == "TRUNCATED":
        apply_truncated(context, output_dir)
    elif error_type == "CORRUPT":
        apply_corrupt(context, output_dir)
    elif error_type == "CONTAMINATED":
        apply_contaminated(context, output_dir, contaminant_records)
    else:
        raise ValueError("Unsupported short-read implant type: {0}".format(error_type))
    if context.implant_count > 1:
        raise ValueError(
            "Short-read sample {name} received multiple implants".format(
                name=context.public_name
            )
        )


def answer_row(context):
    qc = "PASSED" if context.implant.error_type == "NORMAL" else "FAILED"
    return {
        "ID": context.public_name,
        "R1": os.path.basename(context.r1),
        "R2": os.path.basename(context.r2),
        "SPECIES": context.species,
        "QC": qc,
        "ERROR": context.implant.error_type,
        "ST": context.record.get("ST", "Unknown"),
        "AMR": context.record.get("AMR", "Unknown"),
        "Notes": context.implant.notes,
    }


def sample_sheet_row(context):
    return {
        "ID": context.public_name,
        "R1": os.path.basename(context.r1),
        "R2": os.path.basename(context.r2),
        "SPECIES": context.species,
        "QC": "Unknown",
        "ERROR": "Unknown",
        "ST": "Unknown",
        "AMR": "Unknown",
        "Notes": "",
    }


def manifest_row(context):
    return {
        "sample_name": context.public_name,
        "species": context.species,
        "source_r1": os.path.basename(context.source_r1),
        "source_r2": os.path.basename(context.source_r2),
        "error_type": context.implant.error_type,
        "severity": context.implant.severity,
        "notes": context.implant.notes,
        "read_count": count_reads(context.r1),
    }


def create_short_read_error_dataset(
    sample_sheet,
    error_proportion,
    contamination_list_file,
    output_dir,
    random_seed=42,
):
    records = read_records(sample_sheet)
    if not records:
        raise ValueError("sample_sheet contains no records")

    os.makedirs(output_dir, exist_ok=True)
    source_dir = os.path.dirname(sample_sheet) or "."
    with open(contamination_list_file, "r", encoding="utf-8") as handle:
        contaminant_records = [row for row in csv.DictReader(handle)]

    error_plan = build_short_read_error_plan(
        len(records), error_proportion=error_proportion, random_seed=random_seed
    )
    contexts = []
    for index, (record, error_type) in enumerate(zip(records, error_plan)):
        context = make_sample_context(record, index, random_seed, source_dir)
        implant_context(context, error_type, output_dir, contaminant_records)
        contexts.append(context)

    answer_rows = [answer_row(context) for context in contexts]
    sample_rows = [sample_sheet_row(context) for context in contexts]
    manifest_rows = [manifest_row(context) for context in contexts]

    _write_csv(
        os.path.join(output_dir, "answer_sheet.csv"),
        answer_rows,
        ["ID", "R1", "R2", "SPECIES", "QC", "ERROR", "ST", "AMR", "Notes"],
    )
    _write_csv(
        os.path.join(output_dir, "sample_sheet.csv"),
        sample_rows,
        ["ID", "R1", "R2", "SPECIES", "QC", "ERROR", "ST", "AMR", "Notes"],
    )
    _write_csv(
        os.path.join(output_dir, "implant_manifest.csv"),
        manifest_rows,
        [
            "sample_name",
            "species",
            "source_r1",
            "source_r2",
            "error_type",
            "severity",
            "notes",
            "read_count",
        ],
    )
    with open(
        os.path.join(output_dir, "implant_manifest.json"), "w", encoding="utf-8"
    ) as handle:
        json.dump(
            {
                "random_seed": random_seed,
                "error_proportion": error_proportion,
                "summary": dict(Counter(error_plan)),
                "samples": manifest_rows,
            },
            handle,
            indent=2,
        )
