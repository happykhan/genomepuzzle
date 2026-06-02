"""
Rapid long-read and assembly dataset generation.
"""

import csv
import os
import shutil

from genomepuzzle.runtime import require_tool
from genomepuzzle.sample_util import run_badread, run_flye, run_spades
from genomepuzzle.simulate_reads import cleanup_output_dir, fetch_assembly, run_art


def rapid(output_dir, samplelist, random_seed=42):
    """
    Generate a compact assembly benchmark dataset from reference accessions.
    """
    require_tool("docker")
    with open(samplelist, encoding="utf-8") as handle:
        all_records = [row for row in csv.DictReader(handle)]
    os.makedirs(output_dir, exist_ok=True)
    assembly_accessions = [record["accession"] for record in all_records]
    fetch_assembly(assembly_accessions, output_dir)

    for record in all_records:
        assembly_dir = os.path.join(output_dir, "ncbi_dataset", "data", record["accession"])
        assembly_path = [
            os.path.join(assembly_dir, filename)
            for filename in os.listdir(assembly_dir)
            if filename.endswith(".fna")
        ][0]
        original_assembly = os.path.join(output_dir, "{0}_original.fasta".format(record["accession"]))
        shutil.copy2(assembly_path, original_assembly)
        record["original_assembly"] = os.path.basename(original_assembly)

        output_r1 = os.path.join(output_dir, "{0}_R1.fastq.gz".format(record["accession"]))
        output_r2 = os.path.join(output_dir, "{0}_R2.fastq.gz".format(record["accession"]))
        sample = {
            "public_name": record["accession"],
            "platform": "HS25",
            "read_length": 150,
            "coverage": 30,
            "fragment_length": 300,
            "standard_deviation": 50,
            "random_seed": random_seed,
        }
        long_read_path = run_badread(record["accession"], assembly_path, output_dir)["long_reads_path"]
        run_art(sample, output_dir, assembly_path, output_r1, output_r2)
        spades_record = run_spades(output_r1, output_r2, output_dir, record["accession"])
        flye_record = run_flye(
            record["accession"],
            record["assemblystats_totalsequencelength"],
            long_read_path,
            output_dir,
        )
        record.update(spades_record)
        record.update(flye_record)

    with open(os.path.join(output_dir, "sample_sheet.csv"), "w", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_records[0].keys())
        writer.writeheader()
        writer.writerows(all_records)

    cleanup_output_dir(output_dir)
