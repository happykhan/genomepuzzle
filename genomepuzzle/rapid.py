"""
This module provides functions to simulate sequencing reads, run genome assembly tools,
nd process genome data.

Functions:
    run_spades(record, output_dir):

    run_badread(accession, assembly_path, output_dir):

    run_flye(record, genomelen, long_read, short_r1, short_r2, output_dir):

    rapid(output_dir, samplelist, random_seed=42):
"""

import logging
import csv
import os
from genomepuzzle.simulate_reads import run_art, fetch_assembly, cleanup_output_dir
from genomepuzzle.assembly_stats import calculate_assembly_stats


def run_spades(record, output_dir):
    """
    Runs the SPAdes genome assembler using Docker and processes the output.

    Args:
        record (dict): A dictionary containing information about the genome record.
                       Must include the 'accession' key.
        output_dir (str): The directory where the output files will be stored.

    Returns:
        dict: The updated record dictionary with additional keys for SPAdes assembly statistics
              and the path to the final assembly file.

    Raises:
        OSError: If there is an issue with running the Docker command or file operations.

    Notes:
        - This function assumes that Docker is installed and properly configured on the system.
        - The SPAdes Docker image used is `quay.io/biocontainers/spades:3.11.0--py36_0`.
        - The function calculates assembly statistics using the `calculate_assembly_stats` function.
        - The function moves the final assembly file to the output directory and deletes the intermediate SPAdes output directory.
    """
    abs_output_dir = os.path.abspath(output_dir)
    os.system(
        f"docker run --platform linux/x86_64 -v {abs_output_dir}:/data "
        f"quay.io/biocontainers/spades:3.11.0--py36_0 spades.py "
        f"-1 /data/{record['accession']}_R1.fastq.gz "
        f"-2 /data/{record['accession']}_R2.fastq.gz "
        f"-o /data/{record['accession']}_spades"
    )
    output_contigs = f"{output_dir}/{record['accession']}_spades/scaffolds.fasta"
    if os.path.exists(output_contigs):
        stats = calculate_assembly_stats(output_contigs)
        record["spades_total_bases"] = stats["Total assembly size"]
        record["spades_number_of_contigs"] = stats["Number of contigs"]
        record["spades_n50"] = stats["N50"]
        record["spades_gc_content"] = stats["GC content (%)"]
    else:
        logging.error("no spades output")
    # copy the assembly to the output dir
    os.system(f"mv {output_contigs} {output_dir}/{record['accession']}_spades.fasta")
    # delete spades output dir
    os.system(f"rm -r {output_dir}/{record['accession']}_spades")
    record["spades_assembly"] = f"{record['accession']}_scaffolds.fasta"
    return record


def run_badread(accession, assembly_path, output_dir):
    """
    Simulates long reads using the badread tool and saves the output in the specified directory.

    Parameters:
    accession (str): The accession identifier for the output file.
    assembly_path (str): The path to the assembly file to be used as a reference.
    output_dir (str): The directory where the output files will be saved.

    Returns:
    str: The path to the generated long reads file in gzip format.
    """
    abs_output_dir = os.path.abspath(output_dir)
    assembly_file = os.path.basename(assembly_path)
    # copy the assembly to the output dir
    os.system(f"cp {assembly_path} {output_dir}")
    os.system(
        f"docker run --platform linux/x86_64 -v {abs_output_dir}:/data "
        f"quay.io/biocontainers/badread:0.4.1--pyhdfd78af_0 badread simulate "
        f"--reference /data/{assembly_file} --quantity 20x | gzip > {abs_output_dir}/{accession}_long.fastq.gz"
    )
    long_reads_path = f"{output_dir}/{accession}_long.fastq.gz"
    os.remove(f"{output_dir}/{assembly_file}")
    return long_reads_path


def run_flye(record, genomelen, long_read, short_r1, short_r2, output_dir):
    """
    Runs the Flye assembler using Docker and processes the output assembly.

    Args:
        record (dict): A dictionary containing metadata for the genome assembly,
        including the accession number.
        genomelen (int): The estimated genome length.
        long_read (str): The path to the long-read sequencing data file.
        short_r1 (str): The path to the first short-read sequencing data file (R1).
        short_r2 (str): The path to the second short-read sequencing data file (R2).
        output_dir (str): The directory where the output files will be stored.

    Returns:
        dict: The updated record dictionary with additional Flye assembly statistics and file paths.
    """
    accession = record["accession"]
    abs_output_dir = os.path.abspath(output_dir)
    short_r1 = os.path.basename(short_r1)
    short_r2 = os.path.basename(short_r2)
    long_read_file = os.path.basename(long_read)
    os.system(
        f"docker run --platform linux/x86_64 -v {abs_output_dir}:/data "
        f"quay.io/biocontainers/flye:2.3.6--py27ha92aebf_3 flye --nano-raw "
        f"/data/{long_read_file} -o /data/{accession}_flye -g {genomelen}"
    )
    # copy the assembly to the output dir
    output_assembly = f"{output_dir}/{accession}_longassembly.fasta"
    os.system(f"mv {output_dir}/{accession}_flye/scaffolds.fasta {output_assembly}")
    os.system(f"rm -r {output_dir}/{accession}_flye")
    stats = calculate_assembly_stats(output_assembly)
    record["flye_total_bases"] = stats["Total assembly size"]
    record["flye_number_of_contigs"] = stats["Number of contigs"]
    record["flye_n50"] = stats["N50"]
    record["flye_gc_content"] = stats["GC content (%)"]
    record["flye_assembly"] = os.path.basename(output_assembly)
    return record


def rapid(output_dir, samplelist, random_seed=42):
    """
    Simulates sequencing reads, runs assembly tools, and processes genome data.

    Parameters:
    output_dir (str): The directory where output files will be saved.
    samplelist (str): Path to the CSV file containing sample information.
    random_seed (int, optional): Seed for random number generation. Default is 42.

    Raises:
    RuntimeError: If Docker is not installed.

    The function performs the following steps:
    1. Checks if Docker is installed.
    2. Reads the sample list from the provided CSV file.
    3. Creates the output directory if it does not exist.
    4. Fetches assembly data for each sample.
    5. Copies the assembly files to the output directory.
    6. Simulates long reads using badread.
    7. Simulates paired-end reads using ART.
    8. Runs Flye assembler.
    9. Runs SPAdes assembler to get expected assembly metrics.
    10. Writes the processed sample information to a new CSV file in the output directory.
    11. Cleans up the output directory.

    Note:
    - The function assumes that the necessary tools (badread, ART, Flye, SPAdes) are installed and available in the system's PATH.
    - The sample CSV file should contain the following columns: 'accession', 'assemblystats_totalsequencelength'.
    """
    # check docker is installed
    if os.system("docker --version") != 0:
        raise RuntimeError("Docker is not installed")
    all_records = [x for x in csv.DictReader(open(samplelist, encoding="utf-8"))]
    os.makedirs(output_dir, exist_ok=True)
    assembly_accessions = [record["accession"] for record in all_records]
    # fetch assembly data
    fetch_assembly(assembly_accessions, output_dir)
    for record in all_records:
        assembly_dir = os.path.join(
            output_dir, "ncbi_dataset", "data", record["accession"]
        )
        assembly_path = [
            os.path.join(assembly_dir, x)
            for x in os.listdir(assembly_dir)
            if x.endswith(".fna")
        ][0]
        # copy the assembly to the output dir
        os.system(
            f"cp {assembly_path} {output_dir}/{record['accession']}_original.fasta"
        )
        record["original_assembly"] = f"{record['accession']}_original.fasta"
        output_r1 = f"{output_dir}/{record['accession']}_R1.fastq.gz"
        output_r2 = f"{output_dir}/{record['accession']}_R2.fastq.gz"
        sample = {
            "public_name": record["accession"],
            "platform": "HS25",
            "read_length": 150,
            "coverage": 30,
            "fragment_length": 300,
            "standard_deviation": 50,
            "random_seed": random_seed,
        }
        long_read_path = f"{output_dir}/{record['accession']}_long.fastq.gz"
        # simulate long reads with badread
        long_read_path = run_badread(record["accession"], assembly_path, output_dir)
        # simulate reads with art
        run_art(sample, output_dir, assembly_path, output_r1, output_r2)
        # run flye
        record = run_flye(
            record,
            record["assemblystats_totalsequencelength"],
            long_read_path,
            output_r1,
            output_r2,
            output_dir,
        )
        # assembly with spades to get expected assembly metrics
        record = run_spades(record, output_dir)

    # write all_records to a new file in the output_dir
    with open(f"{output_dir}/sample_sheet.csv", "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_records[0].keys())
        writer.writeheader()
        for record in all_records:
            writer.writerow(record)

    cleanup_output_dir(output_dir)
