"""
This module provides various functions to manipulate and introduce errors into sequencing reads.
The functions include degrading quality, truncating reads, subsampling reads,
introducing contamination, and corrupting reads. The module also includes utility functions
for counting reads and passing through files.

Functions:
- degrade_quality: Degrade base quality scores randomly, with poorer quality towards the end of each read.
- truncate_fastq: Truncate reads in a FASTQ file to a specific length.
- subsample_paired_fastq: Subsample paired-end FASTQ files to reduce coverage.
- subsample_paired_read_by_count: Subsample a specific number of reads from paired-end FASTQ files.
- pass_through: Pass through read files to output directory.
- get_contaminated_read_example: Reads a contamination list file and returns a list of contaminants that
    are not in the given species list.
- corrupt_read: Corrupts a read file by writing a null byte at a random position.
- count_reads: Count the number of reads in a FASTQ file.
- contamination: Create contaminated reads by simulating from a contaminant assembly and mixing them with
    input reads.
- introduce_errors: Introduce errors into sequencing reads based on specified error proportions.
"""

import random
import csv
import logging
import os
import shutil
import gzip
from genomepuzzle.simulate_reads import run_art, fetch_assembly, cleanup_output_dir
import subprocess


def degrade_quality(
    input_fastq, output_fastq, min_quality=5, max_quality=30, random_seed=42
):
    """
    Degrade base quality scores randomly, with poorer quality towards the end of each read.

    Parameters:
    - input_fastq: str, path to the input FASTQ file.
    - output_fastq: str, path to the output FASTQ file.
    - min_quality: int, minimum Phred score for base quality at the end of the read.
    - max_quality: int, maximum Phred score for base quality at the start of the read.
    """
    random.seed(random_seed)
    if not 0 <= min_quality <= max_quality <= 93:  # Ensure valid Phred+33 range
        raise ValueError("Quality scores must be in range 0-93.")

    with (
        gzip.open(input_fastq, "rt") as infile,
        gzip.open(output_fastq, "wt") as outfile,
    ):
        while True:
            # Read one FASTQ record (4 lines)
            header = infile.readline().strip()  # Line 1: Sequence ID
            sequence = infile.readline().strip()  # Line 2: Nucleotide sequence
            plus = infile.readline().strip()  # Line 3: +
            # quality = infile.readline().strip()  # Line 4: Quality scores

            # Break if end of file
            if not header:
                break

            read_length = len(sequence)
            degraded_quality = []

            # Generate progressively poorer quality scores
            for i in range(read_length):
                # Calculate a score range that worsens towards the end of the read
                current_max_quality = max_quality - int(
                    (i / read_length) * (max_quality - min_quality)
                )
                random_quality = random.randint(min_quality, current_max_quality)
                degraded_quality.append(
                    chr(random_quality + 33)
                )  # Convert Phred score to ASCII

            # Join the degraded quality scores
            new_quality = "".join(degraded_quality)

            # Write the modified FASTQ record
            outfile.write(f"{header}\n{sequence}\n{plus}\n{new_quality}\n")

    logging.info("Quality degradation complete. Output saved to: %s", output_fastq)


def truncate_fastq(input_fastq, output_fastq, truncate_length=75):
    """
    Truncate reads in a FASTQ file to a specific length.

    Parameters:
    - sample: dict, sample information containing 'r1' and 'r2' keys for read files.
    - output_fastq: str, path to the output truncated FASTQ file.
    - truncate_length: int, length to truncate each read to.
    """

    with (
        open(input_fastq, "r", encoding="utf-8") as infile,
        open(output_fastq, "w", encoding="utf-8") as outfile,
    ):
        while True:
            # Read one complete FASTQ record (4 lines)
            header = infile.readline().strip()  # Line 1: Sequence ID
            sequence = infile.readline().strip()  # Line 2: Nucleotide sequence
            plus = infile.readline().strip()  # Line 3: +
            quality = infile.readline().strip()  # Line 4: Quality scores

            # Break if we reach the end of the file
            if not header:
                break

            # Truncate sequence and quality to the specified length
            truncated_sequence = sequence[:truncate_length]
            truncated_quality = quality[:truncate_length]

            # Write the truncated record to the output file
            outfile.write(
                f"{header}\n{truncated_sequence}\n{plus}\n{truncated_quality}\n"
            )

    logging.info("Truncation complete. Output saved to: %s", output_fastq)


def subsample_paired_fastq(
    input_r1,
    input_r2,
    output_r1,
    output_r2,
    subsample_fraction=0.1,
    random_seed=42,
    seqtk=True,
):
    """
    Subsample paired-end FASTQ files to reduce coverage.

    Parameters:
    - sample: dict, sample information containing 'r1' and 'r2' keys for read files.
    - output_r1: str, output subsampled FASTQ file for read 1.
    - output_r2: str, output subsampled FASTQ file for read 2.
    - subsample_fraction: float, fraction of reads to retain (default = 0.1 for 10%).
    """
    random.seed(random_seed)
    assert 0 < subsample_fraction <= 1, "Subsample fraction must be between 0 and 1."
    if seqtk:
        logging.info("Subsampling using seqtk (r1)...")
        subprocess.run(
            f"gunzip -c {input_r1} | bin/seqtk sample -s {random_seed} - {subsample_fraction} | gzip > {output_r1}",
            shell=True,
            check=True,
        )
        logging.info("Subsampling using seqtk (r2)...")
        subprocess.run(
            f"gunzip -c {input_r2} | bin/seqtk sample -s {random_seed} - {subsample_fraction} | gzip > {output_r2}",
            shell=True,
            check=True,
        )
    else:
        with (
            gzip.open(input_r1, "rt") as r1,
            gzip.open(input_r2, "rt") as r2,
            gzip.open(output_r1, "wt") as out_r1,
            gzip.open(output_r2, "wt") as out_r2,
        ):
            while True:
                # Read 4 lines for R1 (one full FASTQ record)
                r1_record = [r1.readline().strip() for _ in range(4)]
                r2_record = [r2.readline().strip() for _ in range(4)]

                # Stop at end of file
                if not r1_record[0] or not r2_record[0]:
                    break

                # Randomly keep the paired reads based on the subsample fraction
                if random.random() < subsample_fraction:
                    out_r1.write("\n".join(r1_record) + "\n")
                    out_r2.write("\n".join(r2_record) + "\n")

    logging.info("Subsampling complete. Output files saved as:")
    logging.info("  %s", output_r1)
    logging.info("  %s", output_r2)


def subsample_paired_read_by_count(
    input_r1, input_r2, output_r1, output_r2, num_reads=1000, random_seed=42, seqtk=True
):
    """
    Subsample a specific number of reads from paired-end FASTQ files.

    Parameters:
    - input_r1: str, path to the input R1 FASTQ file.
    - input_r2: str, path to the input R2 FASTQ file.
    - output_r1: str, path to the output subsampled R1 FASTQ file.
    - output_r2: str, path to the output subsampled R2 FASTQ file.
    - num_reads: int, number of reads to subsample.
    - random_seed: int, seed for the random number generator.
    """
    random.seed(random_seed)
    if seqtk:
        logging.info("Subsampling by count using seqtk (r1)...")
        subprocess.run(
            f"gunzip -c {input_r1} | bin/seqtk sample -s {random_seed} - {num_reads} | gzip > {output_r1}",
            shell=True,
            check=True,
        )
        logging.info("Subsampling by count using seqtk (r2)...")
        subprocess.run(
            f"gunzip -c {input_r2} | bin/seqtk sample -s {random_seed} - {num_reads} | gzip > {output_r2}",
            shell=True,
            check=True,
        )
    else:
        with (
            gzip.open(input_r1, "rt") as r1,
            gzip.open(input_r2, "rt") as r2,
            gzip.open(output_r1, "wt") as out_r1,
            gzip.open(output_r2, "wt") as out_r2,
        ):
            reads_r1 = []
            reads_r2 = []

            while True:
                r1_record = [r1.readline().strip() for _ in range(4)]
                r2_record = [r2.readline().strip() for _ in range(4)]

                if not r1_record[0] or not r2_record[0]:
                    break

                reads_r1.append(r1_record)
                reads_r2.append(r2_record)

            selected_indices = random.sample(range(len(reads_r1)), int(num_reads))

            for idx in selected_indices:
                out_r1.write("\n".join(reads_r1[idx]) + "\n")
                out_r2.write("\n".join(reads_r2[idx]) + "\n")

    logging.info("Subsampling by count complete. Output files saved as:")
    logging.info("  %s", output_r1)
    logging.info("  %s", output_r2)


def pass_through(r1_path, r2_path, r1_output, r2_output):
    """
    Pass through read files to output directory.

    Parameters:
    - r1_path: str, path to the read 1 file.
    - r2_path: str, path to the read 2 file.
    - r1_output: str, path to the output read 1 file.
    - r2_output: str, path to the output read 2 file.

    Returns:
    - None
    """
    shutil.copy(r1_path, r1_output)
    shutil.copy(r2_path, r2_output)


def get_contaminated_read_example(species_list, contamination_list_file):
    """
    Reads a contamination list file and returns a list of contaminants that
    are not in the given species list.

    Args:
        species_list (list): A list of species names to exclude from the contamination list.
        contamination_list_file (str): The file path to the contamination list CSV file.

    Returns:
        list: A list of dictionaries representing the contaminants that are
        not in the species list and have an assembly.
    """
    contaminant_list = []
    for x in csv.DictReader(open(contamination_list_file, "r", encoding="utf-8")):
        if x["SPECIES"] not in species_list and x.get("ASSEMBLY"):
            contaminant_list.append(x)
    return contaminant_list


def corrupt_read(sample, output_dir, random_seed=42):
    """
    Corrupts a read file by writing a null byte at a random position.

    Args:
        sample (dict): A dictionary containing sample information,
                       where "r1" is the key for the read file path.
        output_dir (str): The directory where the corrupted file will be saved.
        random_seed (int, optional): The seed for the random number generator.
                                     Defaults to 42.

    Returns:
        dict: The original sample dictionary.
    """
    random.seed(random_seed)
    file_to_corrupt = os.path.join(output_dir, os.path.basename(sample["r1"]))
    with open(file_to_corrupt, "r+b") as f:
        f.seek(random.randint(0, os.path.getsize(file_to_corrupt) - 1))
        f.write(b"\x00")
    return sample


def count_reads(fastq_file):
    """
    Count the number of reads in a FASTQ file.

    Parameters:
    - fastq_file: str, path to the FASTQ file.

    Returns:
    - int, number of reads in the FASTQ file.
    """
    count = 0
    with gzip.open(fastq_file, "rt") as f:
        for _ in f:
            count += 1
    return count // 4


def contamination(
    input_r1,
    input_r2,
    output_r1,
    output_r2,
    contaminant_assembly,
    output_dir,
    percentage=50,
    random_seed=42,
):
    """
    Create contaminated reads by simulating from a contaminant assembly and mixing
      them with input reads.

    Parameters:
    input_r1 (str): Path to the input R1 fastq file.
    input_r2 (str): Path to the input R2 fastq file.
    output_r1 (str): Path to the output R1 fastq file.
    output_r2 (str): Path to the output R2 fastq file.
    contaminant_assembly (str): Path to the contaminant assembly file.
    output_dir (str): Directory to store the temporary contaminant files.
    percentage (int, optional): Percentage of contaminant reads to mix with input reads.
        Default is 50.
    random_seed (int, optional): Seed for random number generator. Default is 42.

    Returns:
    None
    """
    fetch_assembly([contaminant_assembly], output_dir)
    ass_dir = os.path.join(output_dir, "ncbi_dataset", "data")
    ass_dir = [
        os.path.join(ass_dir, x)
        for x in os.listdir(ass_dir)
        if x.startswith(contaminant_assembly)
    ][0]

    contaminant_assembly_path = [
        os.path.join(ass_dir, x)
        for x in os.listdir(ass_dir)
        if x.startswith(contaminant_assembly)
    ][0]
    # output_dir/ncbi_dataset/data/dataset_catalog.json

    sample = {
        "public_name": "contaminant",
        "platform": "HS25",
        "read_length": 150,
        "coverage": 10,
        "fragment_length": 300,
        "standard_deviation": 30,
        "random_seed": random_seed,
        "SPECIES": "contaminant",
        "ASSEMBLY": contaminant_assembly_path,
        "SAMPLE_NAME": "contaminant",
    }
    contaminant_output_r1 = os.path.join(output_dir, "contaminant_r1.fastq.gz")
    contaminant_output_r2 = os.path.join(output_dir, "contaminant_r2.fastq.gz")
    run_art(
        sample,
        output_dir,
        contaminant_assembly_path,
        contaminant_output_r1,
        contaminant_output_r2,
    )
    # mix the reads with the contaminant simulated reads
    # how many reads are in input_r1 ?
    num_reads = int(round(count_reads(input_r1)))
    adj_num_reads = (percentage / 100) * num_reads
    logging.info("Number of reads in input_r1: %d", num_reads)
    logging.info("Number of reads to mix in: %d", adj_num_reads)
    # pick subsample adj_num_reads from contaminant_output_r1 and contaminant_output_r2
    subsample_contaminant_output_r1 = os.path.join(
        output_dir, "sub_contaminant_r1.fastq.gz"
    )
    subsample_contaminant_output_r2 = os.path.join(
        output_dir, "sub_contaminant_r2.fastq.gz"
    )
    logging.info(
        "Subsampling %d reads from contaminant reads (%d percent)",
        adj_num_reads,
        percentage,
    )
    subsample_paired_read_by_count(
        input_r1,
        input_r2,
        subsample_contaminant_output_r1,
        subsample_contaminant_output_r2,
        adj_num_reads,
        random_seed,
    )
    logging.info("Subsampling %d percent from real reads", (100 - percentage))
    subsample_paired_fastq(
        input_r1,
        input_r2,
        output_r1,
        output_r2,
        (100 - percentage) / 100,
        random_seed=random_seed,
    )
    logging.info("Appending the contaminant reads to the output files - r1")
    subprocess.run(
        f"cat {subsample_contaminant_output_r1} >> {output_r1}", shell=True, check=True
    )
    logging.info("Appending the contaminant reads to the output files - r2")
    subprocess.run(
        f"cat {subsample_contaminant_output_r2} >> {output_r2}", shell=True, check=True
    )
    logging.info("Removing temporary files...")
    os.remove(contaminant_output_r1)
    os.remove(subsample_contaminant_output_r1)
    os.remove(contaminant_output_r2)
    os.remove(subsample_contaminant_output_r2)
    cleanup_output_dir(output_dir)


def introduce_errors(
    sample_sheet, error_proportion, contamination_list_file, output_dir, random_seed=42
):
    """
    Introduce errors into sequencing reads based on specified error proportions.

    Parameters:
    sample_sheet (str): Path to the CSV file containing sample information.
    error_proportion (float): Proportion of samples to introduce errors into.
    contamination_list_file (str): Path to the file containing contamination details.
    output_dir (str): Directory to save the output files with introduced errors.
    random_seed (int, optional): Seed for random number generator to ensure
                                 reproducibility. Default is 42.

    Returns:
    None

    This function reads the sample information from the provided sample_sheet,
    introduces various types of errors (e.g., truncated reads, corrupted reads,
    contamination, poor quality, low coverage) into a proportion of the samples,
    and writes the modified samples to the output directory. The types of errors
    and their proportions are determined randomly based on the provided
    error_proportion. The function also updates the sample sheet with details of
    the introduced errors and saves it to the output directory.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    random.seed(random_seed)
    types_of_errors = [
        "TRUNCATED",
        "CORRUPT",
        "CONTAMINATED",
        "POOR_QUALITY",
        "LOW_COVERAGE",
    ]

    full_sample_list = []
    with open(sample_sheet, "r", encoding="utf-8") as f:
        full_sample_list = list(csv.DictReader(f))
    # pick with samples to alter given the error_proportion
    # get list of species in sample_sheet
    species_list = set([sample["SPECIES"] for sample in full_sample_list])
    contaminant_list = get_contaminated_read_example(
        species_list, contamination_list_file
    )
    if len(contaminant_list) == 0:
        raise ValueError("No contaminants found in the contamination list.")
    error_samples = random.sample(
        full_sample_list, int(len(full_sample_list) * error_proportion)
    )
    for sample in error_samples:
        if random.random() < 0.4:
            sample["ERROR"] = "CONTAMINATED"
        elif random.random() < 0.2:
            sample["ERROR"] = "LOW_COVERAGE"
        else:
            sample["ERROR"] = random.choice(types_of_errors)
        sample["QC"] = "FAILED"
        sample["Notes"] = "Error introduced"
    used_contaminants = []
    for i, sample in enumerate(full_sample_list):
        count = i+1
        zfill_count = str(count).zfill(2)
        sample['public_name'] = f"sample{zfill_count}"
        input_dir = os.path.dirname(sample_sheet)
        input_r1 = os.path.join(input_dir, os.path.basename(sample["r1"]))
        input_r2 = os.path.join(input_dir, os.path.basename(sample["r2"]))
        output_r1 = os.path.join(output_dir, sample['public_name'] + "_R1.fastq.gz")
        output_r2 = os.path.join(output_dir, sample['public_name'] + "_R2.fastq.gz")
        error_type = sample["ERROR"]
        if error_type == "NONE":
            pass_through(input_r1, input_r2, output_r1, output_r2)
            sample["Notes"] = "No changes."
        elif error_type == "TRUNCATED":
            truncate_length = random.randint(5, 100)
            file_to_truncate = input_r1
            if random.choice([True, False]):
                file_to_truncate = input_r1
                sample["Notes"] = f"Truncated read 1 to {truncate_length} bp"
            else:
                sample["Notes"] = f"Truncated read 2 to {truncate_length} bp"
                file_to_truncate = input_r2
            truncate_fastq(file_to_truncate, output_dir, truncate_length)
        elif error_type == "CORRUPT":
            corruption = random.choice(["r1", "r2", "both"])
            if corruption == "both":
                corrupt_read(sample["r1"], output_dir)
                corrupt_read(sample["r2"], output_dir)
                sample["Notes"] = "Corrupted both reads"
            elif corruption == "r1":
                corrupt_read(sample["r1"], output_dir)
                sample["Notes"] = "Corrupted read 1"
            elif corruption == "r2":
                corrupt_read(sample["r2"], output_dir)
                sample["Notes"] = "Corrupted read 2"
        elif error_type == "CONTAMINATED":
            contaminant = random.choice(contaminant_list)
            used_contaminants.append(contaminant)
            contaminant_list.remove(contaminant)
            if len(contaminant_list) == 0:
                contaminant_list = used_contaminants
            contaminant_assembly_file = contaminant["ASSEMBLY"]
            percentage = random.randint(20, 80)
            contamination(
                input_r1,
                input_r2,
                output_r1,
                output_r2,
                contaminant_assembly_file,
                output_dir,
                percentage,
                random_seed,
            )
            sample["Notes"] = (
                f"Contaminated with {contaminant.get('SAMPLE_NAME')} at {percentage}%"
            )
        elif error_type == "POOR_QUALITY":
            min_quality = random.randint(5, 9)
            max_quality = random.randint(10, 20)
            degrade_quality(
                input_r1,
                output_r1,
                min_quality=min_quality,
                max_quality=max_quality,
                random_seed=random_seed,
            )
            degrade_quality(
                input_r2,
                output_r2,
                min_quality=min_quality,
                max_quality=max_quality,
                random_seed=random_seed,
            )
            sample["Notes"] += f" Quality degraded to {min_quality}-{max_quality}"
        elif error_type == "LOW_COVERAGE":
            subsample_fraction = random.uniform(0.1, 0.4)
            sample["Notes"] = f"Low coverage at {round(subsample_fraction, 2)}"
            subsample_paired_fastq(
                input_r1,
                input_r2,
                output_r1,
                output_r2,
                subsample_fraction,
                random_seed=random_seed,
            )
        else:
            logging.error("Error type %s not implemented", error_type)
            raise ValueError(f"Error type {error_type} not implemented")
        # write sample sheet
        sample['r1'] = os.path.basename(output_r1)
        sample['r2'] = os.path.basename(output_r2)
    write_answer_sheet(output_dir, full_sample_list)
    write_sample_sheet(output_dir, full_sample_list)

def write_answer_sheet(output_dir, full_sample_list):
    sample_sheet = os.path.join(output_dir, "answer_sheet.csv")
    with open(sample_sheet, "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=full_sample_list[0].keys())
        writer.writeheader()
        for sample in full_sample_list:
            writer.writerow(sample)

def write_sample_sheet(output_dir, full_sample_list):
    # Write a version of sample_sheet with only columns public_name, r1, r2, species
    filtered_sample_sheet = os.path.join(output_dir, "sample_sheet.csv")
    with open(filtered_sample_sheet, "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["ID", "R1", "R2", "SPECIES", "QC", "ERROR", "ST", "AMR", "Notes"])
        writer.writeheader()
        for sample in full_sample_list:
            writer.writerow({
                "ID": sample["public_name"],
                "R1": sample["r1"],
                "R2": sample["r2"],
                "SPECIES": sample["SPECIES"],
                "QC": 'Unknown',
                "ERROR": 'Unknown',
                "ST": 'Unknown',
                "AMR": 'Unknown',
                'Notes':''
            })