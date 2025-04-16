"""
This script provides command-line interface for generating the testing dataset for EQA
"""

import argparse
import logging
import sys
import os
import csv
import random
import shutil
import hashlib

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from genomepuzzle.sample.kleborate_sample import KleborateSample
from genomepuzzle.sample.assembly_sample import AssemblySample
from genomepuzzle.sample.contaminated_sample import ContaminatedSample
from genomepuzzle.sample.low_coverage_sample import LowCoverageSample
from genomepuzzle.sample.degrade_sample import DegradedSample
from genomepuzzle.sample.corrupted_sample import CorruptedSample
from genomepuzzle.sample.truncated_sample import TruncatedSample
from genomepuzzle.run_tree2read import create_config_file, run_docker

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def create_kleborate_dataset(
    sample_dict, outputdir, number_of_samples=10, random_seed=42
):
    """
    Create a Kleborate dataset by selecting random samples from a given sample dictionary.
    This function groups samples by their sequence type (ST) and randomly selects one sample
    from each ST group. It picks up to a specified number of samples (default is 10) and
    returns a list of KleborateSample objects.
    Args:
        sample_dict (list of dict): A list of dictionaries where each dictionary represents
                                    a sample with at least an 'st' key for sequence type.
        outputdir (str): The directory where the output files will be stored.
        number_of_samples (int, optional): The maximum number of samples to select. Default is 10.
        random_seed (int, optional): The seed for the random number generator to ensure
            reproducibility. Default is 42.
    Returns:
        list of KleborateSample: A list of KleborateSample objects created from the
        selected samples.
    """
    # Group by ST (st column) and pick one random sample from each ST. Pick up to 10 samples.
    # Convert sample_dict to a list of dictionaries
    kleborate_samples = []
    st_dict = {}
    for sample in sample_dict:
        st = sample["st"]
        if st not in st_dict:
            st_dict[st] = []
        st_dict[st].append(sample)

    # Pick one random sample from each ST. Pick up to 10 samples.
    # pick 10 STs at random
    random.seed(random_seed)
    st_keys = list(st_dict.keys())
    random.shuffle(st_keys)
    st_keys = st_keys[:number_of_samples]
    for st in st_keys:
        samples = st_dict[st]
        sample = random.choice(samples)
        kleborate_samples.append(sample)
    final_samples = [KleborateSample(sample, outputdir) for sample in kleborate_samples]
    return final_samples


def make_typing_output_files(kleb_sample_list, output_dir, test_dir):
    # Create answer sheet
    # for each sample, get the ST and resistance profile and
    #  o-antigen and capsule type; write the answers to a file
    # make kleborate output directory
    os.makedirs(os.path.join(output_dir, test_dir), exist_ok=True)
    fields = []
    all_results = []
    short_fields = []
    short_results = []
    for sample in kleb_sample_list:
        sample.move_assembly_anon(os.path.join(output_dir, test_dir))
        kleborate_output = sample.get_results()
        kleborate_output.update(sample.to_dict())
        fields = list(kleborate_output.keys())
        short_fields = list(sample.to_short_dict().keys())
        all_results.append(kleborate_output)
        short_results.append(sample.to_short_dict(public_name=True))
        if not os.path.exists(sample.assembly_file):
            logging.error("Assembly file %s not found.", sample.assembly_file)
    with open(
        os.path.join(output_dir, test_dir, "answer_sheet.csv"),
        mode="w",
        encoding="utf-8",
    ) as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_results)
    # Create a blank answer sheet for the user to fill in.
    with open(
        os.path.join(output_dir, test_dir, "sample_sheet.csv"),
        mode="w",
        encoding="utf-8",
    ) as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_results)
    with open(
        os.path.join(output_dir, test_dir, "sample_sheet.csv"),
        mode="w",
        encoding="utf-8",
    ) as outfile:
        writer = csv.DictWriter(outfile, fieldnames=short_fields)
        writer.writeheader()
        writer.writerows(short_results)


def create_outbreak_analysis_dataset(ref, tree_file, outbreak_dir):
    print(f"Selected tree file: {tree_file}")
    if not os.path.exists(outbreak_dir):
        os.makedirs(outbreak_dir, exist_ok=True)
        # pick a random tree from outbreak_trees dir
        tree_path = os.path.join("outbreak_trees", tree_file)
        config_path = create_config_file(
            outbreak_dir, "practice_run.cfg", tree_path, ref.assembly_file
        )
        run_docker(config_path, outbreak_dir)
        print("Outbreak practice dataset created in %s", outbreak_dir)
    # Move tree2read files to another dir
    os.makedirs(os.path.join(outbreak_dir, "tree2readinput"), exist_ok=True)
    for file in os.listdir(outbreak_dir):
        if file.endswith(".nwk") or file.endswith(".fasta") or file.endswith(".cfg"):
            if os.path.exists(os.path.join(outbreak_dir, file)):
                shutil.move(
                    os.path.join(outbreak_dir, file),
                    os.path.join(outbreak_dir, "tree2readinput", file),
                )
    # Move the fastq files into place
    fastq_dir = os.path.join(outbreak_dir, "simulated", "fastq")
    fasta_dir = os.path.join(outbreak_dir, "simulated", "fasta_files")
    # Each sample reads are in a folder named after the sample
    sample_name_to_public = {}
    for sample in os.listdir(fastq_dir):
        if os.path.isdir(os.path.join(fastq_dir, sample)):
            # move fastq to this dir instead
            fasta = os.path.join(fasta_dir, f"{sample}.fasta")
            r1 = os.path.join(fastq_dir, sample, f"{sample}_1.fq.gz")
            r2 = os.path.join(fastq_dir, sample, f"{sample}_2.fq.gz")
            with open(fasta, "rb") as f:
                file_data = f.read()
                md5_hash = hashlib.md5(file_data).hexdigest()
            public_name = f"Sample_{md5_hash[0:8]}12"
            sample_name_to_public[sample] = public_name
            if os.path.exists(r1) and os.path.exists(r2):
                shutil.move(
                    r1, os.path.join(outbreak_dir, f"{public_name}_R1.fastq.gz")
                )
                shutil.move(
                    r2, os.path.join(outbreak_dir, f"{public_name}_R2.fastq.gz")
                )
    # Read metadata and create answer and sample sheet
    meta_file = tree_file.replace(".nwk", "_clusters.csv")
    meta_path = os.path.join("outbreak_trees", meta_file)
    sample_sheet_path = os.path.join(outbreak_dir, "sample_sheet.csv")
    answer_sheet_path = os.path.join(outbreak_dir, "answer_sheet.csv")
    with open(meta_path, mode="r", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        meta_list = list(reader)
    # Create the answer sheet
    for row in meta_list:
        row["public_name"] = sample_name_to_public[row["Sample"]]
        row["R1"] = f"{row['public_name']}_R1.fastq.gz"
        row["R2"] = f"{row['public_name']}_R2.fastq.gz"
    with open(answer_sheet_path, mode="w", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=meta_list[0].keys())
        writer.writeheader()
        writer.writerows(meta_list)
    # Create a blank sample sheet
    with open(sample_sheet_path, mode="w", encoding="utf-8") as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames=[
                "Sample",
                "Cluster",
                "Host",
                "Location",
                "Collection Date",
                "Phenotype",
                "R1",
                "R2",
                "SPECIES",
            ],
        )
        writer.writeheader()
        with open(meta_path, mode="r", encoding="utf-8") as meta_file:
            meta_reader = csv.DictReader(meta_file)
            for row in meta_reader:
                row["Sample"] = sample_name_to_public[row["Sample"]]
                row["Cluster"] = ""
                writer.writerow(row)


def generate_dataset(samplelist, species, output_dir, random_seed=42):
    # Open samplelist and get list of available genomes.
    # filter by species, and split into contaminant and valid.
    # Read the sample list
    sample_list = []
    with open(samplelist, mode="r", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        sample_list = list(reader)

    # Filter by species
    species_dict = [x for x in sample_list if x["species"] == species]

    # ------ Create KLEBORATE dataset ------ #

    kleb_sample_list = create_kleborate_dataset(
        species_dict, output_dir, 5, random_seed
    )
    make_typing_output_files(kleb_sample_list, output_dir, "kleborate_test")
    logging.info(
        "Kleborate test dataset created in %s",
        os.path.join(output_dir, "kleborate_test"),
    )
    real_kleb_sample_list = create_kleborate_dataset(
        species_dict, output_dir, 20, random_seed + random_seed
    )
    make_typing_output_files(real_kleb_sample_list, output_dir, "real_kleborate_test")
    logging.info(
        "Kleborate REAL test dataset created in %s",
        os.path.join(output_dir, "real_kleborate_test"),
    )

    # ------ Create ASSEMBLY dataset ------ #

    random.seed(random_seed)
    contaminant_dict = [x for x in sample_list if x["species"] != species]
    contaminant = AssemblySample(random.choice(contaminant_dict), output_dir)
    contaminant.fetch_assembly()
    contaminant.fetch_reads()
    PRACTICE_ERROR_LIST = ["CONTAMINATED"] + (len(kleb_sample_list) - 1) * ["NORMAL"]
    random.shuffle(PRACTICE_ERROR_LIST)
    assembly_practice_dir = os.path.join(output_dir, "assembly_practice")
    os.makedirs(assembly_practice_dir, exist_ok=True)
    assembly_dataset(
        kleb_sample_list,
        PRACTICE_ERROR_LIST,
        contaminant,
        assembly_practice_dir,
        random_seed,
    )

    # Simulated reads of 20 genomes with some (5-8) having some kind of error.
    random.seed(random_seed + random_seed)
    real_contaminant = AssemblySample(random.choice(contaminant_dict), output_dir)
    real_contaminant.fetch_assembly()
    real_contaminant.fetch_reads()
    assembly_real_dir = os.path.join(output_dir, "assembly_real")
    REAL_ERROR_LIST = [
        "CONTAMINATED",
        "CONTAMINATED",
        "POOR_QUALITY",
        "TRUNCATED",
        "LOW_COVERAGE",
    ] + (len(real_kleb_sample_list) - 5) * ["NORMAL"]
    random.shuffle(REAL_ERROR_LIST)
    assembly_dataset(
        real_kleb_sample_list,
        REAL_ERROR_LIST,
        real_contaminant,
        assembly_real_dir,
        random_seed,
    )

    # ------ Create OUTBREAK dataset ------ #

    outbreak_practice_dir = os.path.join(output_dir, "outbreak_practice")
    random.seed(random_seed)
    ref = random.choice(kleb_sample_list)
    tree_files = [
        f
        for f in os.listdir("outbreak_trees")
        if f.endswith(".nwk")
        and f.replace(".nwk", "_clusters.csv") in os.listdir("outbreak_trees")
    ]
    tree_file = random.choice(tree_files)
    create_outbreak_analysis_dataset(ref, tree_file, outbreak_practice_dir)
    logging.info("Outbreak practice dataset created in %s", outbreak_practice_dir)
    outbreak_real_dir = os.path.join(output_dir, "outbreak_real")
    random.seed(random_seed + random_seed)
    ref = random.choice(real_kleb_sample_list)
    tree_file = random.choice(tree_files)
    create_outbreak_analysis_dataset(ref, tree_file, outbreak_real_dir)
    logging.info("Outbreak REAL dataset created in %s", outbreak_real_dir)


def choose_assembly_type(error_type, sample, output_dir):
    SampleClass = None
    if error_type == "NORMAL":
        SampleClass = AssemblySample
    elif error_type == "CONTAMINATED":
        SampleClass = ContaminatedSample
    elif error_type == "LOW_COVERAGE":
        SampleClass = LowCoverageSample
    elif error_type == "POOR_QUALITY":
        SampleClass = DegradedSample
    elif error_type == "TRUNCATED":
        SampleClass = TruncatedSample
    elif error_type == "CORRUPTED":
        SampleClass = CorruptedSample
    else:
        raise ValueError(f"Unknown error type {error_type}")
    assembly_sample = SampleClass(sample.to_dict(), output_dir)
    with open(assembly_sample.assembly_file, "rb") as f:
        file_data = f.read()
        md5_hash = hashlib.md5(file_data).hexdigest()
        assembly_sample.public_name = (
            f"Sample_{md5_hash[0:8]}{SampleClass.__name__[2:6]}"
        )
    return assembly_sample


def assembly_dataset(samplelist, error_list, contaminant, output_dir, random_seed=42):
    # pick a random int between 40 and 80 for contamination level
    random.seed(random_seed)
    contaminant.fetch_reads()
    sample_sheet_path = os.path.join(output_dir, "sample_sheet.csv")
    answer_sheet_path = os.path.join(output_dir, "answer_sheet.csv")
    os.makedirs(output_dir, exist_ok=True)
    samples_dict_list = []
    sample_short_dict_list = []
    for sample, error_type in zip(samplelist, error_list):
        assembly_sample = choose_assembly_type(error_type, sample, output_dir)
        assembly_sample.fetch_reads()
        if hasattr(assembly_sample, "contamination_sample_r1"):
            assembly_sample.contamination_sample_r1 = contaminant.r1
            assembly_sample.contamination_sample_r2 = contaminant.r2
        assembly_sample.create_modified_reads(
            output_prefix=os.path.join(output_dir, assembly_sample.public_name),
            random_seed=random_seed,
        )
        samples_dict_list.append(assembly_sample.to_dict())
        sample_short_dict_list.append(assembly_sample.to_short_dict())
    # Write the answer sheet
    with open(answer_sheet_path, mode="w", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=samples_dict_list[0].keys())
        writer.writeheader()
        writer.writerows(samples_dict_list)
    # Write the sample sheet
    with open(sample_sheet_path, mode="w", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=sample_short_dict_list[0].keys())
        writer.writeheader()
        writer.writerows(sample_short_dict_list)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate testing dataset for ghru2.")
    # Subparser for simulate reads
    parser.add_argument(
        "--samplelist",
        type=str,
        help="List of samples to use",
        default="samplelist.csv",
    )
    parser.add_argument(
        "--species",
        type=str,
        help="Species to use for generating samples",
        default="K. pneumoniae",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the generated samples",
        default="ghru_output_dataset",
    )
    parser.add_argument(
        "--random_seed", type=int, help="Random seed for reproducibility", default=42
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    generate_dataset(args.samplelist, args.species, args.output_dir, args.random_seed)
