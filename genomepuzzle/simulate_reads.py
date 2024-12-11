import logging
import os
import csv
import subprocess
import random
from genomepuzzle.util import check_input_table
from genomepuzzle.util import check_operating_system


def cleanup_output_dir(output_dir):
    """
    Cleans up the specified output directory by removing specific files and directories.

    Args:
        output_dir (str): The path to the output directory to be cleaned.

    Removes:
        - ncbi_dataset.zip: A zip file in the output directory.
        - ncbi_dataset: A directory in the output directory.
        - md5sum.txt: A file in the output directory, if it exists.
        - README.md: A file in the output directory, if it exists.

    Logs:
        - Information about the cleanup process, including which files and directories were removed.
    """
    logging.info("Cleaning up output directory %s", output_dir)
    os.remove(os.path.join(output_dir, "ncbi_dataset.zip"))
    ncbi_dataset_dir = os.path.join(output_dir, "ncbi_dataset")
    if os.path.exists(ncbi_dataset_dir):
        subprocess.run(f"rm -rf {ncbi_dataset_dir}", shell=True, check=True)
        logging.info("Removed directory %s", ncbi_dataset_dir)
    # remove md5sum.txt if it exists
    md5sum_file = os.path.join(output_dir, "md5sum.txt")
    if os.path.exists(md5sum_file):
        os.remove(md5sum_file)
        logging.info("Removed file %s", md5sum_file)
    # remove README.md if it exists
    readme_file = os.path.join(output_dir, "README.md")
    if os.path.exists(readme_file):
        os.remove(readme_file)
        logging.info("Removed file %s", readme_file)


def fetch_assembly(accessions, output_dir):
    """
    Fetches genome assembly files for given accessions and extracts them to the specified output directory.

    Args:
        accessions (list of str): List of genome accession numbers to download.
        output_dir (str): Directory where the downloaded and extracted files will be stored.

    Raises:
        subprocess.CalledProcessError: If there is an error during the download or extraction process.

    Logs:
        Info: When the assembly files are downloaded, already present, or successfully unzipped.
        Error: If there is a failure during the unzipping process.
    """
    if "ncbi_dataset.zip" not in os.listdir(output_dir):
        accessions = " ".join(accessions)
        command = f"./bin/datasets download genome accession {accessions}"
        subprocess.run(command, shell=True, check=True)
        # Move the downloaded file to the output directory
        subprocess.run(f"mv ncbi_dataset.zip {output_dir}", shell=True, check=True)
        logging.info("Downloaded assembly files to %s", output_dir)
    else:
        logging.info("Assembly files already downloaded in %s", output_dir)
    # Unzip the downloaded file
    try:
        subprocess.run(
            f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}",
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as e:
        logging.error("Failed to unzip assembly files: %s", e)
    logging.info("Unzipped assembly files in %s", output_dir)


def run_art(sample, output_dir, reference_genome, output_r1, output_r2):
    """
    Simulates sequencing reads using the ART Illumina tool and processes the output files.

    Args:
        sample (dict): A dictionary containing sample information with keys:
            - 'public_name' (str): The public name of the sample.
            - 'platform' (str): The sequencing platform.
            - 'read_length' (int): The length of the reads.
            - 'coverage' (float): The desired coverage.
            - 'fragment_length' (int): The fragment length.
            - 'standard_deviation' (float): The standard deviation of the fragment length.
            - 'random_seed' (int): The random seed for reproducibility.
        output_dir (str): The directory where the output files will be stored.
        reference_genome (str): The path to the reference genome file.
        output_r1 (str): The desired path for the gzipped R1 output file.
        output_r2 (str): The desired path for the gzipped R2 output file.

    Raises:
        subprocess.CalledProcessError: If the ART Illumina command or gzip command fails.
    """
    art_r1 = os.path.join(output_dir, f"{sample['public_name']}_R1.fq")
    art_r2 = os.path.join(output_dir, f"{sample['public_name']}_R2.fq")
    command = (
        f"bin/art_illumina -ss {sample['platform']} -i {reference_genome} "
        f"-l {sample['read_length']} -f {sample['coverage']} "
        f"-o {os.path.join(output_dir, f'{sample['public_name']}_R')} "
        f"-p -m {sample['fragment_length']} -s {sample['standard_deviation']} "
        f"--rndSeed {sample['random_seed']} -na"
    )
    # rename the output files to the desired names

    subprocess.run(command, shell=True, check=True)
    logging.info("Running command: %s", command)
    # gzip output_r1
    logging.info("gzipping %s...", art_r1)
    subprocess.run(f"gzip -f {art_r1}", shell=True, check=True)
    # gzip output_r2
    logging.info("gzipping %s...", art_r2)
    subprocess.run(f"gzip -f {art_r2}", shell=True, check=True)
    # rename the files
    os.rename(art_r1 + ".gz", output_r1)
    os.rename(art_r2 + ".gz", output_r2)


def download_reads(sample, output_dir, output_r1, output_r2):
    """
    Downloads and processes sequencing reads for a given sample.

    This function fetches sequencing reads using the `fasterq-dump` tool, compresses
    the downloaded files using gzip, and renames the compressed files to specified output names.

    Args:
        sample (dict): A dictionary containing sample information. Must include the keys:
            - 'SHORT_READS': The identifier for the short reads to be downloaded.
            - 'SAMPLE_NAME': The name of the sample.
        output_dir (str): The directory where the downloaded and processed files will be stored.
        output_r1 (str): The desired output file path for the first read file (R1).
        output_r2 (str): The desired output file path for the second read file (R2).

    Raises:
        subprocess.CalledProcessError: If any of the subprocess commands fail.
        OSError: If renaming the files fails.
    """
    fastqdump_output_r1 = os.path.join(output_dir, f"{sample['SHORT_READS']}_1.fastq")
    fastqdump_output_r2 = os.path.join(output_dir, f"{sample['SHORT_READS']}_2.fastq")
    logging.info("Fetching reads for %s", sample["SAMPLE_NAME"])
    command = (
        f"bin/fasterq-dump --outdir {output_dir} --split-files {sample['SHORT_READS']}"
    )
    subprocess.run(command, shell=True, check=True)
    logging.info("Downloaded reads for %s", sample["SAMPLE_NAME"])

    logging.info("gzipping %s...", fastqdump_output_r1)
    subprocess.run(f"gzip -f {fastqdump_output_r1}", shell=True, check=True)
    logging.info("gzipping %s...", fastqdump_output_r2)
    subprocess.run(f"gzip -f {fastqdump_output_r2}", shell=True, check=True)
    # rename the files
    os.rename(fastqdump_output_r1 + ".gz", output_r1)
    os.rename(fastqdump_output_r2 + ".gz", output_r2)


def fetch_reads(all_sample_list, output_dir):
    """
    Fetches reads for each sample in the provided list and stores them
    in the specified output directory.

    For each sample in `all_sample_list`, if the sample is marked to use
    original reads and has short reads, this function checks if the reads
    are already downloaded. If not, it downloads the reads and updates the
    sample dictionary with the paths to the reads and additional metadata.

    Args:
        all_sample_list (list): A list of dictionaries, each representing
        a sample with various attributes.
        output_dir (str): The directory where the reads should be stored.

    Returns:
        list: The updated list of sample dictionaries with paths to the
        reads and additional metadata.

    Raises:
        ValueError: If a sample is marked to use original reads but does
        not have a read record.
    """
    for sample in all_sample_list:
        if sample.get("USE_ORIGINAL_READS").upper() == "TRUE":
            if sample.get("SHORT_READS"):
                output_r1 = os.path.join(
                    output_dir, f"{sample['SAMPLE_NAME']}_R1.fastq.gz"
                )
                output_r2 = os.path.join(
                    output_dir, f"{sample['SAMPLE_NAME']}_R2.fastq.gz"
                )
                if not os.path.exists(output_r1) or not os.path.exists(output_r2):
                    download_reads(sample, output_dir, output_r1, output_r2)
                else:
                    logging.info(
                        "Reads already downloaded for %s", sample["SAMPLE_NAME"]
                    )
                sample["r1"] = output_r1
                sample["r2"] = output_r2
                sample["coverage"] = -1
                sample["read_length"] = -1
                sample["platform"] = "Unknown"
                sample["fragment_length"] = -1
                sample["standard_deviation"] = -1
                sample["random_seed"] = -1
                sample["QC"] = "PASSED"
                sample["ERROR"] = "NONE"
                sample["Notes"] = "None"
            else:
                logging.error(
                    "Sample %s does not have a read record.", sample["SAMPLE_NAME"]
                )
                raise ValueError(
                    f"Sample {sample['SAMPLE_NAME']} does not have a read record."
                )
    return all_sample_list


def simulate_reads(num_samples, sample_list, species, output_dir, random_seed=42):
    """
    Simulate sequencing reads for a given species.

    This function generates simulated sequencing reads for a specified number of samples
    from a given species. It creates an output directory if it doesn't exist, fetches
    reads or assemblies as needed, and generates the reads using ART. The function also
    writes a sample sheet and cleans up the output directory.

    Args:
        num_samples (int): The number of samples to generate.
        sample_list (str): Path to the input sample list file.
        species (str): The species for which to generate samples.
        output_dir (str): The directory where output files will be saved.
        random_seed (int, optional): The seed for random number generation. Defaults to 42.

    Returns:
        None
    """
    check_operating_system()
    # import sample list as dict
    all_sample_list = check_input_table(sample_list)
    random.seed(random_seed)
    # create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Created output directory at %s", output_dir)

    # Create a reduced sample list
    species_sample_list = [x for x in all_sample_list if x.get("SPECIES") == species]
    print(f"Generating {num_samples} samples for {species} in {output_dir}")
    if num_samples < len(species_sample_list):
        species_sample_list = random.sample(species_sample_list, num_samples)
    elif num_samples > len(species_sample_list):
        species_sample_list.extend(
            random.choice(species_sample_list)
            for _ in range(num_samples - len(species_sample_list))
        )

    species_sample_list = fetch_reads(species_sample_list, output_dir)
    # download the assembly genomes for use from the sample list
    # (where sequence reads are not available)
    assembly_list = [x["ASSEMBLY"] for x in species_sample_list if not x.get("r1")]
    if assembly_list:
        fetch_assembly(assembly_list, output_dir)
        # generate samples
        for sample in species_sample_list:
            sample["public_name"] = sample["SAMPLE_NAME"]
            reference_genome = [
                os.path.join(output_dir, "ncbi_dataset/data/", x)
                for x in os.listdir(os.path.join(output_dir, "ncbi_dataset/data/"))
                if x.startswith(sample["ASSEMBLY"])
            ][0]
            reference_genome = [
                os.path.join(reference_genome, x)
                for x in os.listdir(reference_genome)
                if x.endswith(".fna")
            ][0]
            output_r1 = os.path.join(output_dir, f"{sample['public_name']}_R1.fastq.gz")
            output_r2 = os.path.join(output_dir, f"{sample['public_name']}_R2.fastq.gz")

            sample["r1"] = output_r1
            sample["r2"] = output_r2
            # coverage is a random number between 40 and 60
            sample["coverage"] = random.randint(40, 60)
            sample["read_length"] = 150
            sample["platform"] = "HS25"
            sample["fragment_length"] = 200
            sample["standard_deviation"] = 10
            sample["random_seed"] = 42
            sample["QC"] = "PASSED"
            sample["ERROR"] = "NONE"
            sample["Notes"] = "None"
            if not os.path.exists(output_r1) or not os.path.exists(output_r2):
                # use ART
                run_art(sample, output_dir, reference_genome, output_r1, output_r2)
            else:
                logging.info("Reads already generated for %s", sample["SAMPLE_NAME"])
    # write sample sheet
    sample_sheet = os.path.join(output_dir, "sample_sheet.csv")
    with open(sample_sheet, "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=species_sample_list[0].keys())
        writer.writeheader()
        for sample in species_sample_list:
            sample["r1"] = os.path.basename(sample["r1"])
            sample["r2"] = os.path.basename(sample["r2"])
            writer.writerow(sample)
    # clean up output directory
    cleanup_output_dir(output_dir)
    logging.info(
        "Finished generating %d samples for %s in %s", num_samples, species, output_dir
    )
