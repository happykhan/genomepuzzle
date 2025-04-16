import logging
import subprocess
import os
import tempfile
import re
import random
import gzip
import time
import hashlib

class BasicSample:
    def __init__(
        self,
        input_dict: dict,
        output_dir: str,
        random_seed: int = 42,
        bin_dir: str = "bin",
        work_dir: str = "work",
    ):
        self._check_sample(input_dict)
        os.makedirs(work_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        self.sample_name = input_dict.get("sample_name")
        self.short_reads = input_dict.get("short_reads")
        self.public_name = input_dict.get("public_name", self.sample_name)
        self.assembly = input_dict.get("assembly")
        self.species = input_dict.get("species")
        self.st = input_dict.get("st")
        self.bin_dir = bin_dir
        self.random_seed = random_seed
        self.use_original_reads = input_dict.get("use_original_reads")
        self.output_dir = output_dir
        self.assembly_file = input_dict.get("assembly_file", None)
        self.r1 = input_dict.get("r1", None)
        self.r2 = input_dict.get("r2", None)
        self.coverage = input_dict.get("coverage", random.randint(40, 60))
        self.read_length = 150
        self.platform = "HS25"
        self.fragment_length = 200
        self.standard_deviation = 10
        self.random_seed = 42
        self.qc = "PASSED"
        self.error = "NONE"
        self.notes = "None"
        self.read_count = input_dict.get("read_count", None)
        self.work_dir = work_dir

    def to_dict(self):
        return {
            "sample_name": self.sample_name,
            "short_reads": self.short_reads,
            "public_name": self.public_name,
            "assembly": self.assembly,
            "species": self.species,
            "st": self.st,
            "bin_dir": self.bin_dir,
            "random_seed": self.random_seed,
            "use_original_reads": self.use_original_reads,
            "output_dir": self.output_dir,
            "assembly_file": self.assembly_file,
            "r1": self.r1,
            "r2": self.r2,
            "coverage": self.coverage,
            "read_length": self.read_length,
            "platform": self.platform,
            "fragment_length": self.fragment_length,
            "standard_deviation": self.standard_deviation,
            "qc": self.qc,
            "error": self.error,
            "notes": self.notes,
            "read_count": self.read_count,
        }

    def to_short_dict(self, public_name=False):
        if public_name:
            sample_name = self.public_name
            fasta = f"{self.public_name}.fasta"
        else:
            sample_name = self.sample_name   
            fasta = os.path.basename(self.assembly_file)
        return {
            "sample_name": sample_name,
            "species": self.species,
            "r1": self.r1,
            "r2": self.r2,
            "qc": "",
            "error": "",
            "st": "",
            "notes": "",
            "fasta": fasta,
        }

    def _check_sample(self, input_dict):
        # Check input_dict is a dictionary and has the right keys
        if not isinstance(input_dict, dict):
            raise ValueError("Input must be a dictionary")
        if not all(
            key in input_dict
            for key in [
                "sample_name",
                "short_reads",
                "assembly",
                "species",
                "st",
                "use_original_reads",
            ]
        ):
            raise ValueError(
                "Input dictionary must have the keys: sample_name, short_reads, assembly, species, st, use_original_reads"
            )
        # assembly should look like GCA_902155925.1
        if not re.match(r"GC\w_\d+(\.\d+)?", input_dict["assembly"]):
            raise ValueError("Assembly must match the pattern GC\\w_\\d+.\\d")
        # species should be a string
        if not isinstance(input_dict["species"], str):
            raise ValueError("Species must be a string")
        # ST should be ST followed by a number
        if not re.match(r"ST\d+", input_dict["st"]):
            raise ValueError("ST must match the pattern ST\\d+")
        # use_original_reads should be a FALSE or TRUE
        if input_dict["use_original_reads"] not in ["TRUE", "FALSE"]:
            raise ValueError("use_original_reads must be a boolean")


    def move_assembly_anon(self, output_dir):
        """
        Move the assembly file to the output directory with an anonymous name. Rename the internal contig names to be anonymous.
        """
        if not self.assembly_file:
            self.fetch_assembly()
        # Create a md5 hash of the assembly file
        # Calculate the md5 hash of the assembly file
        with open(self.assembly_file, "rb") as f:
            file_data = f.read()
            md5_hash = hashlib.md5(file_data).hexdigest()
        self.public_name = f"Sample_{md5_hash[0:8]}"
        output_path = os.path.join(output_dir, self.public_name + ".fasta")
        # Rename the contigs to be anonymous
        with open(self.assembly_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        with open(output_path, "w", encoding="utf-8") as f:
            count = 1
            for line in lines:
                if line.startswith(">"):
                    contig_name = f"{self.public_name}_contig_{str(count).zfill(4)}"
                    f.write(f">{contig_name}\n")
                    count += 1
                else:
                    f.write(line)
        return self.public_name


    def fetch_assembly(
        self,
    ):
        assembly_file_path = os.path.join(self.work_dir, self.sample_name + ".fasta")
        if not os.path.exists(assembly_file_path):
            with tempfile.TemporaryDirectory() as output_dir:
                time.sleep(3)
                command = f"./bin/datasets download --api-key b82dffb754d98d8b7bec3b0639fb3bdd1c09 genome accession {self.assembly}"
                subprocess.run(command, shell=True, check=True)
                subprocess.run(f"mv ncbi_dataset.zip {output_dir}", shell=True, check=True)
                logging.info("Downloaded assembly files to %s", output_dir)
                
                # Unzip the downloaded file
                try:
                    subprocess.run(
                        f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}",
                        shell=True,
                        check=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                    logging.info("Unzipped assembly files in %s", output_dir)
                    # Get the assembly file
                    assembly_dir = [os.path.join(output_dir + "/ncbi_dataset/data", f) for f in os.listdir(output_dir + "/ncbi_dataset/data") if f.startswith(self.assembly) ][0]
                    
                    assembly_file = [
                        os.path.join(assembly_dir, f)
                        for f in os.listdir(assembly_dir)
                        if f.endswith(".fna") or f.endswith(".fasta")
                    ][0]
                    # Move assembly file to self.output_dir
                    subprocess.run(
                        f"mv {assembly_file} {assembly_file_path}", shell=True, check=True
                    )
                    self.assembly_file = assembly_file_path
                except subprocess.CalledProcessError as e:
                    logging.error("Failed to unzip assembly files: %s", e)
        else:
            self.assembly_file = assembly_file_path
            logging.info("Assembly file already exists. Skipping assembly download.")

    def fetch_reads(self):
        # Fetch the reads
        if self.use_original_reads == "TRUE":
            # Use the original reads
            self._download_reads()
        else:
            # Simulate reads
            self._run_art()
        # Get read count
        self.read_count = self.count_reads(self.r1)

    def _run_art(self, force=False, num_threads=8):
        output_dir = self.work_dir
        art_r1 = os.path.join(output_dir, f"{self.sample_name}_R1.fq")
        art_r2 = os.path.join(output_dir, f"{self.sample_name}_R2.fq")
        output_r1 = os.path.join(output_dir, f"{self.sample_name}_R1.fastq.gz")
        output_r2 = os.path.join(output_dir, f"{self.sample_name}_R2.fastq.gz")
        self.r1 = output_r1
        self.r2 = output_r2
        if not os.path.exists(output_r1) or not os.path.exists(output_r2) or force:
            # Check assembly file path exists
            if not self.assembly_file:
                raise ValueError("Assembly file path not set")
            if not os.path.exists(self.assembly_file) or not os.path.isfile(self.assembly_file):
                raise ValueError(f"Assembly file {self.assembly_file} does not exist or is not a file")
            command = (
                f"{self.bin_dir}/art_illumina -ss {self.platform} -i {self.assembly_file} "
                f"-l {self.read_length} -f {self.coverage} "
                f"-o {os.path.join(output_dir, f'{self.sample_name}_R')} "
                f"-p -m {self.fragment_length} -s {self.standard_deviation} "
                f"--rndSeed {self.random_seed} -na"
            )
            # rename the output files to the desired names

            subprocess.run(command, shell=True, check=True)
            logging.info("Running command: %s", command)
            # gzip output_r1
            logging.info("gzipping %s...", art_r1)
            subprocess.run(f"pigz -f -p {num_threads} {art_r1}", shell=True, check=True)
            # gzip output_r2
            logging.info("gzipping %s...", art_r2)
            subprocess.run(f"pigz -f -p {num_threads} {art_r2}", shell=True, check=True)
            # rename the files
            os.rename(art_r1 + ".gz", output_r1)
            os.rename(art_r2 + ".gz", output_r2)
        else:
            logging.info("Reads already exist. Skipping read simulation.")

    def _download_reads(self):
        output_dir = self.output_dir
        fastqdump_output_r1 = os.path.join(output_dir, f"{self.sample_name}_1.fastq")
        fastqdump_output_r2 = os.path.join(output_dir, f"{self.sample_name}_2.fastq")
        output_r1 = os.path.join(output_dir, f"{self.sample_name}_R1.fastq.gz")
        output_r2 = os.path.join(output_dir, f"{self.sample_name}_R2.fastq.gz")
        logging.info("Fetching reads for %s", self.sample_name)
        command = f"{self.bin_dir}/fasterq-dump --outdir {output_dir} --split-files {self.short_reads}"
        subprocess.run(command, shell=True, check=True)
        logging.info("Downloaded reads for %s", self.sample_name)

        logging.info("gzipping %s...", fastqdump_output_r1)
        subprocess.run(f"gzip -f {fastqdump_output_r1}", shell=True, check=True)
        logging.info("gzipping %s...", fastqdump_output_r2)
        subprocess.run(f"gzip -f {fastqdump_output_r2}", shell=True, check=True)
        # rename the files
        os.rename(fastqdump_output_r1 + ".gz", output_r1)
        os.rename(fastqdump_output_r2 + ".gz", output_r2)
        self.r1 = output_r1
        self.r2 = output_r2

    @staticmethod
    def count_reads(fastq_file):
        """
        Count the number of reads in a FASTQ file.

        Parameters:
        - fastq_file: str, path to the FASTQ file.

        Returns:
        - int, number of reads in the FASTQ file.
        """
        count = 0
        # Check file exists and is a file
        if not os.path.exists(fastq_file) or not os.path.isfile(fastq_file):
            raise ValueError(f"File {fastq_file} does not exist or is not a file")
        with gzip.open(fastq_file, "rt") as f:
            for _ in f:
                count += 1
        return count // 4


    @staticmethod
    def _subsample_paired_read_by_count(input_r1, input_r2, output_r1, output_r2, num_reads=1000, random_seed=42):
        """
        Subsamples paired-end reads by a specified count using seqtk.

        This function takes paired-end read files (input_r1 and input_r2), subsamples
        them to a specified number of reads (num_reads) using the seqtk tool, and
        writes the subsampled reads to output files (output_r1 and output_r2). The
        subsampling process is deterministic if a random seed (random_seed) is provided.

        Args:
            input_r1 (str): Path to the input file for read 1.
            input_r2 (str): Path to the input file for read 2.
            output_r1 (str): Path to the output file for subsampled read 1.
            output_r2 (str): Path to the output file for subsampled read 2.
            num_reads (int, optional): Number of reads to subsample. Defaults to 1000.
            random_seed (int, optional): Seed for the random number generator to ensure
                                         reproducibility. Defaults to 42.

        Raises:
            subprocess.CalledProcessError: If the seqtk command fails.

        Logs:
            Info: Logs the progress and completion of the subsampling process.
        """
 
        random.seed(random_seed)
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

        logging.info("Subsampling by count complete. Output files saved as:")
        logging.info("  %s", output_r1)
        logging.info("  %s", output_r2)
