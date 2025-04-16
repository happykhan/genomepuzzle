"""
This module contains the ContaminatedSample class which extends BasicSample to 
create reads contaminated with some other sample.
"""
import logging
import os
import subprocess
import random
from genomepuzzle.sample.basic_sample import BasicSample

class ContaminatedSample(BasicSample):
    """
    A class used to represent a contaminated sample, inheriting from BasicSample.

    Attributes
    ----------
    modified_r1 : str
        Path to the modified R1 read file.
    modified_r2 : str
        Path to the modified R2 read file.
    contamination_percent : float
        Percentage of contamination in the sample.
    contamination_sample_r1 : str
        Path to the contamination sample R1 read file.
    contamination_sample_r2 : str
        Path to the contamination sample R2 read file.

    Methods
    -------
    to_dict():
        Returns a dictionary representation of the sample with contamination details.
    to_short_dict(public_name=False):
        Returns a short dictionary representation of the sample with contamination details.
    create_modified_reads(output_prefix, percentage=60, random_seed=None):
        Creates modified read files by mixing original reads with contaminant reads.
    """

    def __init__(
        self,
        input_dict: dict,
        output_dir: str,
        random_seed: int = 42,
        bin_dir: str = "bin",
    ):
        super().__init__(
            input_dict, output_dir, random_seed=random_seed, bin_dir=bin_dir
        )
        self.modified_r1 = None
        self.modified_r2 = None
        self.contamination_percent = None
        self.contamination_sample_r1 = None
        self.contamination_sample_r2 = None

    def to_dict(self):
        original_dict = super().to_dict()
        original_dict["qc"] = "FAILED"
        original_dict["error"] = "CONTAMINATED"
        original_dict["notes"] = (
            f"Contamined sample with {self.contamination_percent}% using {os.path.basename(self.contamination_sample_r1)} and {os.path.basename(self.contamination_sample_r2)}"
        )
        return original_dict

    def to_short_dict(self, public_name=False):
        original_dict = super().to_short_dict(public_name)
        original_dict["r1"] = os.path.basename(self.modified_r1)
        original_dict["r2"] = os.path.basename(self.modified_r2)
        original_dict["sample_name"] = self.public_name
        original_dict['fasta'] = f'{self.public_name}.fasta'
        return original_dict

    def create_modified_reads(self, output_prefix, percentage=60, random_seed=None):
        # if random seed, override with a random percent
        if random_seed:
            random.seed(random_seed)
            percentage = random.randint(60, 90)
        self.contamination_percent = percentage
        self.modified_r1 = f"{output_prefix}_R1.fastq.gz"
        self.modified_r2 = f"{output_prefix}_R2.fastq.gz"
        if not os.path.exists(self.modified_r1) and not os.path.exists(
            self.modified_r2
        ):
            # mix the reads with the contaminant simulated reads
            # how many reads are in input_r1 ?
            num_reads = int(round(self.count_reads(self.r1)))
            contaminated_num_reads = (percentage / 100) * num_reads
            original_num_reads = num_reads - contaminated_num_reads
            logging.info("Number of reads in input_r1: %d", num_reads)
            logging.info(
                "Number of good reads to downsample to: %d", original_num_reads
            )
            logging.info("Number of bad reads to add in: %d", contaminated_num_reads)
            # pick subsample adj_num_reads from contaminant_output_r1 and contaminant_output_r2
            logging.info(
                "Subsampling %d reads from real reads (%d percent)",
                original_num_reads,
                100 - percentage,
            )
            subsample_output_r1 = os.path.join(
                self.output_dir, "subsample_real_output_r1.fastq.gz"
            )
            subsample_output_r2 = os.path.join(
                self.output_dir, "subsample_real_output_r2.fastq.gz"
            )
            self._subsample_paired_read_by_count(
                self.r1,
                self.r2,
                subsample_output_r1,
                subsample_output_r2,
                original_num_reads,
                self.random_seed,
            )
            read_subsample_count = self.count_reads(subsample_output_r1)
            logging.info(
                "Number of reads in subsample_output_r1: %d", read_subsample_count
            )
            logging.info(
                "Subsampling %d reads from contaminant reads (%d percent)",
                contaminated_num_reads,
                percentage,
            )
            subsample_contaminant_output_r1 = os.path.join(
                self.output_dir, "subsample_contaminant_output_r1.fastq.gz"
            )
            subsample_contaminant_output_r2 = os.path.join(
                self.output_dir, "subsample_contaminant_output_r2.fastq.gz"
            )
            self._subsample_paired_read_by_count(
                self.contamination_sample_r1,
                self.contamination_sample_r2,
                subsample_contaminant_output_r1,
                subsample_contaminant_output_r2,
                contaminated_num_reads,
                self.random_seed,
            )
            read_subsample_count = self.count_reads(subsample_contaminant_output_r1)
            logging.info(
                "Number of reads in subsample_contaminant_output_r1: %d",
                read_subsample_count,
            )

            # concatenate the subsampled reads
            logging.info("Appending the contaminant reads to the output files - r1")
            subprocess.run(
                f"cat {subsample_contaminant_output_r1} {subsample_output_r1} >> {self.modified_r1}",
                shell=True,
                check=True,
            )
            logging.info("Appending the contaminant reads to the output files - r2")
            subprocess.run(
                f"cat {subsample_contaminant_output_r2} {subsample_output_r2} >> {self.modified_r2}",
                shell=True,
                check=True,
            )
            read_subsample_count = self.count_reads(self.modified_r1)
            logging.info("Number of reads in final output: %d", read_subsample_count)
            os.remove(subsample_output_r1)
            os.remove(subsample_output_r2)
            os.remove(subsample_contaminant_output_r1)
            os.remove(subsample_contaminant_output_r2)
        else:
            logging.info("Modified reads already exist, skipping creation")
        return self.modified_r1, self.modified_r2
