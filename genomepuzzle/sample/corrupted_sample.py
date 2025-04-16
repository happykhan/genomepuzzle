"""
This module contains the LowCoverageSample class which extends BasicSample to 
create modified reads with corrupted bytes.
"""

import logging
import os
import shutil
import random
from genomepuzzle.sample.basic_sample import BasicSample

class CorruptedSample(BasicSample):
    """
    CorruptedSample class for creating and managing corrupted genome reads.

    This class inherits from BasicSample and provides functionality to create modified
    reads by corrupting a random byte in each read file.

    Attributes:
        modified_r1 (str): Path to the modified R1 read file.
        modified_r2 (str): Path to the modified R2 read file.

    Methods:
        __init__(input_dict: dict, output_dir: str, random_seed: int = 42, bin_dir: str = "bin"):
            Initializes the CorruptedSample instance with the given parameters.
        
        create_modified_reads(output_prefix: str, random_seed: int = 42) -> tuple:
        
        _corrupt_random_byte(read_file: str):
            Corrupts a random byte in the specified read file.
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

    def to_dict(self):
        original_dict = super().to_dict()
        original_dict['qc'] = 'FAILED'
        original_dict['error'] = 'CORRUPTED'
        original_dict['notes'] = 'Corrupted sample with random bytes'
        return original_dict

    def to_short_dict(self, public_name=False):
        original_dict = super().to_short_dict(public_name)
        original_dict["r1"] = os.path.basename(self.modified_r1)
        original_dict["r2"] = os.path.basename(self.modified_r2)
        original_dict['sample_name'] = self.public_name
        original_dict['fasta'] = f'{self.public_name}.fasta'
        return original_dict

    def create_modified_reads(self, output_prefix, random_seed=42):
        """
        Creates modified reads by copying the original reads and corrupting a random byte in each.

        Args:
            output_prefix (str): The prefix for the output file names.
            random_seed (int, optional): The seed for the random number generator. Defaults to 42.

        Returns:
            tuple: A tuple containing the paths to the modified R1 and R2 files.

        """
        random.seed(random_seed)
        self.modified_r1 = f"{output_prefix}_R1.fastq.gz"
        self.modified_r2 = f"{output_prefix}_R2.fastq.gz"
        shutil.copy(self.r1, self.modified_r1)
        shutil.copy(self.r2, self.modified_r2)
        logging.info("Corrupting a random byte in %s...", output_prefix)
        self._corrupt_random_byte(self.modified_r1)
        self._corrupt_random_byte(self.modified_r2)
        logging.info("Corruption complete.")
        return self.modified_r1, self.modified_r2

    def _corrupt_random_byte(self, read_file):
        with open(read_file, "r+b") as f:
            f.seek(random.randint(0, os.path.getsize(read_file) - 1))
            f.write(b"\x00")
