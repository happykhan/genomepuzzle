"""
This module contains the AssemblySample class which extends 
BasicSample but is just the normal BasicSample class.
"""

import os
import shutil
from genomepuzzle.sample.basic_sample import BasicSample

class AssemblySample(BasicSample):
    
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
        original_dict['qc'] = 'PASSED'
        original_dict['error'] = 'No error'
        original_dict['notes'] = 'Its fine.'
        return super().to_dict()
    
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
        self.modified_r1 = f"{output_prefix}_R1.fastq.gz"
        self.modified_r2 = f"{output_prefix}_R2.fastq.gz"
        shutil.copy(self.r1, self.modified_r1)
        shutil.copy(self.r2, self.modified_r2)
        return self.modified_r1, self.modified_r2
