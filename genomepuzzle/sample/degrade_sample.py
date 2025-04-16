import random
import gzip
import logging
from genomepuzzle.sample.basic_sample import BasicSample
import os
import subprocess


class DegradedSample(BasicSample):
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
        self.decrement = None

    def to_dict(self):
        original_dict = super().to_dict()
        original_dict["qc"] = "FAILED"
        original_dict["error"] = "DEGRADED"
        original_dict["notes"] = (
            f"Degraded sample so base quality is decremented {self.decrement} points"
        )
        return original_dict

    def to_short_dict(self, public_name=False):
        original_dict = super().to_short_dict(public_name)
        original_dict["r1"] = os.path.basename(self.modified_r1)
        original_dict["r2"] = os.path.basename(self.modified_r2)
        original_dict["sample_name"] = self.public_name
        original_dict['fasta'] = f'{self.public_name}.fasta'
        return original_dict

    def create_modified_reads(self, output_prefix, decrement=20, random_seed=None):
        """
        Create modified reads with degraded quality scores.

        This method generates two new FASTQ files with modified quality scores
        within the specified range. The quality scores are degraded randomly
        using a specified seed for reproducibility.

        Args:
            output_prefix (str): The prefix for the output FASTQ files.
            min_quality (int): The minimum quality score for degradation (inclusive).
            max_quality (int): The maximum quality score for degradation (inclusive).
            random_seed (int): The seed for the random number generator to ensure reproducibility.

        Raises:
            ValueError: If the quality scores are not in the range 0-93.

        Returns:
            None
        """
        if random_seed is None:
            random.seed(random_seed)
            decrement = random.randint(10, 30)
        self.decrement = decrement
        self.modified_r1 = f"{output_prefix}_R1.fastq.gz"
        self.modified_r2 = f"{output_prefix}_R2.fastq.gz"

        logging.info("Degrading read quality in %s by %d ...", output_prefix, decrement)
        self._degrade_fastq_file(self.r1, self.modified_r1, decrement)
        self._degrade_fastq_file(self.r2, self.modified_r2, decrement)
        logging.info("Quality degradation complete. Output saved to: %s", output_prefix)
        return self.modified_r1, self.modified_r2

    def _degrade_fastq_file(self, fastq_file, output_fastq, decrement, random_seed=42):
        command = f"bin/seqtk seq -Q{decrement} {fastq_file} | gzip > {output_fastq}",
        subprocess.run(
            command,
            shell=True,
            check=True,
        )
