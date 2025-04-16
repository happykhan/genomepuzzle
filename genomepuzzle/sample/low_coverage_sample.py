from genomepuzzle.sample.basic_sample import BasicSample
import logging
import random
import os

class LowCoverageSample(BasicSample):
    def __init__(
        self,
        input_dict: dict,
        output_dir: str,
        random_seed: int = 42,
        bin_dir: str = "bin",
    ):
        super().__init__(input_dict, output_dir, random_seed=random_seed, bin_dir=bin_dir)
        self.modified_r1 = None
        self.modified_r2 = None
        self.final_coverage = None

    def to_dict(self):
        original_dict = super().to_dict()
        original_dict['qc'] = 'FAILED'
        original_dict['error'] = 'LOW_COVERAGE'
        original_dict['notes'] = f'Low coverage sample with {self.final_coverage}%'
        return original_dict

    def to_short_dict(self, public_name=False):
        original_dict = super().to_short_dict(public_name)
        original_dict["r1"] = os.path.basename(self.modified_r1)
        original_dict["r2"] = os.path.basename(self.modified_r2)
        original_dict['sample_name'] = self.public_name
        original_dict['fasta'] = f'{self.public_name}.fasta'
        return original_dict


    def create_modified_reads(self, output_prefix, coverage=10, random_seed=None):
        if random_seed:
            random.seed(random_seed)
            coverage = random.randint(1, 20)
        self.final_coverage = coverage
        self.modified_r1 = f"{output_prefix}_R1.fastq.gz"
        self.modified_r2 = f"{output_prefix}_R2.fastq.gz"
        # subsample the reads
        num_reads = int(round(self.count_reads(self.r1)))
        logging.info("Number of reads in input_r1: %d", num_reads)
        logging.info("Subsampling %d reads from real reads (%d percent)", num_reads, coverage)
        self._subsample_paired_read_by_count(
            self.r1, self.r2, self.modified_r1, self.modified_r2, num_reads * (coverage / 100), self.random_seed
        )
        read_subsample_count = self.count_reads(self.modified_r1)
        logging.info("Number of reads in subsample_output_r1: %d", read_subsample_count)
        self.read_count = read_subsample_count
        return self.modified_r1, self.modified_r2
       