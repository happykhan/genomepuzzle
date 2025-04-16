import logging
import gzip
import random
import os
from genomepuzzle.sample.basic_sample import BasicSample


class TruncatedSample(BasicSample):
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
        self.truncate_length = None

    def to_dict(self):
        original_dict = super().to_dict()
        original_dict["qc"] = "FAILED"
        original_dict["error"] = "TRUNCATED"
        original_dict["notes"] = f"Truncated reads to {self.truncate_length}%"
        return original_dict

    def to_short_dict(self, public_name=False):
        original_dict = super().to_short_dict(public_name)
        original_dict["r1"] = os.path.basename(self.modified_r1)
        original_dict["r2"] = os.path.basename(self.modified_r2)
        original_dict["sample_name"] = self.public_name
        original_dict['fasta'] = f'{self.public_name}.fasta'
        return original_dict

    def create_modified_reads(
        self, output_prefix, truncate_length=40, random_seed=None
    ):
        self.modified_r1 = f"{output_prefix}_R1.fastq.gz"
        self.modified_r2 = f"{output_prefix}_R2.fastq.gz"
        if random_seed:
            random.seed(random_seed)
            truncate_length = random.randint(1, 50)
        self.truncate_length = truncate_length
        logging.info(
            "Truncating reads in %s to %d nucleotides...",
            self.modified_r1,
            truncate_length,
        )
        for input_fastq, output_fastq in [
            (self.r1, self.modified_r1),
            (self.r2, self.modified_r2),
        ]:
            with (
                gzip.open(input_fastq, "rt") as infile,
                gzip.open(output_fastq, "wt") as outfile,
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
        return self.modified_r1, self.modified_r2
