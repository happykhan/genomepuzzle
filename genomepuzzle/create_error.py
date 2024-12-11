
import random
import csv 
import logging
import os
import shutil
import gzip 

def degrade_quality(input_fastq, output_fastq, min_quality=5, max_quality=30, random_seed=42):
    """
    Degrade base quality scores randomly, with poorer quality towards the end of each read.

    Parameters:
    - input_fastq: str, path to the input FASTQ file.
    - output_fastq: str, path to the output FASTQ file.
    - min_quality: int, minimum Phred score for base quality at the end of the read.
    - max_quality: int, maximum Phred score for base quality at the start of the read.
    """
    random.seed(random_seed)
    if not (0 <= min_quality <= max_quality <= 93):  # Ensure valid Phred+33 range
        raise ValueError("Quality scores must be in range 0-93.")

    with gzip.open(input_fastq, "rt") as infile, gzip.open(output_fastq, "wt") as outfile:
        while True:
            # Read one FASTQ record (4 lines)
            header = infile.readline().strip()  # Line 1: Sequence ID
            sequence = infile.readline().strip()  # Line 2: Nucleotide sequence
            plus = infile.readline().strip()  # Line 3: +
            quality = infile.readline().strip()  # Line 4: Quality scores

            # Break if end of file
            if not header:
                break

            read_length = len(sequence)
            degraded_quality = []

            # Generate progressively poorer quality scores
            for i in range(read_length):
                # Calculate a score range that worsens towards the end of the read
                current_max_quality = max_quality - int((i / read_length) * (max_quality - min_quality))
                random_quality = random.randint(min_quality, current_max_quality)
                degraded_quality.append(chr(random_quality + 33))  # Convert Phred score to ASCII

            # Join the degraded quality scores
            new_quality = "".join(degraded_quality)

            # Write the modified FASTQ record
            outfile.write(f"{header}\n{sequence}\n{plus}\n{new_quality}\n")

    print(f"Quality degradation complete. Output saved to: {output_fastq}")


def truncate_fastq(input_fastq, output_fastq, truncate_length=75):
    """
    Truncate reads in a FASTQ file to a specific length.

    Parameters:
    - sample: dict, sample information containing 'r1' and 'r2' keys for read files.
    - output_fastq: str, path to the output truncated FASTQ file.
    - truncate_length: int, length to truncate each read to.
    """

    with open(input_fastq, "r") as infile, open(output_fastq, "w") as outfile:
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
            outfile.write(f"{header}\n{truncated_sequence}\n{plus}\n{truncated_quality}\n")

    print(f"Truncation complete. Output saved to: {output_fastq}")

def subsample_paired_fastq(input_r1, input_r2, output_r1, output_r2, subsample_fraction=0.1, random_seed=42):
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

    with open(input_r1, "r") as r1, open(input_r2, "r") as r2, \
         open(output_r1, "w") as out_r1, open(output_r2, "w") as out_r2:

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

    print("Subsampling complete. Output files saved as:")
    print(f"  {output_r1}")
    print(f"  {output_r2}")

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
    # get list of contaminants
    contaminant_list = []
    for x in csv.DictReader(open(contamination_list_file, 'r')):
        if x['SPECIES'] not in species_list and x.get('ASSEMBLY'):
            contaminant_list.append(x)

def corrupt_read(sample, output_dir, random_seed=42):
    random.seed(random_seed)
    file_to_corrupt = os.path.join(output_dir, os.path.basename(sample['r1']))
    with open(file_to_corrupt, 'r+b') as f:
        f.seek(random.randint(0, os.path.getsize(file_to_corrupt) - 1))
        f.write(b'\x00')
    return sample


def introduce_errors(sample_sheet, error_proportion, contamination_list_file, output_dir, random_seed=42):
    # Reads are created and in sample_sheet.csv
    # Now we need to introduce errors in the reads
    # error_proportion is the proportion of samples that will have errors
    # get sample details from file with dict reader
    # create output dir if doesnt exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    random.seed(random_seed)
    types_of_errors = ['TRUNCATED', 'CORRUPT', 'CONTAMINATED', 'POOR_QUALITY', 'LOW_COVERAGE']

    full_sample_list = [ ]   
    with open(sample_sheet, 'r') as f:
        full_sample_list = list(csv.DictReader(f))
    # pick with samples to alter given the error_proportion
    # get list of species in sample_sheet
    species_list = set([sample['SPECIES'] for sample in full_sample_list])

    error_samples = random.sample(full_sample_list, int(len(full_sample_list) * error_proportion))
    for sample in error_samples:
        if random.random() < 0.5:
            sample['ERROR'] = 'CONTAMINATED'
        elif random.random() < 0.3:
            sample['ERROR'] = 'LOW_COVERAGE'
        else:
            sample['ERROR'] = random.choice(types_of_errors)
        sample['QC'] = 'FAILED'
        sample['Notes'] = 'Error introduced'
    
    for sample in full_sample_list:
        input_dir = os.path.dirname(sample_sheet)
        input_r1 = os.path.join(input_dir, os.path.basename(sample['r1']))
        input_r2 = os.path.join(input_dir, os.path.basename(sample['r2']))
        output_r1 = os.path.join(output_dir, os.path.basename(sample['r1']))
        output_r2 = os.path.join(output_dir, os.path.basename(sample['r2']))
        error_type = sample['ERROR']
        if error_type == 'NONE':
            pass_through(input_r1, input_r2, output_r1, output_r2)
            sample['Notes'] += " No changes."
        elif error_type == 'TRUNCATED':
            truncate_length = random.randint(5, 100)
            file_to_truncate = input_r1
            if random.choice([True, False]):
                file_to_truncate = input_r1
                sample['Notes'] = 'Truncated read 1 to {} bp'.format(truncate_length)
            else:
                sample['Notes'] = 'Truncated read 2 to {} bp'.format(truncate_length)  
                file_to_truncate = input_r2         
            sample = truncate_fastq(file_to_truncate, output_dir, truncate_length)
        elif error_type == 'CORRUPT':
            corruption = random.choice(["r1", "r2", "both"])
            if corruption == "both":
                corrupt_read(sample['r1'])
                corrupt_read(sample['r2'])
                sample['Notes'] = 'Corrupted both reads'
            elif corruption == "r1":
                corrupt_read(sample['r1'])
                sample['Notes'] = 'Corrupted read 1'
            elif corruption == "r2":
                corrupt_read(sample['r2'])   
                sample['Notes'] = 'Corrupted read 2'            
        elif error_type == 'CONTAMINATED':
            sample = pass_through(sample, sample_sheet, output_dir, random_seed)
        elif error_type == 'POOR_QUALITY':
            min_quality = random.randint(5, 9)
            max_quality = random.randint(10, 20)
            degrade_quality(input_r1, output_r1, min_quality=min_quality, max_quality=max_quality, random_seed=random_seed)
            degrade_quality(input_r2, output_r2, min_quality=min_quality, max_quality=max_quality, random_seed=random_seed)
            sample['Notes'] += f" Quality degraded to {min_quality}-{max_quality}"
        elif error_type == 'LOW_COVERAGE':
            subsample_fraction = random.uniform(0.1, 0.4)
            subsample_paired_fastq(input_r1, input_r2, output_r1, output_r2, subsample_fraction, random_seed=random_seed)
        else:
            logging.error(f"Error type {error_type} not implemented")
            raise ValueError(f"Error type {error_type} not implemented")
        # write sample sheet
    sample_sheet = os.path.join(output_dir, 'sample_sheet.csv')
    with open(sample_sheet, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=full_sample_list[0].keys())
        writer.writeheader()
        for sample in full_sample_list:
            writer.writerow(sample)

