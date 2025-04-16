import random
import os
import csv
import logging
from genomepuzzle.sample import Sample
from genomepuzzle.sample_util import run_spades, run_badread, run_flye
from genomepuzzle.assembly_stats import calculate_assembly_stats


def add_contamination():
    pass

def make_ref_sample(sample, output_dir):
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
    return sample, output_r1, output_r2, reference_genome


def contamination_menu(
    num_samples: int,
    samplelist: str,
    species: str,
    contamination_type: str,
    output_dir: str,
    assemble: str,
    random_seed: int,
) -> None:
    """Generate contamination samples"""
    random.seed(random_seed)
    with open(samplelist, mode='r', encoding='utf-8') as infile:
        reader = csv.DictReader(infile)
        samples = [row for row in reader]
    main_sample_dict = random.choice([s for s in samples if s["species"] == species])
    if contamination_type == "Species":
        contamination_sample_dict = random.choice([s for s in samples if s["species"] != species])
    elif contamination_type == "ST":
        contamination_sample_dict = random.choice([s for s in samples if s["species"] == species and s["st"] != main_sample_dict["st"]])
    else:
        raise ValueError("Invalid type")
    main_sample = Sample(main_sample_dict, output_dir)
    contamination_sample = Sample(contamination_sample_dict, output_dir)
    contamination_sample.coverage = 20
    main_sample.coverage = 20
    main_sample.fetch_assembly()
    main_sample.fetch_reads()
    contamination_sample.fetch_assembly()
    contamination_sample.fetch_reads()
    # Generate num_samples samples with the contamination.
    proportion_increment = 100 / num_samples
    output_table = [] 
    answer_table = [] 
    for i in range(num_samples):
        # Proportion increment should be 100 / num_samples
        percentage = proportion_increment * (i + 1)
        public_name = f"sample{str(i+1).zfill(2)}"
        output_r1 = os.path.join(output_dir, f"{public_name}_R1.fastq.gz")
        output_r2 = os.path.join(output_dir, f"{public_name}_R2.fastq.gz")
        main_sample.add_contamination(
            contamination_sample.r1,
            contamination_sample.r2,
            percentage,
            output_r1,
            output_r2
        )
        sheet = main_sample.to_short_dict()
        sheet['sample_name'] = public_name
        sheet['r1'] = output_r1
        sheet['r2'] = output_r2
        answer = main_sample.to_dict()
        answer['public_name'] = public_name
        answer['r1'] = output_r1
        answer['r2'] = output_r2
        if percentage >= 10:
            answer['qc'] = "FAILED"
            answer['error'] = "CONTAMINATION"
        answer['notes'] = f"Contamination with {contamination_sample.species} - {contamination_sample.sample_name} at {percentage}%"
        output_table.append(sheet)
        answer_table.append(answer)
        # Save the samples in the output_dir
        if assemble:
            # Need to create fresh genome assemblies
            spades_results = run_spades(output_r1, output_r2, output_dir, public_name)
            answer.update(spades_results)

    output_table_path = os.path.join(output_dir, "sample_sheet.csv")
    answer_table_path = os.path.join(output_dir, "answer_sheet.csv")
    with open(output_table_path, mode='w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_table[0].keys())
        writer.writeheader()
        writer.writerows(output_table)

    with open(answer_table_path, mode='w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=answer_table[0].keys())
        writer.writeheader()
        writer.writerows(answer_table)
    # Remove original reads from the dir
    os.remove(main_sample.r1)
    os.remove(main_sample.r2)
    os.remove(contamination_sample.r1)
    os.remove(contamination_sample.r2)