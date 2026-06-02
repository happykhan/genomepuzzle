import csv
import logging
import os
import random

from genomepuzzle.create_error import contamination as contaminate_reads
from genomepuzzle.sample.assembly_sample import AssemblySample
from genomepuzzle.sample_util import run_spades


def _load_samples(samplelist):
    with open(samplelist, mode="r", encoding="utf-8") as infile:
        return [row for row in csv.DictReader(infile)]


def _choose_samples(samples, species, contamination_type, rng):
    main_sample_dict = rng.choice([sample for sample in samples if sample["species"] == species])
    if contamination_type == "Species":
        contaminant_pool = [sample for sample in samples if sample["species"] != species]
    elif contamination_type == "ST":
        contaminant_pool = [
            sample
            for sample in samples
            if sample["species"] == species and sample["st"] != main_sample_dict["st"]
        ]
    else:
        raise ValueError("Invalid contamination type")
    if not contaminant_pool:
        raise ValueError("No contamination sample available for the requested settings")
    return main_sample_dict, rng.choice(contaminant_pool)


def contamination_menu(
    num_samples,
    samplelist,
    species,
    contamination_type,
    output_dir,
    assemble,
    random_seed,
):
    """Generate contamination samples."""
    rng = random.Random(random_seed)
    os.makedirs(output_dir, exist_ok=True)
    samples = _load_samples(samplelist)
    main_sample_dict, contaminant_sample_dict = _choose_samples(
        samples, species, contamination_type, rng
    )

    main_sample = AssemblySample(main_sample_dict, output_dir, random_seed=random_seed)
    contaminant_sample = AssemblySample(
        contaminant_sample_dict, output_dir, random_seed=random_seed + 1
    )
    main_sample.fetch_assembly()
    main_sample.fetch_reads()
    contaminant_sample.fetch_assembly()

    proportion_increment = 100.0 / num_samples
    output_table = []
    answer_table = []
    for index in range(num_samples):
        percentage = proportion_increment * (index + 1)
        public_name = "sample{count}".format(count=str(index + 1).zfill(2))
        output_r1 = os.path.join(output_dir, "{name}_R1.fastq.gz".format(name=public_name))
        output_r2 = os.path.join(output_dir, "{name}_R2.fastq.gz".format(name=public_name))
        contaminate_reads(
            main_sample.r1,
            main_sample.r2,
            output_r1,
            output_r2,
            contaminant_sample.assembly,
            output_dir,
            percentage=percentage,
            random_seed=random_seed + index,
        )

        answer = main_sample.to_dict()
        answer["public_name"] = public_name
        answer["r1"] = os.path.basename(output_r1)
        answer["r2"] = os.path.basename(output_r2)
        answer["qc"] = "FAILED" if percentage >= 10 else "PASSED"
        answer["error"] = "CONTAMINATION"
        answer["notes"] = (
            "Contamination with {species} - {sample_name} at {percentage:.1f}%".format(
                species=contaminant_sample.species,
                sample_name=contaminant_sample.sample_name,
                percentage=percentage,
            )
        )
        if assemble:
            answer.update(run_spades(output_r1, output_r2, output_dir, public_name))
        answer_table.append(answer)

        sample_row = {
            "sample_name": public_name,
            "species": main_sample.species,
            "r1": os.path.basename(output_r1),
            "r2": os.path.basename(output_r2),
            "qc": "",
            "error": "",
            "st": "",
            "notes": "",
            "fasta": "{name}.fasta".format(name=public_name),
        }
        output_table.append(sample_row)

    output_table_path = os.path.join(output_dir, "sample_sheet.csv")
    answer_table_path = os.path.join(output_dir, "answer_sheet.csv")
    with open(output_table_path, mode="w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_table[0].keys())
        writer.writeheader()
        writer.writerows(output_table)
    with open(answer_table_path, mode="w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=answer_table[0].keys())
        writer.writeheader()
        writer.writerows(answer_table)

    logging.info(
        "Generated %d contamination samples in %s using %s as contaminant",
        num_samples,
        output_dir,
        contaminant_sample.assembly,
    )
