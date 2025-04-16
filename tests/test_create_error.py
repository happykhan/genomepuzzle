import os
from genomepuzzle.create_error import truncate_fastq, contamination, write_sample_sheet


def create_fastq_file(file_path, records):
    with open(file_path, "w", encoding="utf-8") as f:
        for record in records:
            f.write(
                f"@{record['header']}\n{record['sequence']}\n+\n{record['quality']}\n"
            )


def test_truncate_fastq():
    input_fastq = "test_input.fastq"
    output_fastq = "test_output.fastq"
    truncate_length = 5

    records = [
        {"header": "SEQ_ID_1", "sequence": "ACTGACTGACTG", "quality": "IIIIIIIIIIII"},
        {"header": "SEQ_ID_2", "sequence": "TGACTGACTGAC", "quality": "JJJJJJJJJJJJ"},
    ]
    create_fastq_file(input_fastq, records)

    truncate_fastq(input_fastq, output_fastq, truncate_length)

    with open(output_fastq, "r", encoding="utf-8") as f:
        lines = f.readlines()

    assert lines[1].strip() == "ACTGA"
    assert lines[3].strip() == "IIIII"
    assert lines[5].strip() == "TGACT"
    assert lines[7].strip() == "JJJJJ"

    os.remove(input_fastq)
    os.remove(output_fastq)


def test_truncate_fastq_longer_than_sequence():
    input_fastq = "test_input_long.fastq"
    output_fastq = "test_output_long.fastq"
    truncate_length = 20

    records = [
        {"header": "SEQ_ID_1", "sequence": "ACTGACTGACTG", "quality": "IIIIIIIIIIII"},
        {"header": "SEQ_ID_2", "sequence": "TGACTGACTGAC", "quality": "JJJJJJJJJJJJ"},
    ]
    create_fastq_file(input_fastq, records)

    truncate_fastq(input_fastq, output_fastq, truncate_length)

    with open(output_fastq, "r") as f:
        lines = f.readlines()

    assert lines[1].strip() == "ACTGACTGACTG"
    assert lines[3].strip() == "IIIIIIIIIIII"
    assert lines[5].strip() == "TGACTGACTGAC"
    assert lines[7].strip() == "JJJJJJJJJJJJ"

    os.remove(input_fastq)
    os.remove(output_fastq)


def test_contamination():
    input_r1 = "test_input_r1.fastq"
    input_r2 = "test_input_r2.fastq"
    output_r1 = "test_output_r1.fastq"
    output_r2 = "test_output_r2.fastq"
    contaminant_assembly = "test_contaminant_assembly.fasta"
    output_dir = "."
    percentage = 50

    # Create dummy input FASTQ files
    records_r1 = [
        {"header": "SEQ_ID_1", "sequence": "ACTGACTGACTG", "quality": "IIIIIIIIIIII"},
        {"header": "SEQ_ID_2", "sequence": "TGACTGACTGAC", "quality": "JJJJJJJJJJJJ"},
    ]
    records_r2 = [
        {"header": "SEQ_ID_1", "sequence": "GTCAGTCAGTCA", "quality": "KKKKKKKKKKKK"},
        {"header": "SEQ_ID_2", "sequence": "CAGTCAGTCAGT", "quality": "LLLLLLLLLLLL"},
    ]
    create_fastq_file(input_r1, records_r1)
    create_fastq_file(input_r2, records_r2)

    # Create dummy contaminant assembly file
    with open(contaminant_assembly, "w", encoding="utf-8") as f:
        seq = "ACTATACGAGAGATACATACATACATATA" * 30
        f.write(f">contaminant\n{seq}")

    contamination(
        input_r1,
        input_r2,
        output_r1,
        output_r2,
        contaminant_assembly,
        output_dir,
        percentage,
        random_seed=42,
    )
    # delete the dummy files
    os.remove(input_r1)
    os.remove(input_r2)
    os.remove(contaminant_assembly)
    os.remove(output_r1)
    os.remove(output_r2)


def test_write_sample_sheet():
    output_dir = "."
    full_sample_list = [
        {
            "public_name": "sample01",
            "r1": "sample01_R1.fastq.gz",
            "r2": "sample01_R2.fastq.gz",
            "SPECIES": "species1",
            "QC": "Unknown",
            "ERROR": "Unknown",
            "ST": "Unknown",
            "AMR": "Unknown",
            "Notes": "",
        },
        {
            "public_name": "sample02",
            "r1": "sample02_R1.fastq.gz",
            "r2": "sample02_R2.fastq.gz",
            "SPECIES": "species2",
            "QC": "Unknown",
            "ERROR": "Unknown",
            "ST": "Unknown",
            "AMR": "Unknown",
            "Notes": "",
        },
    ]

    write_sample_sheet(output_dir, full_sample_list)
    expected_output = [
        "ID,R1,R2,SPECIES,QC,ERROR,ST,AMR,Notes\n",
        "sample01,sample01_R1.fastq.gz,sample01_R2.fastq.gz,species1,Unknown,Unknown,Unknown,Unknown,\n",
        "sample02,sample02_R1.fastq.gz,sample02_R2.fastq.gz,species2,Unknown,Unknown,Unknown,Unknown,\n",
    ]

    with open(os.path.join(output_dir, "sample_sheet.csv"), "r", encoding="utf-8") as f:
        lines = f.readlines()

    assert lines == expected_output

    os.remove(os.path.join(output_dir, "sample_sheet.csv"))
