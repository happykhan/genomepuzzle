import pytest
import os
from genomepuzzle.create_error import truncate_fastq

def create_fastq_file(file_path, records):
    with open(file_path, 'w') as f:
        for record in records:
            f.write(f"@{record['header']}\n{record['sequence']}\n+\n{record['quality']}\n")

def test_truncate_fastq():
    input_fastq = 'test_input.fastq'
    output_fastq = 'test_output.fastq'
    truncate_length = 5

    records = [
        {'header': 'SEQ_ID_1', 'sequence': 'ACTGACTGACTG', 'quality': 'IIIIIIIIIIII'},
        {'header': 'SEQ_ID_2', 'sequence': 'TGACTGACTGAC', 'quality': 'JJJJJJJJJJJJ'}
    ]
    create_fastq_file(input_fastq, records)

    truncate_fastq(input_fastq, output_fastq, truncate_length)

    with open(output_fastq, 'r') as f:
        lines = f.readlines()

    assert lines[1].strip() == 'ACTGA'
    assert lines[3].strip() == 'IIIII'
    assert lines[5].strip() == 'TGACT'
    assert lines[7].strip() == 'JJJJJ'

    os.remove(input_fastq)
    os.remove(output_fastq)

def test_truncate_fastq_longer_than_sequence():
    input_fastq = 'test_input_long.fastq'
    output_fastq = 'test_output_long.fastq'
    truncate_length = 20

    records = [
        {'header': 'SEQ_ID_1', 'sequence': 'ACTGACTGACTG', 'quality': 'IIIIIIIIIIII'},
        {'header': 'SEQ_ID_2', 'sequence': 'TGACTGACTGAC', 'quality': 'JJJJJJJJJJJJ'}
    ]
    create_fastq_file(input_fastq, records)

    truncate_fastq(input_fastq, output_fastq, truncate_length)

    with open(output_fastq, 'r') as f:
        lines = f.readlines()

    assert lines[1].strip() == 'ACTGACTGACTG'
    assert lines[3].strip() == 'IIIIIIIIIIII'
    assert lines[5].strip() == 'TGACTGACTGAC'
    assert lines[7].strip() == 'JJJJJJJJJJJJ'

    os.remove(input_fastq)
    os.remove(output_fastq)