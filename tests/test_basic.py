import pytest
import os
import csv
from genomepuzzle.util import check_input_table

def create_sample_list(file_path, fieldnames, rows):
    with open(file_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

def test_check_input_table_valid():
    sample_list_file = 'valid_sample_list.csv'
    fieldnames = ['SAMPLE_NAME', 'SHORT_READS', 'ASSEMBLY', 'SPECIES', 'ST', 'AMR', 'USE_ORIGINAL_READS']
    rows = [
        {'SAMPLE_NAME': 'sample1', 'SHORT_READS': 'reads1', 'ASSEMBLY': 'assembly1', 'SPECIES': 'species1', 'ST': 'st1', 'AMR': 'amr1', 'USE_ORIGINAL_READS': 'TRUE'},
        {'SAMPLE_NAME': 'sample2', 'SHORT_READS': 'reads2', 'ASSEMBLY': 'assembly2', 'SPECIES': 'species2', 'ST': 'st2', 'AMR': 'amr2', 'USE_ORIGINAL_READS': 'TRUE'}
    ]
    create_sample_list(sample_list_file, fieldnames, rows)
    assert len(check_input_table(sample_list_file)) == 2
    os.remove(sample_list_file)

def test_check_input_table_invalid():
    sample_list_file = 'invalid_sample_list.csv'
    fieldnames = ['SAMPLE_NAME', 'SHORT_READS', 'ASSEMBLY', 'SPECIES', 'USE_ORIGINAL_READS']
    rows = [
        {'SAMPLE_NAME': 'sample1', 'SHORT_READS': 'reads1', 'ASSEMBLY': 'assembly1', 'SPECIES': 'species1', 'USE_ORIGINAL_READS': 'TRUE'},
        {'SAMPLE_NAME': 'sample2', 'SHORT_READS': 'reads2', 'ASSEMBLY': 'assembly2', 'SPECIES': 'species2', 'USE_ORIGINAL_READS': 'TRUE'}
    ]
    create_sample_list(sample_list_file, fieldnames, rows)
    with pytest.raises(ValueError, match="Input file is missing required columns: AMR, ST"):
        check_input_table(sample_list_file)
    os.remove(sample_list_file)