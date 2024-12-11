import pytest
import sys
from io import StringIO
from genomepuzzle.main import parse_arguments

def test_parse_arguments_simulate():
    test_args = [
        'main.py', 'simulate', '--num_samples', '5', '--samplelist', 'test_samplelist.csv',
        '--species', 'E. coli', '--output_dir', 'test_output', '--random_seed', '123'
    ]
    sys.argv = test_args
    args = parse_arguments()
    assert args.command == 'simulate'
    assert args.num_samples == 5
    assert args.samplelist == 'test_samplelist.csv'
    assert args.species == 'E. coli'
    assert args.output_dir == 'test_output'
    assert args.random_seed == 123

def test_parse_arguments_errors():
    test_args = [
        'main.py', 'errors', '--samplelist', 'test_output/sample_sheet.csv', '--error_proportion', '0.7',
        '--random_seed', '123', '--contamination_list', 'test_contaminants.csv', '--output_dir', 'test_final_output'
    ]
    sys.argv = test_args
    args = parse_arguments()
    assert args.command == 'errors'
    assert args.samplelist == 'test_output/sample_sheet.csv'
    assert args.error_proportion == 0.7
    assert args.random_seed == 123
    assert args.contamination_list == 'test_contaminants.csv'
    assert args.output_dir == 'test_final_output'

def test_parse_arguments_no_command():
    test_args = ['main.py']
    sys.argv = test_args
    with pytest.raises(SystemExit):
        parse_arguments()