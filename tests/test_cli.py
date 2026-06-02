import pytest
import sys
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
    args = parse_arguments()
    assert args.command is None


def test_parse_arguments_hybrid():
    test_args = [
        'main.py', 'hybrid', '--samplelist', 'hybrid.csv',
        '--contamination_list', 'contaminants.csv', '--mode', 'practice',
        '--output_dir', 'hybrid_output', '--random_seed', '321'
    ]
    sys.argv = test_args
    args = parse_arguments()
    assert args.command == 'hybrid'
    assert args.samplelist == 'hybrid.csv'
    assert args.contamination_list == 'contaminants.csv'
    assert args.mode == 'practice'
    assert args.output_dir == 'hybrid_output'
    assert args.random_seed == 321


def test_parse_arguments_contamination():
    test_args = [
        'main.py', 'contamination', '--num_samples', '6', '--samplelist', 'samples.csv',
        '--species', 'K. pneumoniae', '--type', 'Species',
        '--output_dir', 'contam_output', '--random_seed', '9'
    ]
    sys.argv = test_args
    args = parse_arguments()
    assert args.command == 'contamination'
    assert args.num_samples == 6
    assert args.samplelist == 'samples.csv'
    assert args.species == 'K. pneumoniae'
    assert args.type == 'Species'
    assert args.output_dir == 'contam_output'
    assert args.random_seed == 9
    assert args.assemble is False
