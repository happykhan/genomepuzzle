import os
import pytest
import subprocess
from genomepuzzle.simulate_reads import run_art

def test_run_art(mocker):
    sample = {
        'platform': 'HS25',
        'read_length': 150,
        'coverage': 50,
        'public_name': 'sample1',
        'fragment_length': 200,
        'standard_deviation': 10,
        'random_seed': 42
    }
    output_dir = 'test_output'
    reference_genome = 'test_reference.fna'
    output_r1 = os.path.join(output_dir, 'sample1_R1.fq')
    output_r2 = os.path.join(output_dir, 'sample1_R2.fq')

    # Mock subprocess.run to avoid actually running the commands
    mocker.patch('subprocess.run')

    # Run the function
    run_art(sample, output_dir, reference_genome, output_r1, output_r2)

    # Check that the correct commands were run
    expected_command = f"bin/art_illumina -ss HS25 -i test_reference.fna -l 150 -f 50 -o {os.path.join(output_dir, 'sample1_R')} -p -m 200 -s 10 --rndSeed 42 -na"
    subprocess.run.assert_any_call(expected_command, shell=True, check=True)
    subprocess.run.assert_any_call(f"gzip -f {output_r1}", shell=True, check=True)
    subprocess.run.assert_any_call(f"gzip -f {output_r2}", shell=True, check=True)

def test_run_art_invalid_command(mocker):
    sample = {
        'platform': 'HS25',
        'read_length': 150,
        'coverage': 50,
        'public_name': 'sample1',
        'fragment_length': 200,
        'standard_deviation': 10,
        'random_seed': 42
    }
    output_dir = 'test_output'
    reference_genome = 'test_reference.fna'
    output_r1 = os.path.join(output_dir, 'sample1_R1.fq')
    output_r2 = os.path.join(output_dir, 'sample1_R2.fq')

    # Mock subprocess.run to raise an exception
    mocker.patch('subprocess.run', side_effect=subprocess.CalledProcessError(1, 'cmd'))

    # Run the function and check that it raises an exception
    with pytest.raises(subprocess.CalledProcessError):
        run_art(sample, output_dir, reference_genome, output_r1, output_r2)