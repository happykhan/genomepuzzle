import os
import pytest
import subprocess
from genomepuzzle.simulate_reads import fetch_assembly

def test_fetch_assembly_download(mocker):
    output_dir = 'test_output'
    accessions = ['GCF_000001405.39', 'GCF_000001635.27']

    # Mock os.listdir to simulate the absence of 'ncbi_dataset.zip'
    mocker.patch('os.listdir', return_value=[])

    # Mock subprocess.run to avoid actually running the commands
    mocker.patch('subprocess.run')

    # Run the function
    fetch_assembly(accessions, output_dir)

    # Check that the correct commands were run
    expected_download_command = f"./bin/datasets download genome accession {' '.join(accessions)}"
    expected_move_command = f"mv ncbi_dataset.zip {output_dir}"
    expected_unzip_command = f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}"

    subprocess.run.assert_any_call(expected_download_command, shell=True, check=True)
    subprocess.run.assert_any_call(expected_move_command, shell=True, check=True)
    subprocess.run.assert_any_call(expected_unzip_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def test_fetch_assembly_already_downloaded(mocker):
    output_dir = 'test_output'
    accessions = ['GCF_000001405.39', 'GCF_000001635.27']

    # Mock os.listdir to simulate the presence of 'ncbi_dataset.zip'
    mocker.patch('os.listdir', return_value=['ncbi_dataset.zip'])

    # Mock subprocess.run to avoid actually running the commands
    mocker.patch('subprocess.run')

    # Run the function
    fetch_assembly(accessions, output_dir)

    # Check that the download and move commands were not run
    expected_download_command = f"./bin/datasets download genome accession {' '.join(accessions)}"
    expected_move_command = f"mv ncbi_dataset.zip {output_dir}"
    expected_unzip_command = f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}"

    subprocess.run.assert_any_call(expected_unzip_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run.assert_not_called()

def test_fetch_assembly_unzip_failure(mocker):
    output_dir = 'test_output'
    accessions = ['GCF_000001405.39', 'GCF_000001635.27']

    # Mock os.listdir to simulate the absence of 'ncbi_dataset.zip'
    mocker.patch('os.listdir', return_value=[])

    # Mock subprocess.run to raise an exception for the unzip command
    mocker.patch('subprocess.run', side_effect=[None, None, subprocess.CalledProcessError(1, 'cmd')])

    # Run the function and check that it logs an error
    with pytest.raises(subprocess.CalledProcessError):
        fetch_assembly(accessions, output_dir)