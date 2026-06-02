import os
import pytest
import subprocess
from genomepuzzle.simulate_reads import fetch_assembly
from genomepuzzle.runtime import resolve_tool

def test_fetch_assembly_download(mocker):
    output_dir = 'test_output'
    accessions = ['GCF_000001405.39', 'GCF_000001635.27']

    # Mock os.listdir to simulate the absence of 'ncbi_dataset.zip'
    mocker.patch('os.listdir', return_value=[])
    mocker.patch('os.path.exists', return_value=True)
    mocker.patch('shutil.move')

    # Mock subprocess.run to avoid actually running the commands
    mocker.patch('subprocess.run')

    # Run the function
    fetch_assembly(accessions, output_dir)

    # Check that the correct commands were run
    expected_download_command = [resolve_tool('datasets'), 'download', 'genome', 'accession'] + accessions
    expected_unzip_command = ['unzip', '-o', os.path.join(output_dir, 'ncbi_dataset.zip'), '-d', output_dir]

    subprocess.run.assert_any_call(expected_download_command, check=True)
    subprocess.run.assert_any_call(expected_unzip_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def test_fetch_assembly_unzip_failure(mocker):
    output_dir = 'test_output'
    accessions = ['GCF_000001405.39', 'GCF_000001635.27']

    # Mock os.listdir to simulate the absence of 'ncbi_dataset.zip'
    mocker.patch('os.listdir', return_value=[])
    mocker.patch('os.path.exists', return_value=True)
    mocker.patch('shutil.move')

    # Mock subprocess.run to raise an exception for the unzip command
    mocker.patch('subprocess.run', side_effect=[None, subprocess.CalledProcessError(1, 'cmd')])

    # Run the function and check that it logs an error
    fetch_assembly(accessions, output_dir)
