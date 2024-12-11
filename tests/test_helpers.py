import os
import pytest
import shutil
from genomepuzzle.simulate_reads import cleanup_output_dir, fetch_assembly
import subprocess
from unittest.mock import patch, call


@pytest.fixture
def setup_output_dir(tmp_path):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    (output_dir / 'ncbi_dataset.zip').touch()
    (output_dir / 'ncbi_dataset').mkdir()
    (output_dir / 'md5sum.txt').touch()
    (output_dir / 'README.md').touch()
    yield output_dir

def test_cleanup_output_dir(setup_output_dir):
    output_dir = setup_output_dir
    cleanup_output_dir(output_dir)
    
    assert not (output_dir / 'ncbi_dataset.zip').exists()
    assert not (output_dir / 'ncbi_dataset').exists()
    assert not (output_dir / 'md5sum.txt').exists()
    assert not (output_dir / 'README.md').exists()

def test_cleanup_output_dir_partial_files(setup_output_dir):
    output_dir = setup_output_dir
    (output_dir / 'ncbi_dataset.zip').unlink()
    cleanup_output_dir(output_dir)
    
    assert not (output_dir / 'ncbi_dataset').exists()
    assert not (output_dir / 'md5sum.txt').exists()
    assert not (output_dir / 'README.md').exists()

def test_cleanup_output_dir_no_files(setup_output_dir):
    output_dir = setup_output_dir
    shutil.rmtree(output_dir / 'ncbi_dataset')
    (output_dir / 'ncbi_dataset.zip').unlink()
    (output_dir / 'md5sum.txt').unlink()
    (output_dir / 'README.md').unlink()
    cleanup_output_dir(output_dir)
    
    assert not (output_dir / 'ncbi_dataset.zip').exists()
    assert not (output_dir / 'ncbi_dataset').exists()
    assert not (output_dir / 'md5sum.txt').exists()
    assert not (output_dir / 'README.md').exists()

@pytest.fixture
def setup_output_dir(tmp_path):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir

def test_fetch_assembly_download(setup_output_dir):
    accessions = ["GCF_000001405.39", "GCF_000001635.27"]
    output_dir = setup_output_dir

    with patch("subprocess.run") as mock_run:
        fetch_assembly(accessions, output_dir)
        mock_run.assert_has_calls([
            call("./bin/datasets download genome accession GCF_000001405.39 GCF_000001635.27", shell=True, check=True),
            call(f"mv ncbi_dataset.zip {output_dir}", shell=True, check=True),
            call(f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        ])

def test_fetch_assembly_already_downloaded(setup_output_dir):
    output_dir = setup_output_dir
    (output_dir / "ncbi_dataset.zip").touch()

    with patch("subprocess.run") as mock_run:
        fetch_assembly(["GCF_000001405.39"], output_dir)
        mock_run.assert_called_once_with(f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)