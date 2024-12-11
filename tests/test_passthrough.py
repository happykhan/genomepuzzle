import pytest
from genomepuzzle.create_error import pass_through

def test_pass_through(tmp_path):
    # Create temporary input files
    r1_content = "@SEQ_ID\nGATTTGGGGTTTCCCAGTCACGACGTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n"
    r2_content = "@SEQ_ID\nGATTTGGGGTTTCCCAGTCACGACGTT\n+\n!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n"
    r1_path = tmp_path / "input_r1.fastq"
    r2_path = tmp_path / "input_r2.fastq"
    r1_output = tmp_path / "output_r1.fastq"
    r2_output = tmp_path / "output_r2.fastq"

    with open(r1_path, "w") as f:
        f.write(r1_content)
    with open(r2_path, "w") as f:
        f.write(r2_content)

    # Run the function
    pass_through(str(r1_path), str(r2_path), str(r1_output), str(r2_output))

    # Check that the output files exist and have the correct content
    assert r1_output.exists()
    assert r2_output.exists()

    with open(r1_output, "r") as f:
        assert f.read() == r1_content

    with open(r2_output, "r") as f:
        assert f.read() == r2_content

def test_pass_through_nonexistent_file(tmp_path):
    r1_path = tmp_path / "nonexistent_r1.fastq"
    r2_path = tmp_path / "nonexistent_r2.fastq"
    r1_output = tmp_path / "output_r1.fastq"
    r2_output = tmp_path / "output_r2.fastq"

    with pytest.raises(FileNotFoundError):
        pass_through(str(r1_path), str(r2_path), str(r1_output), str(r2_output))