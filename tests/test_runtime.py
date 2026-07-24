import subprocess
from unittest.mock import Mock, call, patch

from genomepuzzle.runtime import seqtk_sample


def test_seqtk_sample_connects_files_directly_to_avoid_pipe_deadlock(tmp_path):
    input_path = tmp_path / "input.fastq"
    output_path = tmp_path / "output.fastq.gz"
    input_path.write_text("@read\nACGT\n+\nIIII\n")
    sampler = Mock()
    sampler.wait.return_value = 0
    compressor = Mock()
    compressor.wait.return_value = 0

    with patch(
        "genomepuzzle.runtime.require_tool",
        side_effect=["/pixi/bin/seqtk", "/pixi/bin/pigz"],
    ):
        with patch(
            "genomepuzzle.runtime.subprocess.Popen",
            side_effect=[sampler, compressor],
        ) as popen:
            seqtk_sample(input_path, output_path, 42, 0.1)

    assert popen.call_args_list[0] == call(
        [
            "/pixi/bin/seqtk",
            "sample",
            "-s",
            "42",
            str(input_path),
            "0.1",
        ],
        stdout=subprocess.PIPE,
    )
    assert popen.call_args_list[1].kwargs["stdin"] is sampler.stdout
    assert popen.call_args_list[1].args[0] == [
        "/pixi/bin/pigz",
        "-n",
        "-p",
        "1",
        "-c",
    ]
    sampler.stdout.close.assert_called_once()
