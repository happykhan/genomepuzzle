import json
from pathlib import Path

import pytest

from genomepuzzle.reads_release import package_read_release
from genomepuzzle.release import load_release_spec


def _spec(path: Path, exercise: str) -> None:
    path.write_text(
        """
release_id = "{0}-2026"
exercise = "{0}"
mode = "practice"

[[samples]]
source_id = "source-a"
public_id = "Sample_a"
implant = "LOW_COVERAGE"
""".format(exercise).strip()
        + "\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize("exercise", ["assembly", "hybrid"])
def test_package_read_release(tmp_path, exercise):
    source = tmp_path / "source"
    source.mkdir()
    (source / "source-a_R1.fastq.gz").write_bytes(b"r1")
    (source / "source-a_R2.fastq.gz").write_bytes(b"r2")
    if exercise == "hybrid":
        (source / "source-a_long.fastq.gz").write_bytes(b"long")
    spec_path = tmp_path / "release.toml"
    _spec(spec_path, exercise)
    output = tmp_path / "output"

    manifests = package_read_release(
        load_release_spec(spec_path),
        source,
        {"source-a": {"qc": "fail", "species": "K. pneumoniae"}},
        output,
    )

    public_text = Path(manifests["public_manifest"]).read_text()
    private = json.loads(Path(manifests["private_manifest"]).read_text())
    assert "source-a" not in public_text
    assert private["samples"][0]["source_id"] == "source-a"
    assert private["samples"][0]["expected_answers"]["qc"] == "fail"
    assert (output / "public/sample_sheet.csv").is_file()
    assert (output / "COMPLETE").is_file()
    expected_file_count = 3 if exercise == "hybrid" else 2
    assert len(list((output / "public/files").iterdir())) == expected_file_count


def test_package_read_release_requires_answers(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    spec_path = tmp_path / "release.toml"
    _spec(spec_path, "assembly")
    with pytest.raises(ValueError, match="answers missing"):
        package_read_release(
            load_release_spec(spec_path), source, {}, tmp_path / "output"
        )
