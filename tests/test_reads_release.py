import gzip
import json
import threading
import time
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
    for mate in (1, 2):
        with gzip.open(
            source / "source-a_R{0}.fastq.gz".format(mate), "wt"
        ) as handle:
            handle.write("@source-a/{0}\nACGT\n+\nIIII\n".format(mate))
    if exercise == "hybrid":
        with gzip.open(source / "source-a_long.fastq.gz", "wt") as handle:
            handle.write("@source-a-long\nACGT\n+\nIIII\n")
    spec_path = tmp_path / "release.toml"
    _spec(spec_path, exercise)
    output = tmp_path / "output"

    manifests = package_read_release(
        load_release_spec(spec_path),
        source,
        {"source-a": {"qc": "fail", "species": "K. pneumoniae"}},
        output,
        implant_validations={
            "source-a": {"status": "passed", "checks": ["fixture_validation"]}
        },
    )

    public_text = Path(manifests["public_manifest"]).read_text()
    private = json.loads(Path(manifests["private_manifest"]).read_text())
    assert "source-a" not in public_text
    assert private["samples"][0]["source_id"] == "source-a"
    assert private["samples"][0]["expected_answers"]["qc_status"] == "FAIL"
    assert (
        private["samples"][0]["expected_answers"]["failure_reason"]
        == "LOW_COVERAGE"
    )
    assert (output / "public/sample_sheet.csv").is_file()
    assert (output / "COMPLETE.json").is_file()
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


def test_package_read_release_uses_slurm_cpus_for_independent_samples(
    tmp_path, monkeypatch
):
    source = tmp_path / "source"
    source.mkdir()
    sample_blocks = []
    answers = {}
    for index in range(4):
        source_id = "source-{0}".format(index)
        sample_blocks.append(
            """
[[samples]]
source_id = "{source_id}"
public_id = "Sample_{index}"
""".format(source_id=source_id, index=index)
        )
        answers[source_id] = {"qc": "pass", "species": "K. pneumoniae"}
        for mate in (1, 2):
            with gzip.open(
                source / "{0}_R{1}.fastq.gz".format(source_id, mate), "wt"
            ) as handle:
                handle.write("@{0}/{1}\nACGT\n+\nIIII\n".format(source_id, mate))
    spec_path = tmp_path / "release.toml"
    spec_path.write_text(
        'release_id = "parallel-package"\n'
        'exercise = "assembly"\n'
        'mode = "practice"\n'
        + "\n".join(sample_blocks),
        encoding="utf-8",
    )

    from genomepuzzle import reads_release

    original = reads_release.anonymize_paired_fastq
    thread_ids = set()

    def observed_anonymize(*args, **kwargs):
        thread_ids.add(threading.get_ident())
        time.sleep(0.05)
        return original(*args, **kwargs)

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setattr(reads_release, "anonymize_paired_fastq", observed_anonymize)

    package_read_release(
        load_release_spec(spec_path),
        source,
        answers,
        tmp_path / "output",
    )

    assert len(thread_ids) > 1


@pytest.mark.parametrize(
    ("fault_type", "present_roles", "zero_role"),
    [
        ("MISSING_R2", {"read_1"}, None),
        ("ZERO_BYTE_R1", {"read_1", "read_2"}, "read_1"),
    ],
)
def test_package_read_release_accepts_only_declared_catastrophic_file_faults(
    tmp_path, fault_type, present_roles, zero_role
):
    source = tmp_path / "source"
    source.mkdir()
    r1 = source / "source-a_R1.fastq.gz"
    r2 = source / "source-a_R2.fastq.gz"
    if fault_type != "MISSING_R1":
        if zero_role == "read_1":
            r1.write_bytes(b"")
        else:
            with gzip.open(r1, "wt") as handle:
                handle.write("@source-a/1\nACGT\n+\nIIII\n")
    if fault_type != "MISSING_R2":
        if zero_role == "read_2":
            r2.write_bytes(b"")
        else:
            with gzip.open(r2, "wt") as handle:
                handle.write("@source-a/2\nACGT\n+\nIIII\n")
    spec_path = tmp_path / "release.toml"
    spec_path.write_text(
        """
release_id = "catastrophic-file-fault"
exercise = "assembly"
mode = "practice"

[[samples]]
source_id = "source-a"
public_id = "Sample_fault"
implant = "{fault_type}"
""".format(fault_type=fault_type).strip()
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "output"
    package_read_release(
        load_release_spec(spec_path),
        source,
        {"source-a": {"species": "Klebsiella pneumoniae"}},
        output,
        implant_validations={
            "source-a": {
                "status": "passed",
                "checks": ["catastrophic_file_fault_materialized"],
            }
        },
    )
    manifest = json.loads((output / "public/manifest.json").read_text())
    files = manifest["samples"][0]["files"]
    assert set(files) == present_roles
    if zero_role:
        assert files[zero_role]["size"] == 0
    fault = json.loads(
        (output / "private/implant_manifest.json").read_text()
    )["samples"][0]
    assert fault["fault_type"] == fault_type
    assert fault["failure_reason"] in {"EMPTY_FILE", "MISSING_MATE"}


def test_hybrid_release_can_deliberately_omit_long_reads(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for mate in (1, 2):
        with gzip.open(source / "source-a_R{0}.fastq.gz".format(mate), "wt") as handle:
            handle.write("@source-a/{0}\nACGT\n+\nIIII\n".format(mate))
    spec_path = tmp_path / "release.toml"
    spec_path.write_text(
        """
release_id = "missing-long-reads"
exercise = "hybrid"
mode = "practice"

[[samples]]
source_id = "source-a"
public_id = "Sample_fault"
implant = "MISSING_LONG_READS"
""".strip()
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "output"
    package_read_release(
        load_release_spec(spec_path),
        source,
        {"source-a": {"species": "Klebsiella pneumoniae"}},
        output,
        implant_validations={
            "source-a": {
                "status": "passed",
                "checks": ["missing_long_reads_materialized"],
            }
        },
    )
    manifest = json.loads((output / "public/manifest.json").read_text())
    assert set(manifest["samples"][0]["files"]) == {"read_1", "read_2"}
