import gzip
import json
import threading
import time

from genomepuzzle.read_generation import generate_read_release
from genomepuzzle.release import load_release_spec


def _fake_short(reference, r1, r2, **_):
    reference_name = reference.read_text().splitlines()[0][1:]
    for path, mate in ((r1, 1), (r2, 2)):
        with gzip.open(path, "wt") as handle:
            for index in range(20):
                handle.write(
                    "@{0}-{1}/{2}\nACGT\n+\nIIII\n".format(
                        reference_name, index, mate
                    )
                )


def test_spec_driven_assembly_generation_needs_no_answer_json(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    (source / "source-a.fasta").write_text(">source\n" + "A" * 500 + "\n")
    spec_path = tmp_path / "assembly.toml"
    spec_path.write_text(
        """
release_id = "assembly-generated"
exercise = "assembly"
mode = "practice"

[[samples]]
source_id = "source-a"
public_id = "Sample_a"
implant = "LOW_COVERAGE"
[samples.implant_parameters]
read_fraction = 0.02
[samples.expected_answers]
species = "Klebsiella pneumoniae"
""".strip()
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "genomepuzzle.read_generation._simulate_short_reads", _fake_short
    )

    output = tmp_path / "release"
    generate_read_release(load_release_spec(spec_path), source, output)
    provenance = json.loads((output / "private/provenance.json").read_text())

    assert (output / "COMPLETE.json").is_file()
    validation = provenance["samples"][0]["provenance"]["validation"]
    assert validation["status"] == "passed"
    assert validation["read_fraction"] == 0.02
    assert validation["expected_short_coverage"] == 0.6
    with gzip.open(output / "public/files/Sample_a_R1.fastq.gz", "rt") as handle:
        header = handle.readline()
    assert "Sample_a" in header
    assert "source-a" not in header


def test_too_few_reads_uses_explicit_private_pair_count(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    (source / "source-a.fasta").write_text(">source\n" + "A" * 500 + "\n")
    spec_path = tmp_path / "assembly.toml"
    spec_path.write_text(
        """
release_id = "assembly-truncated"
exercise = "assembly"
mode = "practice"

[[samples]]
source_id = "source-a"
public_id = "Sample_a"
implant = "TRUNCATE_TO_READ_PAIRS"
[samples.implant_parameters]
retained_read_pairs = 7
[samples.expected_answers]
species = "Klebsiella pneumoniae"
""".strip()
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "genomepuzzle.read_generation._simulate_short_reads", _fake_short
    )

    output = tmp_path / "release"
    generate_read_release(load_release_spec(spec_path), source, output)
    provenance = json.loads((output / "private/provenance.json").read_text())
    fault = provenance["samples"][0]["fault"]
    validation = provenance["samples"][0]["provenance"]["validation"]

    assert fault == {
        "failure_reason": "TOO_FEW_READS",
        "fault_type": "TRUNCATE_TO_READ_PAIRS",
        "parameters": {"retained_read_pairs": 7},
    }
    assert validation["short_read_pairs"] == 7
    for mate in ("R1", "R2"):
        with gzip.open(
            output / "public/files/Sample_a_{0}.fastq.gz".format(mate), "rt"
        ) as handle:
            assert sum(1 for _line in handle) == 28


def test_slurm_allocation_generates_independent_samples_concurrently(
    tmp_path, monkeypatch
):
    spec_path = tmp_path / "assembly.toml"
    samples = "\n".join(
        """
[[samples]]
source_id = "source-{index}"
public_id = "Sample_{index}"
[samples.expected_answers]
species = "Klebsiella pneumoniae"
""".format(index=index)
        for index in range(4)
    )
    spec_path.write_text(
        'release_id = "parallel"\nexercise = "assembly"\nmode = "practice"\n'
        + samples,
        encoding="utf-8",
    )
    thread_ids = set()

    def fake_generate(spec, sample, source_path, work_dir):
        thread_ids.add(threading.get_ident())
        time.sleep(0.05)
        return sample.source_id, {"species": "Klebsiella pneumoniae"}, {
            "status": "passed"
        }

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setattr(
        "genomepuzzle.read_generation._generate_sample", fake_generate
    )
    monkeypatch.setattr(
        "genomepuzzle.read_generation.package_read_release",
        lambda *args, **kwargs: {"status": "ok"},
    )

    result = generate_read_release(
        load_release_spec(spec_path), tmp_path / "sources", tmp_path / "release"
    )

    assert result == {"status": "ok"}
    assert len(thread_ids) > 1


def test_contamination_is_mixed_by_exact_final_read_counts(tmp_path, monkeypatch):
    source = tmp_path / "sources"
    source.mkdir()
    (source / "target.fasta").write_text(">target\n" + "A" * 500 + "\n")
    (source / "ecoli.fasta").write_text(">ecoli\n" + "C" * 400 + "\n")
    spec_path = tmp_path / "assembly.toml"
    spec_path.write_text(
        """
release_id = "exact-mixture"
exercise = "assembly"
mode = "practice"

[[samples]]
source_id = "target"
public_id = "Sample_dirty"
implant = "CONTAMINATED"
[samples.implant_parameters]
contaminant_source_id = "ecoli"
contamination_fraction = 0.50
[samples.expected_answers]
species = "Klebsiella pneumoniae"
""".strip()
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "genomepuzzle.read_generation._simulate_short_reads", _fake_short
    )
    output = tmp_path / "release"
    generate_read_release(load_release_spec(spec_path), source, output)
    provenance = json.loads((output / "private/provenance.json").read_text())
    validation = provenance["samples"][0]["provenance"]["validation"]
    assert validation["clean_pairs"] == 20
    assert validation["contaminant_pairs"] == 20
    assert validation["achieved_contamination_fraction"] == 0.5
