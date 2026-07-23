import gzip
import json

from genomepuzzle.read_generation import generate_read_release
from genomepuzzle.release import load_release_spec


def _fake_short(reference, r1, r2, **_):
    for path, mate in ((r1, 1), (r2, 2)):
        with gzip.open(path, "wt") as handle:
            for index in range(20):
                handle.write(
                    "@source-{0}-{1}\nACGT\n+\nIIII\n".format(mate, index)
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
read_fraction = 0.5
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
    assert validation["read_fraction"] == 0.5
