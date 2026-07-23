import json
from pathlib import Path

import pytest

from genomepuzzle.release import load_release_spec
from genomepuzzle.typing import build_typing_release, fragment_records, read_fasta


def _write_spec(path: Path) -> None:
    path.write_text(
        """
schema_version = "1.0"
release_id = "typing-2026"
exercise = "typing"
mode = "practice"
master_seed = 17

[[samples]]
source_id = "clean"
public_id = "Sample_clean"

[[samples]]
source_id = "target"
public_id = "Sample_mixed"
implant = "MIXED_CONTIGS"
[samples.implant_parameters]
contaminant_source_id = "contaminant"
contamination_fraction = 0.25

[[samples]]
source_id = "target"
identity_key = "target-fragmented"
public_id = "Sample_fragmented"
implant = "FRAGMENTED"
[samples.implant_parameters]
fragment_size = 100
""".strip()
        + "\n",
        encoding="utf-8",
    )


def test_fragment_records_validates_size():
    with pytest.raises(ValueError, match="at least 100"):
        fragment_records([("one", "ACGT")], 20)


def test_build_typing_release_is_anonymous_and_tracks_truth(tmp_path):
    sources = tmp_path / "sources"
    sources.mkdir()
    (sources / "clean.fasta").write_text(">ACCESSION secret\n" + "A" * 250 + "\n")
    (sources / "target.fasta").write_text(">target\n" + "C" * 250 + "\n")
    (sources / "contaminant.fasta").write_text(">contaminant\n" + "G" * 100 + "\n")
    spec_path = tmp_path / "typing.toml"
    _write_spec(spec_path)
    release_dir = tmp_path / "release"

    manifests = build_typing_release(
        load_release_spec(spec_path),
        sources,
        release_dir,
        analyser=lambda path: {"analysis_status": "complete", "st": path.stem},
    )

    public = json.loads(Path(manifests["public_manifest"]).read_text())
    private = json.loads(Path(manifests["private_manifest"]).read_text())
    assert [row["sample_id"] for row in public["samples"]] == [
        "Sample_clean",
        "Sample_mixed",
        "Sample_fragmented",
    ]
    assert "source_id" not in json.dumps(public)
    assert "contaminant" not in json.dumps(public)
    assert private["samples"][1]["implant"]["type"] == "MIXED_CONTIGS"
    assert private["samples"][1]["provenance"]["contaminant_source_file"].endswith(
        "contaminant.fasta"
    )
    assert (release_dir / "COMPLETE.json").is_file()

    clean_records = read_fasta(release_dir / "public/files/Sample_clean.fasta")
    mixed_records = read_fasta(release_dir / "public/files/Sample_mixed.fasta")
    fragmented_records = read_fasta(
        release_dir / "public/files/Sample_fragmented.fasta"
    )
    assert all(name.startswith("Sample_clean_contig_") for name, _ in clean_records)
    assert sum(len(seq) for _, seq in mixed_records) == 312
    mixed_validation = private["samples"][1]["provenance"]["validation"]
    assert mixed_validation["contaminant_bases"] == 62
    assert mixed_validation["achieved_contamination_fraction"] == 0.248
    assert [len(seq) for _, seq in fragmented_records] == [100, 100, 50]


def test_build_typing_release_rejects_nonempty_destination(tmp_path):
    sources = tmp_path / "sources"
    sources.mkdir()
    (sources / "clean.fasta").write_text(">a\n" + "A" * 100 + "\n")
    (sources / "target.fasta").write_text(">b\n" + "C" * 100 + "\n")
    (sources / "contaminant.fasta").write_text(">c\n" + "G" * 100 + "\n")
    spec_path = tmp_path / "typing.toml"
    _write_spec(spec_path)
    release_dir = tmp_path / "release"
    release_dir.mkdir()
    (release_dir / "existing").write_text("do not overwrite")
    with pytest.raises(ValueError, match="not empty"):
        build_typing_release(
            load_release_spec(spec_path), sources, release_dir, analyser=lambda _: {}
        )
