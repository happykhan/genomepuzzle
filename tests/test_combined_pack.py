import gzip
import json
from pathlib import Path

import pytest

from genomepuzzle.combined_pack import (
    build_combined_pack,
    load_combined_pack_spec,
    validate_combined_pack,
)
from genomepuzzle.reads_release import package_read_release
from genomepuzzle.release import load_release_spec


def _component_release(tmp_path: Path, exercise: str) -> Path:
    source = tmp_path / "{0}-source".format(exercise)
    source.mkdir()
    for mate in (1, 2):
        with gzip.open(source / "source-a_R{0}.fastq.gz".format(mate), "wt") as handle:
            handle.write("@source-a/{0}\nACGT\n+\nIIII\n".format(mate))
    if exercise == "hybrid":
        with gzip.open(source / "source-a_long.fastq.gz", "wt") as handle:
            handle.write("@source-a-long\nACGT\n+\nIIII\n")
    spec_path = tmp_path / "{0}.toml".format(exercise)
    spec_path.write_text(
        """
release_id = "challenge-2-{exercise}"
exercise = "{exercise}"
mode = "challenge"

[[samples]]
source_id = "source-a"
public_id = "Sample_{exercise}"
""".format(exercise=exercise).strip()
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "{0}-release".format(exercise)
    package_read_release(
        load_release_spec(spec_path),
        source,
        {"source-a": {"qc": "pass", "species": "Klebsiella pneumoniae"}},
        output,
    )
    return output


def _pack_spec(
    tmp_path: Path,
    assembly_release: Path,
    hybrid_release: Path,
) -> Path:
    assembly_digest = json.loads(
        (assembly_release / "COMPLETE.json").read_text(encoding="utf-8")
    )["bundle_sha256"]
    hybrid_digest = json.loads(
        (hybrid_release / "COMPLETE.json").read_text(encoding="utf-8")
    )["bundle_sha256"]
    spec = tmp_path / "combined.toml"
    spec.write_text(
        """
schema_version = "1.0"
pack_id = "combined-eqa"
title = "Combined assembly EQA"
description = "One participant exercise with two evidence tracks."
mode = "challenge"
release_date = "2026-08-14"
instructions = [
  "Complete both tracks.",
  "Return both result sheets.",
]

[[tracks]]
name = "short-read"
release_id = "challenge-2-assembly"
exercise = "assembly"
bundle_sha256 = "{assembly_digest}"

[[tracks]]
name = "long-read"
release_id = "challenge-2-hybrid"
exercise = "hybrid"
bundle_sha256 = "{hybrid_digest}"
""".format(
            assembly_digest=assembly_digest,
            hybrid_digest=hybrid_digest,
        ).strip()
        + "\n",
        encoding="utf-8",
    )
    return spec


def test_build_combined_pack_keeps_tracks_and_private_evidence_separate(tmp_path):
    assembly = _component_release(tmp_path, "assembly")
    hybrid = _component_release(tmp_path, "hybrid")
    spec = _pack_spec(tmp_path, assembly, hybrid)

    output = build_combined_pack(
        spec,
        short_read_release=assembly,
        long_read_release=hybrid,
        output_dir=tmp_path / "combined",
    )
    report = validate_combined_pack(output)
    instructions = (output / "public/instructions.md").read_text(encoding="utf-8")

    assert report["samples"] == {"short-read": 1, "long-read": 1}
    assert (output / "public/short-read/sample_sheet.csv").is_file()
    assert (output / "public/long-read/sample_sheet.csv").is_file()
    assert (output / "private/short-read/answer_key.json").is_file()
    assert (output / "private/long-read/answer_key.json").is_file()
    assert not list((output / "public").rglob("private"))
    assert "Submit both result sheets" in instructions
    assert "hybrid-assembly task" in instructions
    assert "zero-byte file can both be intentional" in instructions


def test_combined_pack_rejects_unexpected_component_digest(tmp_path):
    assembly = _component_release(tmp_path, "assembly")
    hybrid = _component_release(tmp_path, "hybrid")
    spec = _pack_spec(tmp_path, assembly, hybrid)
    text = spec.read_text(encoding="utf-8")
    assembly_digest = json.loads(
        (assembly / "COMPLETE.json").read_text(encoding="utf-8")
    )["bundle_sha256"]
    spec.write_text(
        text.replace(assembly_digest, "0" * 64),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="short-read bundle digest"):
        build_combined_pack(
            spec,
            short_read_release=assembly,
            long_read_release=hybrid,
            output_dir=tmp_path / "combined",
        )


def test_combined_pack_detects_public_tampering(tmp_path):
    assembly = _component_release(tmp_path, "assembly")
    hybrid = _component_release(tmp_path, "hybrid")
    output = build_combined_pack(
        _pack_spec(tmp_path, assembly, hybrid),
        short_read_release=assembly,
        long_read_release=hybrid,
        output_dir=tmp_path / "combined",
    )
    with open(output / "public/instructions.md", "a", encoding="utf-8") as handle:
        handle.write("changed\n")

    with pytest.raises(ValueError, match="digest"):
        validate_combined_pack(output)


def test_combined_pack_spec_requires_exact_track_pair(tmp_path):
    spec = tmp_path / "invalid.toml"
    spec.write_text(
        """
schema_version = "1.0"
pack_id = "combined-eqa"
title = "Combined EQA"
description = "Missing the hybrid track."
mode = "challenge"
release_date = "2026-08-14"
instructions = ["Complete the exercise."]

[[tracks]]
name = "short-read"
release_id = "challenge-2-assembly"
exercise = "assembly"
bundle_sha256 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
""".strip()
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="requires one short-read"):
        load_combined_pack_spec(spec)
