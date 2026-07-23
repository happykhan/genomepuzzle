import json

import pytest

from genomepuzzle.release import (
    ReleaseArtifactSample,
    derive_sample_id,
    load_release_spec,
    resolve_release_samples,
    write_release_manifests,
)


def write_spec(path, sample_blocks, extra=""):
    path.write_text(
        """
schema_version = "1.0"
release_id = "2026-round-1-typing-practice"
exercise = "typing"
mode = "practice"
master_seed = 470
id_salt_env = "GENOMEPUZZLE_ID_SALT"
{extra}
{samples}
""".format(extra=extra, samples="\n".join(sample_blocks)),
        encoding="utf-8",
    )


def test_release_identity_is_stable_and_independent_of_input_order(tmp_path):
    first = tmp_path / "first.toml"
    second = tmp_path / "second.toml"
    samples = [
        '[[samples]]\nsource_id = "GCA_000001.1"',
        '[[samples]]\nsource_id = "GCA_000002.1"\nimplant = "MIXED_CONTIGS"',
    ]
    write_spec(first, samples)
    write_spec(second, list(reversed(samples)))

    first_resolved = {
        row.source_id: (row.sample_id, row.random_seed)
        for row in resolve_release_samples(load_release_spec(first), id_salt="private")
    }
    second_resolved = {
        row.source_id: (row.sample_id, row.random_seed)
        for row in resolve_release_samples(load_release_spec(second), id_salt="private")
    }

    assert first_resolved == second_resolved
    assert "GCA" not in first_resolved["GCA_000001.1"][0]
    assert "MIXED" not in first_resolved["GCA_000002.1"][0]


def test_release_identity_changes_with_release_or_private_salt():
    baseline = derive_sample_id("round-1", "source-a", "salt-a")
    assert baseline != derive_sample_id("round-2", "source-a", "salt-a")
    assert baseline != derive_sample_id("round-1", "source-a", "salt-b")


def test_release_spec_rejects_duplicate_identity_keys(tmp_path):
    path = tmp_path / "duplicate.toml"
    write_spec(
        path,
        [
            '[[samples]]\nsource_id = "source-a"',
            '[[samples]]\nsource_id = "source-a"',
        ],
    )
    with pytest.raises(ValueError, match="duplicate sample identity_key"):
        load_release_spec(path)


def test_explicit_public_ids_do_not_require_a_salt(tmp_path):
    path = tmp_path / "explicit.toml"
    write_spec(
        path,
        ['[[samples]]\nsource_id = "source-a"\npublic_id = "Sample_fixed123"'],
    )
    resolved = resolve_release_samples(load_release_spec(path))
    assert resolved[0].sample_id == "Sample_fixed123"


def test_release_manifests_keep_truth_private(tmp_path):
    spec_path = tmp_path / "release.toml"
    write_spec(
        spec_path,
        ['[[samples]]\nsource_id = "GCA_000001.1"\npublic_id = "Sample_fixed123"'],
    )
    spec = load_release_spec(spec_path)
    release_dir = tmp_path / "release"
    files_dir = release_dir / "public" / "files"
    files_dir.mkdir(parents=True)
    fasta = files_dir / "Sample_fixed123.fasta"
    fasta.write_text(">Sample_fixed123_contig_1\nACGT\n", encoding="utf-8")

    outputs = write_release_manifests(
        release_dir,
        spec,
        [
            ReleaseArtifactSample(
                sample_id="Sample_fixed123",
                source_id="GCA_000001.1",
                random_seed=123,
                files={"assembly": str(fasta)},
                expected_answers={"st": "ST42", "species": "K. pneumoniae"},
                implant="MIXED_CONTIGS",
                implant_parameters={"fraction": 0.1},
                public_metadata={"display_label": "Sample_fixed123"},
                private_provenance={"source_sha256": "abc123"},
            )
        ],
        generator={"git_commit": "deadbeef"},
    )

    public_payload = json.loads(
        (release_dir / "public" / "dataset_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    private_payload = json.loads(
        (release_dir / "private" / "provenance.json").read_text(encoding="utf-8")
    )
    answer_payload = json.loads(
        (release_dir / "private" / "answer_key.json").read_text(encoding="utf-8")
    )
    public_text = json.dumps(public_payload)

    assert "GCA_000001.1" not in public_text
    assert "MIXED_CONTIGS" not in public_text
    assert "ST42" not in public_text
    assert private_payload["samples"][0]["source_id"] == "GCA_000001.1"
    assert answer_payload["samples"][0]["answers"]["st"] == "ST42"
    assert "Sample_fixed123.fasta" in (
        release_dir / "public" / "checksums.sha256"
    ).read_text(encoding="utf-8")
    assert set(outputs) == {
        "public_manifest",
        "private_manifest",
        "answer_key",
        "implant_manifest",
        "checksums",
    }


def test_public_metadata_rejects_private_fields(tmp_path):
    spec_path = tmp_path / "release.toml"
    write_spec(
        spec_path,
        ['[[samples]]\nsource_id = "source-a"\npublic_id = "Sample_fixed123"'],
    )
    release_dir = tmp_path / "release"
    files_dir = release_dir / "public" / "files"
    files_dir.mkdir(parents=True)
    fasta = files_dir / "Sample_fixed123.fasta"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="private field"):
        write_release_manifests(
            release_dir,
            load_release_spec(spec_path),
            [
                ReleaseArtifactSample(
                    sample_id="Sample_fixed123",
                    source_id="source-a",
                    random_seed=1,
                    files={"assembly": str(fasta)},
                    expected_answers={},
                    public_metadata={"nested": {"reference_accession": "secret"}},
                )
            ],
        )


def test_release_manifest_requires_every_specified_sample(tmp_path):
    spec_path = tmp_path / "release.toml"
    write_spec(
        spec_path,
        [
            '[[samples]]\nsource_id = "source-a"\npublic_id = "Sample_fixed123"',
            '[[samples]]\nsource_id = "source-b"\npublic_id = "Sample_fixed456"',
        ],
    )
    release_dir = tmp_path / "release"
    files_dir = release_dir / "public" / "files"
    files_dir.mkdir(parents=True)
    fasta = files_dir / "Sample_fixed123.fasta"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="sample count"):
        write_release_manifests(
            release_dir,
            load_release_spec(spec_path),
            [
                ReleaseArtifactSample(
                    sample_id="Sample_fixed123",
                    source_id="source-a",
                    random_seed=1,
                    files={"assembly": str(fasta)},
                    expected_answers={},
                )
            ],
        )


def test_participant_filename_must_use_public_sample_id(tmp_path):
    spec_path = tmp_path / "release.toml"
    write_spec(
        spec_path,
        ['[[samples]]\nsource_id = "source-a"\npublic_id = "Sample_fixed123"'],
    )
    release_dir = tmp_path / "release"
    files_dir = release_dir / "public" / "files"
    files_dir.mkdir(parents=True)
    fasta = files_dir / "source-a.fasta"
    fasta.write_text(">contig\nACGT\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must begin with sample_id"):
        write_release_manifests(
            release_dir,
            load_release_spec(spec_path),
            [
                ReleaseArtifactSample(
                    sample_id="Sample_fixed123",
                    source_id="source-a",
                    random_seed=1,
                    files={"assembly": str(fasta)},
                    expected_answers={},
                )
            ],
        )
