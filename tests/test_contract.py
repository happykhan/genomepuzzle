import json
from importlib.resources import files

import pytest

from genomepuzzle.contract import (
    BUNDLE_SCHEMA_VERSION,
    exercise_contract,
    failure_reason_for_implant,
    inspect_release,
    retained_read_pairs_for_fault,
    validate_release_bundle,
)
from genomepuzzle.release import load_release_spec
from genomepuzzle.typing import build_typing_release


def _typing_release(tmp_path, analyser):
    source = tmp_path / "sources"
    source.mkdir()
    (source / "source.fasta").write_text(">private_accession\n" + "A" * 200 + "\n")
    spec_path = tmp_path / "release.toml"
    spec_path.write_text(
        """
release_id = "contract-test"
exercise = "typing"
mode = "practice"

[[samples]]
source_id = "source"
public_id = "Sample_public"
""".strip()
        + "\n",
        encoding="utf-8",
    )
    output = tmp_path / "release"
    build_typing_release(
        load_release_spec(spec_path), source, output, analyser=analyser
    )
    return output


def test_release_contract_is_self_describing_and_canonical(tmp_path):
    output = _typing_release(
        tmp_path,
        lambda _: {
            "analysis_status": "complete",
            "kleborate_st": "ST42",
            "species": "Klebsiella pneumoniae",
        },
    )
    index = json.loads((output / "release.json").read_text())
    answers = json.loads((output / "private/answer_key.json").read_text())
    schema = json.loads((output / "public/submission_schema.json").read_text())
    sheet_header = (output / "public/sample_sheet.csv").read_text().splitlines()[0]

    assert index["schema_version"] == BUNDLE_SCHEMA_VERSION
    assert answers["samples"][0]["answers"]["st"] == "42"
    assert "kleborate_st" not in json.dumps(answers)
    assert schema["fields"][0]["name"] == "sample_id"
    assert [field["name"] for field in schema["fields"][:3]] == [
        "sample_id",
        "qc_status",
        "failure_reason",
    ]
    assert answers["samples"][0]["answers"]["qc_status"] == "PASS"
    assert answers["samples"][0]["answers"]["failure_reason"] == "NONE"
    assert sheet_header == ",".join(field["name"] for field in schema["fields"])
    assert inspect_release(output)["status"] == "complete"


def test_completed_bundle_detects_tampering(tmp_path):
    output = _typing_release(
        tmp_path, lambda _: {"analysis_status": "complete", "st": "ST42"}
    )
    with open(output / "public/instructions.md", "a", encoding="utf-8") as handle:
        handle.write("tampered\n")
    with pytest.raises(ValueError, match="digest"):
        validate_release_bundle(output, require_complete=True)


def test_pending_analysis_is_not_sealed(tmp_path):
    output = _typing_release(tmp_path, analyser=None)
    assert not (output / "COMPLETE.json").exists()
    with pytest.raises(ValueError, match="pending"):
        validate_release_bundle(output)


def test_packaged_json_schemas_are_valid_json():
    schema_root = files("genomepuzzle.schemas")
    for name in ("release.schema.json", "manifest.schema.json", "submission.schema.json"):
        payload = json.loads(schema_root.joinpath(name).read_text(encoding="utf-8"))
        assert payload["$schema"].endswith("2020-12/schema")


def test_typing_scoring_is_conditional_on_qc(tmp_path):
    output = _typing_release(
        tmp_path,
        lambda _: {
            "analysis_status": "complete",
            "st": "ST42",
            "species": "Klebsiella pneumoniae",
        },
    )
    policy = json.loads((output / "private/scoring_policy.json").read_text())
    fields = {field["name"]: field for field in policy["fields"]}
    assert fields["qc_status"]["score_when"] is None
    assert fields["failure_reason"]["score_when"] is None
    assert fields["species"]["score_when"] == {
        "field": "qc_status",
        "equals": "PASS",
    }


def test_exercises_use_relevant_canonical_columns():
    columns = {
        exercise: [field["name"] for field in exercise_contract(exercise)["fields"]]
        for exercise in ("typing", "assembly", "hybrid", "outbreak")
    }
    assert columns["typing"] == [
        "sample_id",
        "qc_status",
        "failure_reason",
        "species",
        "st",
        "k_locus",
        "capsule_type",
        "wzi",
        "o_locus",
        "o_type",
        "bla_carb",
        "notes",
    ]
    assembly_columns = [
        "sample_id",
        "qc_status",
        "failure_reason",
        "species",
        "assembler",
        "contig_count",
        "total_length",
        "n50",
        "longest_contig",
        "notes",
    ]
    assert columns["assembly"] == assembly_columns
    assert columns["hybrid"] == assembly_columns
    assert columns["outbreak"] == [
        "sample_id",
        "qc_status",
        "failure_reason",
        "species",
        "cluster",
        "notes",
    ]


def test_private_faults_map_to_small_public_vocabulary():
    assert (
        failure_reason_for_implant("assembly", "TRUNCATED_R2_TO_10_READS")
        == "TOO_FEW_READS"
    )
    assert failure_reason_for_implant("typing", "FRAGMENTED") == "EXTREME_FRAGMENTATION"
    assert failure_reason_for_implant("outbreak", "NORMAL") == "NONE"
    with pytest.raises(ValueError, match="no participant failure reason"):
        failure_reason_for_implant("assembly", "SUBTLE_MATE_COUNT_DIFFERENCE")


def test_catastrophic_read_count_fault_is_explicit_and_bounded():
    assert retained_read_pairs_for_fault("TEN_READ_PAIRS", {}) == 10
    assert (
        retained_read_pairs_for_fault(
            "TRUNCATE_TO_READ_PAIRS", {"retained_read_pairs": 37}
        )
        == 37
    )
    with pytest.raises(ValueError, match="requires integer"):
        retained_read_pairs_for_fault("TRUNCATE_TO_READ_PAIRS", {})
    with pytest.raises(ValueError, match="between 1 and 100"):
        retained_read_pairs_for_fault(
            "TRUNCATE_TO_READ_PAIRS", {"retained_read_pairs": 101}
        )
