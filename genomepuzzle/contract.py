"""Website-facing GenomePuzzle release contract and bundle validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


BUNDLE_SCHEMA_VERSION = "2.0"
SAMPLE_ID_FIELDS = {"sample_id", "sample", "id", "public_name"}
PENDING_VALUES = {"pending", "pending_analysis", "pending_kleborate"}
PRIVATE_KEY_NAMES = {
    "answer",
    "answers",
    "answer_key",
    "expected",
    "expected_answers",
    "implant",
    "implant_type",
    "reference_accession",
    "source",
    "source_id",
    "source_name",
    "truth",
}


def _field(
    name: str,
    label: str,
    description: str,
    *,
    required: bool = True,
    identifier: bool = False,
    scored: bool = True,
    scorer: str = "exact",
    aliases: Sequence[str] | None = None,
) -> dict[str, Any]:
    return {
        "name": name,
        "label": label,
        "description": description,
        "type": "string",
        "required": required,
        "identifier": identifier,
        "scored": scored and not identifier,
        "scorer": "identifier" if identifier else scorer,
        "aliases": list(aliases or []),
    }


EXERCISE_CONTRACTS: dict[str, dict[str, Any]] = {
    "typing": {
        "title": "Genotyping",
        "description": "Recover species, sequence type, surface loci and carbapenemases.",
        "instructions": [
            "Analyse every supplied assembly.",
            "Keep the sample_id values and column names unchanged.",
            "Return the completed sample_sheet.csv as a CSV file.",
        ],
        "fields": [
            _field("sample_id", "Sample", "Public sample identifier.", identifier=True),
            _field("species", "Species", "Organism call.", required=False),
            _field("st", "Sequence type", "MLST sequence type."),
            _field("k_locus", "K locus", "Capsule locus.", required=False),
            _field("capsule_type", "Capsule type", "Capsule serotype.", required=False),
            _field("wzi", "wzi", "wzi allele.", required=False),
            _field("o_locus", "O locus", "O-antigen locus.", required=False),
            _field("o_type", "O type", "O-antigen type.", required=False),
            _field(
                "bla_carb",
                "Carbapenemases",
                "Detected carbapenemase genes separated by semicolons.",
                required=False,
                scorer="unordered_list",
                aliases=["kleborate_bla_carb"],
            ),
        ],
    },
    "assembly": {
        "title": "Short-read assembly",
        "description": "Assemble paired short reads and identify troublesome datasets.",
        "instructions": [
            "Assemble every paired-read dataset with a reproducible workflow.",
            "Report the final taxonomic call and QC decision.",
            "Use notes only for concise supporting interpretation.",
        ],
        "fields": [
            _field("sample_id", "Sample", "Public sample identifier.", identifier=True),
            _field("species", "Species", "Final taxonomic call."),
            _field("qc", "QC", "PASS or FAIL.", aliases=["qc_decision"]),
            _field("error", "Problem", "Detected problem type.", required=False),
            _field("notes", "Notes", "Concise interpretation.", required=False, scored=False),
        ],
    },
    "hybrid": {
        "title": "Hybrid assembly",
        "description": "Combine short and long reads and assess the resulting assembly.",
        "instructions": [
            "Use both paired short reads and long reads for every sample.",
            "Report the assembler, taxonomic call and final QC decision.",
            "Do not infer or report the private source accession.",
        ],
        "fields": [
            _field("sample_id", "Sample", "Public sample identifier.", identifier=True),
            _field("species", "Species", "Final taxonomic call."),
            _field("assembler", "Assembler", "Assembly workflow used.", required=False, scored=False),
            _field("qc", "QC", "PASS or FAIL.", aliases=["qc_decision"]),
            _field("error", "Problem", "Detected problem type.", required=False),
            _field("notes", "Notes", "Concise interpretation.", required=False, scored=False),
        ],
    },
    "outbreak": {
        "title": "Phylogeny and outbreak investigation",
        "description": "Infer relatedness and identify epidemiologically plausible clusters.",
        "instructions": [
            "Perform read QC, mapping, variant calling and phylogenetic analysis.",
            "Assign a cluster label to every sample; label names themselves are arbitrary.",
            "Flag samples that should be excluded from the cluster analysis.",
        ],
        "fields": [
            _field("sample_id", "Sample", "Public sample identifier.", identifier=True),
            _field(
                "cluster",
                "Cluster",
                "Inferred outbreak cluster label.",
                scorer="partition",
            ),
            _field("species", "Species", "Final taxonomic call.", required=False),
            _field(
                "qc_decision",
                "QC decision",
                "include or exclude.",
                aliases=["qc"],
            ),
            _field("notes", "Notes", "Concise interpretation.", required=False, scored=False),
        ],
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def json_dump(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / ".{0}.tmp-{1}".format(path.name, os.getpid())
    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_answer_fields(answers: Mapping[str, Any]) -> dict[str, Any]:
    """Translate analyser-specific names into the stable participant contract."""

    aliases = {
        "kleborate_st": "st",
        "sequence_type": "st",
        "tax_classification": "species",
        "error_type": "error",
    }
    normalized: dict[str, Any] = {}
    for raw_name, value in answers.items():
        name = aliases.get(str(raw_name).strip().lower(), str(raw_name).strip().lower())
        if name in normalized and normalized[name] != value:
            raise ValueError("conflicting answer values for field {0}".format(name))
        normalized[name] = value
    return normalized


def exercise_contract(exercise: str) -> dict[str, Any]:
    try:
        contract = EXERCISE_CONTRACTS[exercise]
    except KeyError as exc:
        raise ValueError("unsupported exercise contract: {0}".format(exercise)) from exc
    return json.loads(json.dumps(contract))


def write_participant_contract(
    release_dir: str | os.PathLike[str],
    *,
    release_id: str,
    exercise: str,
    mode: str,
    sample_ids: Sequence[str],
    title: str | None = None,
    description: str | None = None,
    instructions: Sequence[str] | None = None,
    pass_threshold: float = 0.8,
    scored_fields: set[str] | None = None,
) -> dict[str, str]:
    """Write the submission schema, blank answer sheet, instructions and scorer policy."""

    if not 0 < pass_threshold <= 1:
        raise ValueError("pass_threshold must be greater than 0 and at most 1")
    contract = exercise_contract(exercise)
    public_dir = Path(release_dir) / "public"
    private_dir = Path(release_dir) / "private"
    public_dir.mkdir(parents=True, exist_ok=True)
    private_dir.mkdir(parents=True, exist_ok=True)
    fields = contract["fields"]

    submission_schema = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": release_id,
        "exercise": exercise,
        "mode": mode,
        "fields": fields,
    }
    submission_path = public_dir / "submission_schema.json"
    json_dump(submission_path, submission_schema)

    sheet_path = public_dir / "sample_sheet.csv"
    names = [field["name"] for field in fields]
    with open(sheet_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        for sample_id in sample_ids:
            writer.writerow({"sample_id": sample_id})

    final_title = title or contract["title"]
    final_description = description or contract["description"]
    final_instructions = list(instructions or contract["instructions"])
    instructions_path = public_dir / "instructions.md"
    with open(instructions_path, "w", encoding="utf-8") as handle:
        handle.write("# {0}\n\n{1}\n\n".format(final_title, final_description))
        for index, item in enumerate(final_instructions, start=1):
            handle.write("{0}. {1}\n".format(index, item))

    effective_scored_fields = scored_fields
    scoring_policy = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": release_id,
        "scorer_version": "2.0",
        "pass_threshold": pass_threshold,
        "require_all_samples": True,
        "reject_unexpected_samples": True,
        "manual_review": {
            "enabled": True,
            "flag_on_parse_error": True,
            "flag_within_pass_margin": 0.05,
        },
        "fields": [
            {
                "name": field["name"],
                "scored": field["scored"]
                and (
                    effective_scored_fields is None
                    or field["name"] in effective_scored_fields
                ),
                "scorer": field["scorer"],
                "weight": 1.0,
                "aliases": field["aliases"],
            }
            for field in fields
            if not field["identifier"]
        ],
    }
    scoring_path = private_dir / "scoring_policy.json"
    json_dump(scoring_path, scoring_policy)
    return {
        "submission_schema": str(submission_path),
        "sample_sheet": str(sheet_path),
        "instructions": str(instructions_path),
        "scoring_policy": str(scoring_path),
    }


def write_release_index(
    release_dir: str | os.PathLike[str],
    *,
    release_id: str,
    exercise: str,
    mode: str,
    title: str,
    description: str,
) -> Path:
    path = Path(release_dir) / "release.json"
    payload = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": release_id,
        "exercise": exercise,
        "mode": mode,
        "title": title,
        "description": description,
        "artifacts": {
            "public_manifest": "public/manifest.json",
            "legacy_public_manifest": "public/dataset_manifest.json",
            "sample_sheet": "public/sample_sheet.csv",
            "submission_schema": "public/submission_schema.json",
            "instructions": "public/instructions.md",
            "checksums": "public/checksums.sha256",
            "answer_key": "private/answer_key.json",
            "scoring_policy": "private/scoring_policy.json",
            "provenance": "private/provenance.json",
            "implant_manifest": "private/implant_manifest.json",
            "validation_report": "private/validation_report.json",
        },
    }
    json_dump(path, payload)
    return path


def _load_json(path: Path) -> Any:
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("invalid JSON artifact {0}: {1}".format(path, exc)) from exc


def _assert_no_private_keys(value: Any, path: str = "public") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if str(key).lower() in PRIVATE_KEY_NAMES:
                raise ValueError("private key in public artifact: {0}.{1}".format(path, key))
            _assert_no_private_keys(item, "{0}.{1}".format(path, key))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _assert_no_private_keys(item, "{0}[{1}]".format(path, index))


def _read_first_headers(path: Path, limit: int = 100) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    headers: list[str] = []
    mode = "rt"
    try:
        with opener(path, mode, encoding="utf-8", errors="replace") as handle:
            if path.name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
                for index, line in enumerate(handle):
                    if index % 4 == 0:
                        headers.append(line.strip())
                        if len(headers) >= limit:
                            break
            elif path.name.endswith((".fasta", ".fa", ".fna")):
                for line in handle:
                    if line.startswith(">"):
                        headers.append(line.strip())
                        if len(headers) >= limit:
                            break
    except (OSError, EOFError) as exc:
        raise ValueError("participant sequence file is unreadable: {0}: {1}".format(path, exc)) from exc
    return headers


def _bundle_digest(release_dir: Path) -> str:
    digest = hashlib.sha256()
    excluded = {"COMPLETE", "COMPLETE.json"}
    for path in sorted(item for item in release_dir.rglob("*") if item.is_file()):
        if path.name in excluded:
            continue
        relative = path.relative_to(release_dir).as_posix()
        if relative.startswith("build/"):
            continue
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def validate_release_bundle(
    release_dir: str | os.PathLike[str],
    *,
    require_complete: bool = False,
    write_report: bool = False,
) -> dict[str, Any]:
    """Validate the complete website/private release contract and participant files."""

    root = Path(release_dir).resolve()
    required = [
        root / "release.json",
        root / "public" / "manifest.json",
        root / "public" / "dataset_manifest.json",
        root / "public" / "sample_sheet.csv",
        root / "public" / "submission_schema.json",
        root / "public" / "instructions.md",
        root / "public" / "checksums.sha256",
        root / "private" / "answer_key.json",
        root / "private" / "scoring_policy.json",
        root / "private" / "provenance.json",
        root / "private" / "implant_manifest.json",
    ]
    missing = [str(path.relative_to(root)) for path in required if not path.is_file()]
    if missing:
        raise ValueError("release bundle is missing: {0}".format(", ".join(missing)))

    index = _load_json(root / "release.json")
    manifest = _load_json(root / "public" / "manifest.json")
    legacy_manifest = _load_json(root / "public" / "dataset_manifest.json")
    schema = _load_json(root / "public" / "submission_schema.json")
    answers = _load_json(root / "private" / "answer_key.json")
    policy = _load_json(root / "private" / "scoring_policy.json")
    provenance = _load_json(root / "private" / "provenance.json")
    for payload_name, payload in (
        ("release", index),
        ("manifest", manifest),
        ("submission schema", schema),
        ("answer key", answers),
        ("scoring policy", policy),
        ("provenance", provenance),
    ):
        if not isinstance(payload, dict):
            raise ValueError("{0} must be a JSON object".format(payload_name))
    for payload in (index, manifest, schema, answers, policy, provenance):
        if payload.get("release_id") != index.get("release_id"):
            raise ValueError("release_id differs across bundle artifacts")
    if index.get("schema_version") != BUNDLE_SCHEMA_VERSION:
        raise ValueError("unsupported release bundle schema version")
    if manifest != legacy_manifest:
        raise ValueError("manifest.json and dataset_manifest.json differ")
    _assert_no_private_keys(manifest)

    sample_rows = manifest.get("samples")
    answer_rows = answers.get("samples")
    fields = schema.get("fields")
    if not isinstance(sample_rows, list) or not sample_rows:
        raise ValueError("public manifest has no samples")
    if not isinstance(answer_rows, list) or not isinstance(fields, list):
        raise ValueError("invalid answer key or submission schema")
    sample_ids = [row.get("sample_id") for row in sample_rows]
    if any(not isinstance(value, str) or not value for value in sample_ids):
        raise ValueError("public manifest contains an invalid sample_id")
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("public manifest contains duplicate sample IDs")
    answer_ids = [row.get("sample_id") for row in answer_rows]
    if answer_ids != sample_ids:
        raise ValueError("answer key sample order or membership differs from manifest")

    field_names = [field.get("name") for field in fields]
    if not field_names or field_names[0] != "sample_id" or len(set(field_names)) != len(field_names):
        raise ValueError("submission schema requires unique fields beginning with sample_id")
    scored_fields = {
        item.get("name") for item in policy.get("fields", []) if item.get("scored")
    }
    for row in answer_rows:
        row_answers = row.get("answers")
        if not isinstance(row_answers, dict):
            raise ValueError("answer row must contain an answers object")
        unknown = set(row_answers) - set(field_names) - {"analysis_status"}
        if unknown:
            raise ValueError(
                "answer key contains fields outside submission schema: {0}".format(
                    ", ".join(sorted(unknown))
                )
            )
        pending = {str(value).strip().lower() for value in row_answers.values()} & PENDING_VALUES
        if pending:
            raise ValueError("answer key contains pending reference analysis")
        absent = scored_fields - set(row_answers)
        if absent:
            raise ValueError(
                "answer key is missing scored fields for {0}: {1}".format(
                    row["sample_id"], ", ".join(sorted(absent))
                )
            )

    with open(root / "public" / "sample_sheet.csv", encoding="utf-8", newline="") as handle:
        sheet_rows = list(csv.DictReader(handle))
    if not sheet_rows or list(sheet_rows[0].keys()) != field_names:
        raise ValueError("sample_sheet.csv columns differ from submission schema")
    if [row.get("sample_id") for row in sheet_rows] != sample_ids:
        raise ValueError("sample_sheet.csv samples differ from manifest")

    provenance_sources = {
        str(row.get("source_id", ""))
        for row in provenance.get("samples", [])
        if row.get("source_id")
    }
    for row in provenance.get("samples", []):
        validation = row.get("provenance", {}).get("validation", {})
        if validation.get("status") != "passed":
            raise ValueError(
                "sample {0} has no passing implant validation".format(
                    row.get("sample_id", "<unknown>")
                )
            )
    file_count = 0
    for sample in sample_rows:
        files = sample.get("files")
        if not isinstance(files, dict) or not files:
            raise ValueError("sample {0} has no participant files".format(sample["sample_id"]))
        for details in files.values():
            path = root / "public" / "files" / str(details.get("filename", ""))
            if not path.is_file():
                raise ValueError("missing participant file: {0}".format(path.name))
            if details.get("size") != path.stat().st_size:
                raise ValueError("size mismatch for {0}".format(path.name))
            if details.get("sha256") != sha256_file(path):
                raise ValueError("checksum mismatch for {0}".format(path.name))
            for header in _read_first_headers(path):
                if sample["sample_id"] not in header:
                    raise ValueError("non-anonymous sequence header in {0}".format(path.name))
                anonymous_part = header.replace(sample["sample_id"], "")
                leaked = [
                    source
                    for source in provenance_sources
                    if source and source in anonymous_part
                ]
                if leaked:
                    raise ValueError("source identity leaked in {0} header".format(path.name))
            file_count += 1

    if require_complete:
        complete_path = root / "COMPLETE.json"
        if not complete_path.is_file():
            raise ValueError("release is not complete")
        complete = _load_json(complete_path)
        if complete.get("bundle_sha256") != _bundle_digest(root):
            raise ValueError("completed release digest does not match bundle")

    report = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": index["release_id"],
        "status": "passed",
        "validated_at": utc_now(),
        "checks": {
            "contract": "passed",
            "checksums": "passed",
            "answers": "passed",
            "privacy_headers": "passed",
        },
        "sample_count": len(sample_ids),
        "participant_file_count": file_count,
    }
    if write_report:
        json_dump(root / "private" / "validation_report.json", report)
    return report


def complete_release(release_dir: str | os.PathLike[str]) -> Path:
    """Validate and seal a release with a content digest."""

    root = Path(release_dir).resolve()
    report = validate_release_bundle(root, write_report=True)
    complete_path = root / "COMPLETE.json"
    payload = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": report["release_id"],
        "status": "complete",
        "validated_at": report["validated_at"],
        "bundle_sha256": _bundle_digest(root),
    }
    json_dump(complete_path, payload)
    legacy = root / "COMPLETE"
    if legacy.exists():
        legacy.unlink()
    return complete_path


def inspect_release(release_dir: str | os.PathLike[str]) -> dict[str, Any]:
    root = Path(release_dir).resolve()
    report = validate_release_bundle(root, require_complete=(root / "COMPLETE.json").is_file())
    index = _load_json(root / "release.json")
    implants = _load_json(root / "private" / "implant_manifest.json")
    complete = _load_json(root / "COMPLETE.json") if (root / "COMPLETE.json").is_file() else None
    return {
        "release_id": index["release_id"],
        "exercise": index["exercise"],
        "mode": index["mode"],
        "status": "complete" if complete else "unsealed",
        "samples": report["sample_count"],
        "participant_files": report["participant_file_count"],
        "implants": implants.get("summary", {}),
        "bundle_sha256": complete.get("bundle_sha256") if complete else None,
    }
