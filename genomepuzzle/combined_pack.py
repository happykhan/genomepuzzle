"""Build a combined short-read and hybrid participant release pack."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    import tomli as tomllib

from genomepuzzle.contract import (
    BUNDLE_SCHEMA_VERSION,
    json_dump,
    sha256_file,
    utc_now,
    validate_release_bundle,
)


PACK_SCHEMA_VERSION = "1.0"
TRACKS = {
    "short-read": "assembly",
    "long-read": "hybrid",
}


@dataclass(frozen=True)
class TrackSpec:
    name: str
    release_id: str
    exercise: str
    bundle_sha256: str


@dataclass(frozen=True)
class CombinedPackSpec:
    pack_id: str
    title: str
    description: str
    mode: str
    release_date: str
    instructions: tuple[str, ...]
    tracks: tuple[TrackSpec, ...]


def _required_string(payload: Mapping[str, Any], name: str) -> str:
    value = payload.get(name)
    if not isinstance(value, str) or not value.strip():
        raise ValueError("{0} must be a non-empty string".format(name))
    return value.strip()


def load_combined_pack_spec(path: str | os.PathLike[str]) -> CombinedPackSpec:
    """Load and validate a combined EQA pack specification."""

    with open(path, "rb") as handle:
        payload = tomllib.load(handle)
    if payload.get("schema_version") != PACK_SCHEMA_VERSION:
        raise ValueError("unsupported combined pack schema version")
    mode = _required_string(payload, "mode")
    if mode not in {"practice", "challenge"}:
        raise ValueError("mode must be practice or challenge")
    raw_instructions = payload.get("instructions")
    if (
        not isinstance(raw_instructions, list)
        or not raw_instructions
        or any(not isinstance(item, str) or not item.strip() for item in raw_instructions)
    ):
        raise ValueError("instructions must be a non-empty list of strings")
    raw_tracks = payload.get("tracks")
    if not isinstance(raw_tracks, list):
        raise ValueError("tracks must be an array of tables")
    tracks = tuple(
        TrackSpec(
            name=_required_string(track, "name"),
            release_id=_required_string(track, "release_id"),
            exercise=_required_string(track, "exercise"),
            bundle_sha256=_required_string(track, "bundle_sha256"),
        )
        for track in raw_tracks
        if isinstance(track, Mapping)
    )
    if {track.name: track.exercise for track in tracks} != TRACKS:
        raise ValueError(
            "combined pack requires one short-read assembly track and one long-read hybrid track"
        )
    if any(
        len(track.bundle_sha256) != 64
        or any(character not in "0123456789abcdef" for character in track.bundle_sha256)
        for track in tracks
    ):
        raise ValueError("track bundle_sha256 values must be lower-case SHA-256 digests")
    return CombinedPackSpec(
        pack_id=_required_string(payload, "pack_id"),
        title=_required_string(payload, "title"),
        description=_required_string(payload, "description"),
        mode=mode,
        release_date=_required_string(payload, "release_date"),
        instructions=tuple(item.strip() for item in raw_instructions),
        tracks=tracks,
    )


def _load_json(path: Path) -> dict[str, Any]:
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("invalid JSON artifact {0}: {1}".format(path, exc)) from exc
    if not isinstance(payload, dict):
        raise ValueError("JSON artifact must be an object: {0}".format(path))
    return payload


def _tree_digest(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path == root / "COMPLETE.json":
            continue
        relative = path.relative_to(root).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _submission_summary(track_dir: Path) -> dict[str, Any]:
    schema = _load_json(track_dir / "submission_schema.json")
    fields = schema.get("fields")
    if not isinstance(fields, list):
        raise ValueError("track submission schema has no fields")
    return {
        "columns": [field.get("name") for field in fields],
        "allowed_values": {
            str(field["name"]): field["allowed_values"]
            for field in fields
            if isinstance(field, Mapping) and field.get("allowed_values")
        },
        "scored_columns": [
            field.get("name")
            for field in fields
            if isinstance(field, Mapping) and field.get("scored")
        ],
    }


def _write_participant_instructions(
    path: Path,
    spec: CombinedPackSpec,
    track_details: Mapping[str, Mapping[str, Any]],
) -> None:
    lines = [
        "# {0}".format(spec.title),
        "",
        spec.description,
        "",
        "**Release date:** {0}".format(spec.release_date),
        "",
        "## What to submit",
        "",
        "This is one EQA exercise with two separately assessed tracks. Submit both result sheets; "
        "do not merge their rows or columns into one file.",
        "",
        "| Track | Input evidence | Result template |",
        "| --- | --- | --- |",
        "| Short-read | Paired short-read FASTQ files | `short-read/sample_sheet.csv` |",
        "| Long-read | Paired short-read plus long-read FASTQ files | `long-read/sample_sheet.csv` |",
        "",
        "The long-read track is a hybrid-assembly task. Use and compare both sequencing modalities; "
        "it is not a long-read-only assembly.",
        "",
        "## Instructions",
        "",
    ]
    lines.extend(
        "{0}. {1}".format(index, instruction)
        for index, instruction in enumerate(spec.instructions, start=1)
    )
    lines.extend(
        [
            "",
            "## Expected outputs",
            "",
            "Keep every supplied `sample_id` and header unchanged. Return UTF-8 CSV files with one "
            "row per sample. For `PASS`, complete the analytical fields available from that track. "
            "For `FAIL`, select one permitted categorical reason; unavailable analytical fields may "
            "remain blank.",
            "",
        ]
    )
    for track_name in TRACKS:
        detail = track_details[track_name]
        summary = detail["submission"]
        lines.extend(
            [
                "### {0}".format(
                    "Short-read result sheet"
                    if track_name == "short-read"
                    else "Long-read result sheet"
                ),
                "",
                "Columns: `{0}`".format("`, `".join(summary["columns"])),
                "",
                "Scored columns: `{0}`".format("`, `".join(summary["scored_columns"])),
                "",
            ]
        )
        for field, values in summary["allowed_values"].items():
            lines.append(
                "Allowed `{0}` values: `{1}`".format(field, "`, `".join(values))
            )
        lines.append("")
    lines.extend(
        [
            "## Important notes",
            "",
            "- Sample identifiers are track-specific. Do not assume a short-read sample corresponds "
            "to a long-read sample with a similar position or name.",
            "- A missing file link and a supplied zero-byte file can both be intentional assessment "
            "inputs. Report the observed evidence; do not invent or repair a missing input.",
            "- `qc_status` must be `PASS` or `FAIL`. Use `NONE` as the failure reason only for a "
            "passing sample.",
            "- Record the assembler or workflow and version where requested. These fields support "
            "review even when they are not automatically scored.",
            "- Keep the two completed result sheets. They are separate evidence for one combined "
            "participant exercise.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_public_checksums(public_dir: Path) -> None:
    checksum_path = public_dir / "checksums.sha256"
    rows = []
    for path in sorted(item for item in public_dir.rglob("*") if item.is_file()):
        if path == checksum_path:
            continue
        rows.append(
            "{0}  {1}".format(
                sha256_file(path),
                path.relative_to(public_dir).as_posix(),
            )
        )
    checksum_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def build_combined_pack(
    spec_file: str | os.PathLike[str],
    *,
    short_read_release: str | os.PathLike[str],
    long_read_release: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
) -> Path:
    """Combine two sealed component releases without mixing their evidence."""

    spec = load_combined_pack_spec(spec_file)
    destination = Path(output_dir).resolve()
    if destination.exists() and any(destination.iterdir()):
        raise ValueError("combined pack output directory is not empty")
    sources = {
        "short-read": Path(short_read_release).resolve(),
        "long-read": Path(long_read_release).resolve(),
    }
    track_specs = {track.name: track for track in spec.tracks}
    track_details: dict[str, dict[str, Any]] = {}
    for track_name, exercise in TRACKS.items():
        source = sources[track_name]
        report = validate_release_bundle(source, require_complete=True)
        index = _load_json(source / "release.json")
        complete = _load_json(source / "COMPLETE.json")
        expected = track_specs[track_name]
        if index.get("release_id") != expected.release_id:
            raise ValueError(
                "{0} release_id does not match pack specification".format(track_name)
            )
        if index.get("exercise") != exercise or index.get("mode") != spec.mode:
            raise ValueError("{0} release exercise or mode does not match pack".format(track_name))
        if complete.get("bundle_sha256") != expected.bundle_sha256:
            raise ValueError("{0} bundle digest does not match pack specification".format(track_name))
        track_details[track_name] = {
            "release_id": expected.release_id,
            "exercise": exercise,
            "bundle_sha256": expected.bundle_sha256,
            "sample_count": report["sample_count"],
            "participant_file_count": report["participant_file_count"],
            "submission": _submission_summary(source / "public"),
        }

    destination.mkdir(parents=True, exist_ok=True)
    public_dir = destination / "public"
    private_dir = destination / "private"
    for track_name, source in sources.items():
        shutil.copytree(source / "public", public_dir / track_name)
        shutil.copytree(source / "private", private_dir / track_name)
        shutil.copy2(source / "release.json", private_dir / track_name / "release.json")
        shutil.copy2(source / "COMPLETE.json", private_dir / track_name / "COMPLETE.json")

    _write_participant_instructions(
        public_dir / "instructions.md",
        spec,
        track_details,
    )
    public_index = {
        "schema_version": PACK_SCHEMA_VERSION,
        "pack_id": spec.pack_id,
        "title": spec.title,
        "description": spec.description,
        "mode": spec.mode,
        "release_date": spec.release_date,
        "tracks": [
            {
                "name": name,
                "release_id": track_details[name]["release_id"],
                "exercise": track_details[name]["exercise"],
                "sample_count": track_details[name]["sample_count"],
                "sample_sheet": "{0}/sample_sheet.csv".format(name),
                "submission_schema": "{0}/submission_schema.json".format(name),
            }
            for name in TRACKS
        ],
    }
    json_dump(public_dir / "pack.json", public_index)
    _write_public_checksums(public_dir)
    pack_index = {
        **public_index,
        "component_bundles": {
            name: {
                "bundle_sha256": track_details[name]["bundle_sha256"],
                "participant_file_count": track_details[name][
                    "participant_file_count"
                ],
                "public_directory": "public/{0}".format(name),
                "private_evidence_directory": "private/{0}".format(name),
            }
            for name in TRACKS
        },
    }
    json_dump(destination / "pack.json", pack_index)
    json_dump(
        destination / "COMPLETE.json",
        {
            "schema_version": PACK_SCHEMA_VERSION,
            "pack_id": spec.pack_id,
            "status": "complete",
            "validated_at": utc_now(),
            "pack_sha256": _tree_digest(destination),
        },
    )
    validate_combined_pack(destination)
    return destination


def validate_combined_pack(pack_dir: str | os.PathLike[str]) -> dict[str, Any]:
    """Validate a combined pack and its separated component evidence."""

    root = Path(pack_dir).resolve()
    required = [
        root / "pack.json",
        root / "COMPLETE.json",
        root / "public" / "pack.json",
        root / "public" / "instructions.md",
        root / "public" / "checksums.sha256",
    ]
    for track_name in TRACKS:
        required.extend(
            [
                root / "public" / track_name / "sample_sheet.csv",
                root / "public" / track_name / "submission_schema.json",
                root / "private" / track_name / "answer_key.json",
                root / "private" / track_name / "scoring_policy.json",
                root / "private" / track_name / "release.json",
                root / "private" / track_name / "COMPLETE.json",
            ]
        )
    missing = [str(path.relative_to(root)) for path in required if not path.is_file()]
    if missing:
        raise ValueError("combined pack is missing: {0}".format(", ".join(missing)))
    if any(path.name == "private" for path in (root / "public").rglob("*")):
        raise ValueError("private evidence is present beneath the public directory")
    index = _load_json(root / "pack.json")
    public_index = _load_json(root / "public" / "pack.json")
    complete = _load_json(root / "COMPLETE.json")
    if index.get("schema_version") != PACK_SCHEMA_VERSION:
        raise ValueError("unsupported combined pack schema version")
    if public_index.get("pack_id") != index.get("pack_id"):
        raise ValueError("public and private pack identifiers differ")
    if complete.get("pack_sha256") != _tree_digest(root):
        raise ValueError("completed combined pack digest does not match pack")
    checksum_rows = {}
    with open(root / "public" / "checksums.sha256", encoding="utf-8") as handle:
        for row in handle:
            digest, relative = row.rstrip("\n").split("  ", 1)
            checksum_rows[relative] = digest
    public_files = {
        path.relative_to(root / "public").as_posix(): sha256_file(path)
        for path in (root / "public").rglob("*")
        if path.is_file() and path != root / "public" / "checksums.sha256"
    }
    if checksum_rows != public_files:
        raise ValueError("public checksum inventory differs from combined pack files")
    for track_name, exercise in TRACKS.items():
        with open(
            root / "public" / track_name / "sample_sheet.csv",
            encoding="utf-8",
            newline="",
        ) as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            raise ValueError("{0} result template has no samples".format(track_name))
        component_index = _load_json(root / "private" / track_name / "release.json")
        if component_index.get("exercise") != exercise:
            raise ValueError("{0} component exercise is incorrect".format(track_name))
    return {
        "schema_version": PACK_SCHEMA_VERSION,
        "pack_id": index["pack_id"],
        "status": "passed",
        "tracks": list(TRACKS),
        "samples": {
            track["name"]: track["sample_count"]
            for track in public_index.get("tracks", [])
        },
    }
