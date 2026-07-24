"""Versioned release specifications, identities, and manifest packaging."""

from __future__ import annotations

import hashlib
import hmac
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from genomepuzzle.contract import (
    BUNDLE_SCHEMA_VERSION,
    exercise_contract,
    failure_reason_for_implant,
    json_dump,
    normalize_answer_fields,
    write_participant_contract,
    write_release_index,
)
from genomepuzzle.provenance import runtime_provenance

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10
    import tomli as tomllib


RELEASE_SCHEMA_VERSION = "1.0"
EXERCISES = {"typing", "assembly", "hybrid", "outbreak"}
MODES = {"practice", "challenge"}
NORMAL_IMPLANTS = {"NORMAL", "NONE"}
SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
FORBIDDEN_PUBLIC_KEYS = {
    "answer",
    "answers",
    "answer_key",
    "expected",
    "expected_answers",
    "fault",
    "fault_type",
    "implant",
    "implant_type",
    "reference_accession",
    "source",
    "source_id",
    "source_name",
    "truth",
}


@dataclass(frozen=True)
class ReleaseSampleSpec:
    """One resolved source sample requested by a release specification."""

    source_id: str
    identity_key: str
    implant: str = "NORMAL"
    implant_parameters: Mapping[str, Any] = field(default_factory=dict)
    expected_answers: Mapping[str, Any] = field(default_factory=dict)
    public_id: str | None = None


@dataclass(frozen=True)
class ReleaseSpec:
    """Minimal shared configuration required by every exercise generator."""

    release_id: str
    exercise: str
    mode: str
    master_seed: int
    samples: tuple[ReleaseSampleSpec, ...]
    schema_version: str = RELEASE_SCHEMA_VERSION
    id_salt_env: str = "GENOMEPUZZLE_ID_SALT"
    title: str | None = None
    description: str | None = None
    instructions: tuple[str, ...] = ()
    pass_threshold: float = 0.8
    inputs: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class ResolvedReleaseSample:
    """A source sample after its anonymous identity and seed are resolved."""

    source_id: str
    identity_key: str
    sample_id: str
    random_seed: int
    implant: str
    implant_parameters: Mapping[str, Any] = field(default_factory=dict)
    expected_answers: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ReleaseArtifactSample:
    """Final participant files plus the private truth used to produce them."""

    sample_id: str
    source_id: str
    random_seed: int
    files: Mapping[str, str]
    expected_answers: Mapping[str, Any]
    implant: str = "NORMAL"
    implant_parameters: Mapping[str, Any] = field(default_factory=dict)
    public_metadata: Mapping[str, Any] = field(default_factory=dict)
    private_provenance: Mapping[str, Any] = field(default_factory=dict)


def _require_identifier(value: str, field_name: str) -> str:
    if not value or not SAFE_IDENTIFIER.fullmatch(value):
        raise ValueError(
            "{field} must contain only letters, numbers, '.', '_' or '-'".format(
                field=field_name
            )
        )
    return value


def _require_non_empty(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("{field} must be a non-empty string".format(field=field_name))
    return value.strip()


def load_release_spec(path: str | os.PathLike[str]) -> ReleaseSpec:
    """Load and validate a versioned TOML release specification."""

    with open(path, "rb") as handle:
        raw = tomllib.load(handle)

    schema_version = str(raw.get("schema_version", RELEASE_SCHEMA_VERSION))
    if schema_version != RELEASE_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported release schema version: {0}".format(schema_version)
        )

    release_id = _require_identifier(
        _require_non_empty(raw.get("release_id"), "release_id"), "release_id"
    )
    exercise = _require_non_empty(raw.get("exercise"), "exercise").lower()
    if exercise not in EXERCISES:
        raise ValueError(
            "exercise must be one of: {0}".format(", ".join(sorted(EXERCISES)))
        )
    mode = _require_non_empty(raw.get("mode"), "mode").lower()
    if mode not in MODES:
        raise ValueError("mode must be one of: {0}".format(", ".join(sorted(MODES))))

    master_seed = raw.get("master_seed", 42)
    if not isinstance(master_seed, int) or master_seed < 0:
        raise ValueError("master_seed must be a non-negative integer")

    raw_samples = raw.get("samples")
    if not isinstance(raw_samples, list) or not raw_samples:
        raise ValueError("release specification must contain at least one [[samples]]")

    samples = []
    identity_keys = set()
    public_ids = set()
    for index, raw_sample in enumerate(raw_samples, start=1):
        if not isinstance(raw_sample, dict):
            raise ValueError("samples entry {0} must be a table".format(index))
        source_id = _require_non_empty(
            raw_sample.get("source_id"), "samples[{0}].source_id".format(index)
        )
        identity_key = _require_non_empty(
            raw_sample.get("identity_key", source_id),
            "samples[{0}].identity_key".format(index),
        )
        if identity_key in identity_keys:
            raise ValueError("duplicate sample identity_key: {0}".format(identity_key))
        identity_keys.add(identity_key)

        implant = _require_non_empty(
            raw_sample.get("implant", "NORMAL"),
            "samples[{0}].implant".format(index),
        ).upper()
        parameters = raw_sample.get("implant_parameters", {})
        if not isinstance(parameters, dict):
            raise ValueError(
                "samples[{0}].implant_parameters must be a table".format(index)
            )
        expected_answers = raw_sample.get("expected_answers", {})
        if not isinstance(expected_answers, dict):
            raise ValueError(
                "samples[{0}].expected_answers must be a table".format(index)
            )

        public_id = raw_sample.get("public_id")
        if public_id is not None:
            public_id = _require_identifier(
                _require_non_empty(
                    public_id, "samples[{0}].public_id".format(index)
                ),
                "samples[{0}].public_id".format(index),
            )
            if public_id in public_ids:
                raise ValueError("duplicate public_id: {0}".format(public_id))
            public_ids.add(public_id)

        samples.append(
            ReleaseSampleSpec(
                source_id=source_id,
                identity_key=identity_key,
                implant=implant,
                implant_parameters=parameters,
                expected_answers=expected_answers,
                public_id=public_id,
            )
        )

    id_salt_env = _require_identifier(
        _require_non_empty(
            raw.get("id_salt_env", "GENOMEPUZZLE_ID_SALT"), "id_salt_env"
        ),
        "id_salt_env",
    )
    title = raw.get("title")
    if title is not None:
        title = _require_non_empty(title, "title")
    description = raw.get("description")
    if description is not None:
        description = _require_non_empty(description, "description")
    instructions_raw = raw.get("instructions", [])
    if not isinstance(instructions_raw, list) or any(
        not isinstance(item, str) or not item.strip() for item in instructions_raw
    ):
        raise ValueError("instructions must be an array of non-empty strings")
    pass_threshold = raw.get("pass_threshold", 0.8)
    if (
        not isinstance(pass_threshold, (int, float))
        or isinstance(pass_threshold, bool)
        or not 0 < float(pass_threshold) <= 1
    ):
        raise ValueError("pass_threshold must be greater than 0 and at most 1")
    raw_inputs = raw.get("inputs", {})
    if not isinstance(raw_inputs, dict) or any(
        not isinstance(key, str) or not isinstance(value, str)
        for key, value in raw_inputs.items()
    ):
        raise ValueError("inputs must be a table of string paths")
    return ReleaseSpec(
        schema_version=schema_version,
        release_id=release_id,
        exercise=exercise,
        mode=mode,
        master_seed=master_seed,
        id_salt_env=id_salt_env,
        samples=tuple(samples),
        title=title,
        description=description,
        instructions=tuple(item.strip() for item in instructions_raw),
        pass_threshold=float(pass_threshold),
        inputs=dict(raw_inputs),
    )


def derive_sample_id(
    release_id: str,
    identity_key: str,
    id_salt: str,
    prefix: str = "Sample",
) -> str:
    """Create a stable, non-reversible public identifier using a private salt."""

    _require_identifier(release_id, "release_id")
    _require_identifier(prefix, "prefix")
    identity_key = _require_non_empty(identity_key, "identity_key")
    id_salt = _require_non_empty(id_salt, "id_salt")
    digest = hmac.new(
        id_salt.encode("utf-8"),
        "{0}\0{1}".format(release_id, identity_key).encode("utf-8"),
        hashlib.sha256,
    ).hexdigest()[:12]
    return "{0}_{1}".format(prefix, digest)


def derive_sample_seed(release_id: str, identity_key: str, master_seed: int) -> int:
    """Derive a stable positive 31-bit tool seed independent of input order."""

    if not isinstance(master_seed, int) or master_seed < 0:
        raise ValueError("master_seed must be a non-negative integer")
    payload = "{0}\0{1}\0{2}".format(release_id, identity_key, master_seed)
    seed = int.from_bytes(hashlib.sha256(payload.encode("utf-8")).digest()[:4], "big")
    return seed & 0x7FFFFFFF or 1


def resolve_release_samples(
    spec: ReleaseSpec, id_salt: str | None = None
) -> tuple[ResolvedReleaseSample, ...]:
    """Resolve anonymous IDs and stable seeds for all specified samples."""

    needs_salt = any(sample.public_id is None for sample in spec.samples)
    if needs_salt and id_salt is None:
        id_salt = os.environ.get(spec.id_salt_env)
    if needs_salt and not id_salt:
        raise ValueError(
            "Set {0} or pass id_salt to resolve public sample IDs".format(
                spec.id_salt_env
            )
        )

    resolved = []
    seen_ids = set()
    for sample in spec.samples:
        sample_id = sample.public_id or derive_sample_id(
            spec.release_id, sample.identity_key, id_salt or ""
        )
        if sample_id in seen_ids:
            raise ValueError("public sample ID collision: {0}".format(sample_id))
        seen_ids.add(sample_id)
        resolved.append(
            ResolvedReleaseSample(
                source_id=sample.source_id,
                identity_key=sample.identity_key,
                sample_id=sample_id,
                random_seed=derive_sample_seed(
                    spec.release_id, sample.identity_key, spec.master_seed
                ),
                implant=sample.implant,
                implant_parameters=sample.implant_parameters,
                expected_answers=sample.expected_answers,
            )
        )
    return tuple(resolved)


def sha256_file(path: str | os.PathLike[str]) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_available_release_directory(path: str | os.PathLike[str]) -> Path:
    """Allow a workflow's private build state but reject release overwrites."""

    destination = Path(path)
    if destination.exists():
        unexpected = [item.name for item in destination.iterdir() if item.name != "build"]
        if unexpected:
            raise ValueError("release directory is not empty: {0}".format(destination))
    return destination


def _assert_public_metadata_safe(value: Any, path: str = "metadata") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized_key = str(key).lower()
            if normalized_key in FORBIDDEN_PUBLIC_KEYS:
                raise ValueError(
                    "private field {0}.{1} cannot appear in a public manifest".format(
                        path, key
                    )
                )
            _assert_public_metadata_safe(item, "{0}.{1}".format(path, key))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _assert_public_metadata_safe(item, "{0}[{1}]".format(path, index))


def _json_dump(path: Path, payload: Mapping[str, Any]) -> None:
    json_dump(path, payload)


def write_release_manifests(
    release_dir: str | os.PathLike[str],
    spec: ReleaseSpec,
    samples: Sequence[ReleaseArtifactSample],
    generator: Mapping[str, Any] | None = None,
) -> dict[str, str]:
    """Write paired public/private manifests for final participant files.

    Files are not copied. Callers must place them under ``public/files`` before
    calling this function.
    """

    release_path = Path(release_dir)
    public_dir = release_path / "public"
    private_dir = release_path / "private"
    files_dir = public_dir / "files"
    public_dir.mkdir(parents=True, exist_ok=True)
    private_dir.mkdir(parents=True, exist_ok=True)
    files_dir.mkdir(parents=True, exist_ok=True)

    if not samples:
        raise ValueError("at least one artifact sample is required")
    if len(samples) != len(spec.samples):
        raise ValueError(
            "artifact sample count ({actual}) does not match release specification "
            "({expected})".format(actual=len(samples), expected=len(spec.samples))
        )

    public_rows = []
    private_rows = []
    answer_rows = []
    implant_rows = []
    checksum_lines = []
    seen_ids = set()
    normalized_answers_by_sample: list[dict[str, Any]] = []
    for sample in samples:
        _require_identifier(sample.sample_id, "sample_id")
        if sample.sample_id in seen_ids:
            raise ValueError("duplicate artifact sample_id: {0}".format(sample.sample_id))
        seen_ids.add(sample.sample_id)
        _assert_public_metadata_safe(sample.public_metadata)
        if not sample.files:
            raise ValueError("sample {0} has no participant files".format(sample.sample_id))

        file_manifest = {}
        for role, raw_path in sorted(sample.files.items()):
            _require_identifier(str(role), "file role")
            path = Path(raw_path)
            if not path.is_file():
                raise ValueError("participant file does not exist: {0}".format(path))
            filename = path.name
            if not (
                filename == "{0}.fasta".format(sample.sample_id)
                or filename.startswith("{0}_".format(sample.sample_id))
            ):
                raise ValueError(
                    "participant filename must begin with sample_id: {0}".format(
                        filename
                    )
                )
            digest = sha256_file(path)
            size = path.stat().st_size
            file_manifest[str(role)] = {
                "filename": filename,
                "sha256": digest,
                "size": size,
            }
            try:
                relative_path = path.resolve().relative_to(public_dir.resolve())
            except ValueError as exc:
                raise ValueError(
                    "participant file must be inside {0}: {1}".format(public_dir, path)
                ) from exc
            checksum_lines.append(
                "{0}  {1}".format(digest, relative_path.as_posix())
            )

        public_rows.append(
            {
                "sample_id": sample.sample_id,
                "files": file_manifest,
                "metadata": dict(sample.public_metadata),
            }
        )
        private_file_checksums = {
            role: details["sha256"] for role, details in file_manifest.items()
        }
        raw_answers = dict(sample.expected_answers)
        raw_answers.setdefault(
            "qc_status",
            "PASS" if sample.implant in NORMAL_IMPLANTS else "FAIL",
        )
        raw_answers.setdefault(
            "failure_reason",
            failure_reason_for_implant(spec.exercise, sample.implant),
        )
        normalized_answers = normalize_answer_fields(raw_answers)
        expected_qc_status = (
            "PASS" if sample.implant in NORMAL_IMPLANTS else "FAIL"
        )
        expected_failure_reason = failure_reason_for_implant(
            spec.exercise, sample.implant
        )
        if normalized_answers.get("qc_status") != expected_qc_status:
            raise ValueError(
                "sample {0} qc_status does not match private fault type {1}".format(
                    sample.sample_id, sample.implant
                )
            )
        if normalized_answers.get("failure_reason") != expected_failure_reason:
            raise ValueError(
                "sample {0} failure_reason does not match private fault type {1}".format(
                    sample.sample_id, sample.implant
                )
            )
        normalized_answers_by_sample.append(normalized_answers)
        private_rows.append(
            {
                "sample_id": sample.sample_id,
                "source_id": sample.source_id,
                "random_seed": sample.random_seed,
                "fault": {
                    "fault_type": sample.implant,
                    "failure_reason": expected_failure_reason,
                    "parameters": dict(sample.implant_parameters),
                },
                "expected_answers": normalized_answers,
                "provenance": dict(sample.private_provenance),
                "file_checksums": private_file_checksums,
            }
        )
        answer_rows.append(
            {
                "sample_id": sample.sample_id,
                "answers": normalized_answers,
            }
        )
        implant_rows.append(
            {
                "sample_id": sample.sample_id,
                "fault_type": sample.implant,
                "failure_reason": expected_failure_reason,
                "parameters": dict(sample.implant_parameters),
            }
        )

    expected_ids = {
        sample.public_id
        for sample in spec.samples
        if sample.public_id is not None
    }
    if not expected_ids.issubset(seen_ids):
        raise ValueError(
            "artifact sample IDs omit an explicit release specification ID"
        )

    header = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "release_id": spec.release_id,
        "exercise": spec.exercise,
        "mode": spec.mode,
    }
    public_manifest = dict(header)
    public_manifest["samples"] = public_rows
    contract = exercise_contract(spec.exercise)
    public_manifest["title"] = spec.title or contract["title"]
    public_manifest["description"] = spec.description or contract["description"]
    private_manifest = dict(header)
    private_manifest["generator"] = {
        **runtime_provenance(),
        **dict(generator or {}),
    }
    private_manifest["samples"] = private_rows
    answer_key = dict(header)
    answer_key["samples"] = answer_rows
    implant_manifest = dict(header)
    implant_manifest["samples"] = implant_rows
    implant_manifest["summary"] = {
        "normal": sum(
            row["fault_type"] in NORMAL_IMPLANTS for row in implant_rows
        ),
        "troublesome": sum(
            row["fault_type"] not in NORMAL_IMPLANTS for row in implant_rows
        ),
    }

    public_manifest_path = public_dir / "dataset_manifest.json"
    private_manifest_path = private_dir / "provenance.json"
    answer_key_path = private_dir / "answer_key.json"
    implant_manifest_path = private_dir / "implant_manifest.json"
    checksums_path = public_dir / "checksums.sha256"
    _json_dump(public_manifest_path, public_manifest)
    _json_dump(public_dir / "manifest.json", public_manifest)
    _json_dump(private_manifest_path, private_manifest)
    _json_dump(answer_key_path, answer_key)
    _json_dump(implant_manifest_path, implant_manifest)
    with open(checksums_path, "w", encoding="utf-8") as handle:
        for line in sorted(checksum_lines):
            handle.write(line + "\n")

    passing_answers = [
        answers
        for answers in normalized_answers_by_sample
        if answers.get("qc_status") == "PASS"
    ]
    scoring_answers = passing_answers or normalized_answers_by_sample
    common_answer_fields = set(scoring_answers[0])
    for answers in scoring_answers[1:]:
        common_answer_fields &= set(answers)
    common_answer_fields.discard("analysis_status")
    write_participant_contract(
        release_path,
        release_id=spec.release_id,
        exercise=spec.exercise,
        mode=spec.mode,
        sample_ids=[sample.sample_id for sample in samples],
        title=spec.title,
        description=spec.description,
        instructions=spec.instructions or None,
        pass_threshold=spec.pass_threshold,
        scored_fields=common_answer_fields,
    )
    write_release_index(
        release_path,
        release_id=spec.release_id,
        exercise=spec.exercise,
        mode=spec.mode,
        title=spec.title or contract["title"],
        description=spec.description or contract["description"],
    )

    return {
        "public_manifest": str(public_manifest_path),
        "private_manifest": str(private_manifest_path),
        "answer_key": str(answer_key_path),
        "implant_manifest": str(implant_manifest_path),
        "checksums": str(checksums_path),
    }
