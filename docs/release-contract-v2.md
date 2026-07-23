# GenomePuzzle release contract v2

The release directory is the only interface between GenomePuzzle and
GHRUPuzzles:

```text
release/
├── release.json
├── public/
│   ├── manifest.json
│   ├── dataset_manifest.json
│   ├── sample_sheet.csv
│   ├── submission_schema.json
│   ├── instructions.md
│   ├── checksums.sha256
│   └── files/
├── private/
│   ├── answer_key.json
│   ├── scoring_policy.json
│   ├── implant_manifest.json
│   ├── provenance.json
│   └── validation_report.json
├── build/
│   ├── plan.json
│   ├── stages/
│   ├── scripts/
│   └── logs/
└── COMPLETE.json
```

`release.json` is the canonical import entry point. `manifest.json` contains
only participant-safe metadata and file roles. `dataset_manifest.json` is an
identical compatibility copy for the current GHRUPuzzles publisher and may be
removed in a later schema version.

The JSON Schemas are packaged under `genomepuzzle/schemas/`:

- `release.schema.json`;
- `manifest.schema.json`; and
- `submission.schema.json`.

GHRUPuzzles owns storage and server-side implementations of safe scorer types.
GenomePuzzle owns dataset truth, submission fields, scorer configuration,
instructions, implants, provenance and biological validation.

The website must not duplicate exercise columns, infer answer fields from
whatever keys happen to occur, or perform biological generation.

## File roles

Per-sample file roles are stable:

- `assembly`;
- `read_1`;
- `read_2`; and
- `long_reads`.

The website maps these roles to authenticated download links after import.
Storage URLs are not embedded in the generated bundle.

## Completion

`COMPLETE.json` contains the release ID, validation timestamp, status and
content digest. GHRUPuzzles must reject an upload if the digest or any
per-file checksum differs.

`build/` contains mutable operational state and is not part of the published
content digest. It may be retained for audit or excluded from R2 publication.
