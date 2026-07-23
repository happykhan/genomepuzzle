# 2026 pilot releases

`typing-practice.toml` freezes the first practice typing cohort. It includes
six unmodified assemblies and two deliberately troublesome assemblies:

- `Sample_TP007` is split into 1 kb fragments.
- `Sample_TP008` contains contigs from a second assembly at approximately 12%
  of the target assembly length.

Source accessions and implants belong only in the generator repository and
the private release package. Participant files contain anonymous sample and
contig identifiers.

Build the release from the repository root:

```bash
pixi run genomepuzzle release validate-spec \
  --spec releases/2026-pilot/typing-practice.toml

pixi run genomepuzzle release fetch-assemblies \
  --spec releases/2026-pilot/typing-practice.toml \
  --output-dir generated/2026-pilot/typing-sources

pixi run genomepuzzle release build-typing \
  --spec releases/2026-pilot/typing-practice.toml \
  --source-dir generated/2026-pilot/typing-sources \
  --output-dir generated/2026-pilot/typing-practice
```

Do not publish a build made with `--skip-analysis`: its private answer key is
explicitly marked `pending_kleborate`. Review the final Kleborate calls for
the two implanted samples before publishing. Kleborate is pinned in
`pixi.toml`; run all release commands via `pixi run` to use that version.
