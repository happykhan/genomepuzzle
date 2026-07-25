# Combined assembly EQA release pack

The 14 August 2026 participant pack presents the Challenge 2 short-read
assembly and hybrid assembly releases as one EQA exercise. The two contracts
remain independently scoreable and auditable.

The pack definition is
`releases/2026-website/combined-eqa-challenge-2.toml`. It pins the exact
component release IDs and sealed bundle digests already recorded in the
Challenge 2 publication receipt.

## Build

Use the validated contract 2.1 release directories. This is lightweight
packaging and may run locally:

```bash
pixi run genomepuzzle release build-combined-pack \
  --spec releases/2026-website/combined-eqa-challenge-2.toml \
  --short-read-release generated/2026-website/challenge-2-assembly \
  --long-read-release generated/2026-website/challenge-2-hybrid \
  --output-dir generated/2026-website/challenge-2-combined-assembly-eqa
```

The command refuses an unsealed component, a release ID or exercise mismatch,
and any component digest other than the one pinned in the pack definition.
Validate a retained or transferred pack with:

```bash
pixi run genomepuzzle release validate-combined-pack \
  --pack-dir generated/2026-website/challenge-2-combined-assembly-eqa
```

## Pack layout

```text
challenge-2-combined-assembly-eqa/
├── pack.json
├── COMPLETE.json
├── public/
│   ├── instructions.md
│   ├── pack.json
│   ├── checksums.sha256
│   ├── short-read/
│   │   ├── sample_sheet.csv
│   │   ├── submission_schema.json
│   │   └── files/
│   └── long-read/
│       ├── sample_sheet.csv
│       ├── submission_schema.json
│       └── files/
└── private/
    ├── short-read/
    │   ├── answer_key.json
    │   ├── scoring_policy.json
    │   └── provenance.json
    └── long-read/
        ├── answer_key.json
        ├── scoring_policy.json
        └── provenance.json
```

Only `public/` is participant-facing. The `private/short-read` and
`private/long-read` directories retain separate answer, scoring, fault,
validation, and provenance evidence. Never send or publish the top-level pack
as though every path were public.

## Expected participant outputs

Participants return two UTF-8 CSV files for one EQA exercise:

1. `short-read/sample_sheet.csv`, completed for the short-read samples.
2. `long-read/sample_sheet.csv`, completed for the hybrid samples.

Both templates use:

```csv
sample_id,qc_status,failure_reason,species,assembler,contig_count,total_length,n50,longest_contig,notes
```

`sample_id` identifies the row. `qc_status`, `failure_reason`, and species for
passing samples are the automatically scored evidence. Assembly workflow and
statistics support review but are not exact-value scored because valid
workflows differ.

The sheets must remain separate. Combining the rows is ambiguous because
sample identifiers are release-specific and the allowed failure vocabulary
differs between tracks.

## Participant wording decisions

The generated `public/instructions.md` makes these distinctions explicit:

- “long-read” means the hybrid track, not a long-read-only assembly;
- both short- and long-read evidence must be assessed for a passing hybrid
  sample;
- a missing role or zero-byte file may be an intentional assessment input;
- participants must not repair, invent, or substitute missing evidence;
- `NONE` is valid only with `PASS`;
- sample order or similar-looking names do not imply correspondence between
  tracks; and
- two completed result sheets are separate evidence for one participant
  exercise.

## Release-date preflight

The pack records 14 August 2026 as its release date. The GHRU Puzzles portal
controls download and submission access independently. Before circulating
participant copy, confirm the portal round opens on the intended date and use
one of these unambiguous phrases:

- “The exercise pack will be released on 14 August” when downloads open that
  day; or
- “Participant instructions will be issued on 14 August; the timed exercise
  opens at [time and date]” when the portal opens later.

Do not say “opens on 14 August” unless the configured portal timestamp agrees.
At the time this pack definition was prepared, the website repository's
Challenge 2 migration used `2026-08-16T23:00:00Z`; release coordination must
either change that configuration or state the later opening explicitly.

## Final checks

Before handoff:

1. build against the exact sealed component bundles;
2. run `validate-combined-pack`;
3. verify `public/checksums.sha256` after transfer;
4. open both templates and confirm they contain ten track-specific sample
   rows;
5. confirm the public tree contains no `private/` directory;
6. participant-test the missing long-read role and zero-byte short-read mate;
7. submit each CSV to its matching component release; and
8. reconcile the participant wording with the portal opening timestamp.
