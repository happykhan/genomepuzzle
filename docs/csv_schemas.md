# Participant CSV contract

Every v2 release contains:

- `public/submission_schema.json`, the machine-readable source of truth; and
- `public/sample_sheet.csv`, a blank participant template containing every
  public sample identifier.

GHRUPuzzles must render instructions and validate uploads from this schema
rather than maintaining a second list of columns.

## Shared rules

- The first column is always `sample_id`.
- Participants must not add, remove or rename sample identifiers.
- Column names are lowercase snake case.
- Required, scored, unscored and alias behaviour is explicit in
  `submission_schema.json` and `private/scoring_policy.json`.
- Optional fields remain present in the template.

## Genotyping

```text
sample_id,species,st,k_locus,capsule_type,wzi,o_locus,o_type,bla_carb
```

`st` is the canonical sequence-type name. Analyser-specific names such as
`kleborate_st` never appear in participant or website artifacts.
`bla_carb` is scored as an order-independent list.

## Short-read assembly

```text
sample_id,species,qc,error,notes
```

## Hybrid assembly

```text
sample_id,species,assembler,qc,error,notes
```

Source/reference accessions are deliberately excluded.

## Phylogeny and outbreak

```text
sample_id,cluster,species,qc_decision,notes
```

Cluster labels are arbitrary. Scoring compares the inferred partition rather
than requiring participants to reproduce organiser label names.

## Epidemiological metadata

The private outbreak-generation metadata CSV requires:

- `Sample`: private source/tip identifier;
- `Cluster`: private expected cluster;
- `SPECIES`: private expected species.

Any other columns are included as public contextual metadata unless explicitly
classified as private by the generator. Private source identifiers are always
replaced by the release sample IDs.

Legacy CSV formats are accepted only by commands under `genomepuzzle legacy`;
they are not valid v2 release bundles.
