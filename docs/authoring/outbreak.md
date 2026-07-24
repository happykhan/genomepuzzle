# Phylogeny and outbreak

Outbreak releases provide paired reads and safe epidemiological metadata.
Participants infer relatedness, identify clusters and record samples that
should be excluded or treated cautiously.

## Native simulation

Native generation starts from one frozen reference and private metadata:

```toml
[inputs]
base_genome = "../../sources/outbreak/reference.fasta"
metadata_csv = "outbreak-metadata.csv"
```

GenomePuzzle introduces shared cluster mutations followed by private
sample-specific mutations, then simulates paired reads with ART.

The metadata CSV requires:

```csv
Sample,Cluster,SPECIES,ward,collection_date
tip-01,A,Klebsiella pneumoniae,ICU,2026-01-05
tip-02,A,Klebsiella pneumoniae,ICU,2026-01-07
```

`Sample`, `Cluster` and `SPECIES` remain private truth fields. Additional
columns become participant context after private source identifiers are
replaced with anonymous sample IDs.

## Pre-simulated cohorts

To package an existing cohort, use `source_dir` instead of `base_genome` and
retain `metadata_csv`. This route is useful when reads were produced by a
validated external simulation process.

## Interpretation

Cluster truth and QC truth are deliberately separate. A failed sample can have
a known biological cluster while being excluded from partition scoring.

Cluster labels are scored as a partition. Participants do not need to reproduce
the organiser's arbitrary label names.

## Supported implants

| Implant | Important parameters | Intended signal |
| --- | --- | --- |
| `LOW_COVERAGE` | `read_fraction` | An isolate unsuitable for confident placement |
| `TEN_READ_PAIRS` | — | Exactly ten paired reads |
| `ZERO_BYTE_R1`, `ZERO_BYTE_R2` | — | A required mate is a literal zero-byte file |
| `MISSING_R1`, `MISSING_R2` | — | A required mate is absent |
| `CONTAMINATED` | `contaminant_source_id`, `contamination_fraction` | A 30–90% different-species mixture |
| `WRONG_ORGANISM` | `replacement_source_id` | Complete replacement with another organism |

`NORMAL` and `NONE` leave the simulated or pre-simulated pair intact. A
faulted sample keeps its biological cluster truth, while its expected
`qc_status` becomes `FAIL` and it is omitted from partition scoring.
