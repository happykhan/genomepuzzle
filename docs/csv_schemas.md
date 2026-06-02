# CSV Schemas

This document defines the expected CSV inputs and generated CSV outputs used by `genomepuzzle`.

## `datasets/samplelist.csv`

Used by `short simulate`, `short contamination`, and legacy workflows.

Required columns:
- `SAMPLE_NAME`
- `SHORT_READS`
- `ASSEMBLY`
- `SPECIES`
- `ST`
- `AMR`
- `USE_ORIGINAL_READS`

Semantics:
- `SAMPLE_NAME`: internal source identifier
- `SHORT_READS`: SRA run accession when original reads should be downloaded
- `ASSEMBLY`: NCBI assembly accession
- `SPECIES`: organism label used for filtering and reporting
- `ST`: sequence type label
- `AMR`: antimicrobial resistance label or summary
- `USE_ORIGINAL_READS`: `TRUE` or `FALSE`

## Clean short-read `sample_sheet.csv`

Produced by `short simulate`. Used as input to `short errors`.

Typical columns:
- `SAMPLE_NAME`
- `SHORT_READS`
- `ASSEMBLY`
- `SPECIES`
- `ST`
- `AMR`
- `USE_ORIGINAL_READS`
- `public_name`
- `r1`
- `r2`
- `coverage`
- `read_length`
- `platform`
- `fragment_length`
- `standard_deviation`
- `random_seed`
- `QC`
- `ERROR`
- `Notes`

## Short-read implanted outputs

Produced by `short errors`.

### Public `sample_sheet.csv`

Columns:
- `ID`
- `R1`
- `R2`
- `SPECIES`
- `QC`
- `ERROR`
- `ST`
- `AMR`
- `Notes`

Public values intentionally hide the real answer key.

### Private `answer_sheet.csv`

Columns:
- `ID`
- `R1`
- `R2`
- `SPECIES`
- `QC`
- `ERROR`
- `ST`
- `AMR`
- `Notes`

### Private `implant_manifest.csv`

Columns:
- `sample_name`
- `species`
- `source_r1`
- `source_r2`
- `error_type`
- `severity`
- `notes`
- `read_count`

## `datasets/rapid_data.csv`

Used by `rapid` and `long hybrid`.

Required columns:
- `accession`

Common optional columns:
- `organism_organismname`
- `species`
- `assemblystats_totalsequencelength`
- `short_read_coverage`
- `public_name`

## Hybrid public `sample_sheet.csv`

Produced by `long hybrid`.

Columns:
- `sample_name`
- `reference_accession`
- `species`
- `tax_classification`
- `r1`
- `r2`
- `long_reads`
- `assembler`
- `qc`
- `notes`

## Hybrid private `answer_sheet.csv`

Columns:
- `public_name`
- `species`
- `reference_accession`
- `tax_classification`
- `assembler`
- `qc`
- `notes`
- `error_type`
- `severity`

## Hybrid private `implant_manifest.csv`

Columns:
- `sample_name`
- `reference_accession`
- `species`
- `error_type`
- `severity`
- `notes`
- `short_read_count`
- `long_read_count`
