# Short-read assembly

Assembly releases simulate paired Illumina reads from frozen reference
assemblies using ART. Participants assemble the reads, identify the organism
and report a QC interpretation.

## Supported implants

| Implant | Important parameters | Intended signal |
| --- | --- | --- |
| `LOW_COVERAGE` | `read_fraction` | No more than approximately 1× depth |
| `TEN_READ_PAIRS` | — | Exactly ten paired reads |
| `ZERO_BYTE_R1`, `ZERO_BYTE_R2` | — | A required mate is a literal zero-byte file |
| `MISSING_R1`, `MISSING_R2` | — | A required mate is absent |
| `CONTAMINATED` | `contaminant_source_id`, `contamination_fraction` | A 30–90% mixture with a different species |
| `WRONG_ORGANISM` | `replacement_source_id` | Complete replacement with another organism |

Short-read depth defaults to 30×. Simulation parameters can be set per sample:

```toml
[samples.implant_parameters]
short_coverage = 40
read_length = 150
fragment_length = 300
fragment_sd = 50
```

Additional contaminant or replacement references must be present in the same
`source_dir`. Contamination defaults to 50%; GenomePuzzle records its achieved
fraction from the final reads.

## Outputs

Each public sample receives `read_1` and `read_2` FASTQ roles. The private
release records source checksums, simulation seed, requested and achieved
implant parameters, expected answers and post-implant validation.

Normal samples retain strict paired-FASTQ validation. Missing and zero-byte
roles are accepted only when the exact private `fault_type` declares them.
