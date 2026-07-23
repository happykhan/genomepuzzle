# Short-read assembly

Assembly releases simulate paired Illumina reads from frozen reference
assemblies using ART. Participants assemble the reads, identify the organism
and report a QC interpretation.

## Supported implants

| Implant | Important parameters | Intended signal |
| --- | --- | --- |
| `LOW_COVERAGE` | `read_fraction` | Insufficient depth after simulation |
| `POOR_QUALITY` | `min_quality`, `max_quality` | Uniformly degraded base qualities |
| `TRUNCATED` | `read_length` | Unusually short reads |
| `CONTAMINATED` | `contaminant_source_id`, `contamination_fraction` | Cross-sample or cross-species mixture |

Short-read depth defaults to 30×. Simulation parameters can be set per sample:

```toml
[samples.implant_parameters]
short_coverage = 40
read_length = 150
fragment_length = 300
fragment_sd = 50
```

For a contaminated sample, the contaminant reference must be present in the
same `source_dir`.

## Outputs

Each public sample receives `read_1` and `read_2` FASTQ roles. The private
release records source checksums, simulation seed, requested and achieved
implant parameters, expected answers and post-implant validation.

Do not use malformed gzip files as the primary troublesome case. The implant
should exercise assembly or QC judgement, not merely test whether a pipeline
can open a file.
