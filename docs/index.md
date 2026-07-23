# GenomePuzzle

GenomePuzzle builds reproducible microbial genomics datasets for training and
proficiency assessment. It creates anonymous participant files alongside a
strictly separated private answer key, provenance record, implant manifest and
validation report.

Four exercise types use the same release workflow:

| Exercise | Participant material | Typical task |
| --- | --- | --- |
| Genotyping | Anonymous assemblies | Run Kleborate and interpret typing calls |
| Short-read assembly | Paired Illumina reads | Assemble, identify and assess QC |
| Hybrid assembly | Paired short reads and long reads | Build and evaluate a hybrid assembly |
| Outbreak | Paired reads and public metadata | Reconstruct relatedness and identify clusters |

Every dataset can include realistic troublesome samples. The fault is recorded
privately and must pass an implant-specific validation before the bundle can be
sealed.

## What GenomePuzzle owns

GenomePuzzle is the dataset-generation side of the assessment platform. It
owns:

- deterministic sample identities and simulation seeds;
- biological generation and implanted problems;
- expected analyses and scoring configuration;
- public/private separation;
- checksums, provenance and validation; and
- a versioned bundle consumed by a delivery website.

Authentication, participant submissions, manual review and certificates belong
in the delivery application, not in this repository.

## A release at a glance

A versioned TOML specification becomes a resolved SLURM plan. The cluster
generates and analyses the dataset, then an independent validation stage checks
the public and private artifacts. Passing releases receive a `COMPLETE.json`
seal and can be handed to the delivery website.

Heavy biological work is always submitted to SLURM. Pixi supplies the exact
software environment, so there are no bundled binaries, Docker fallbacks or
silent use of whatever happens to be on `PATH`.

[Install GenomePuzzle](getting-started/installation.md){ .md-button .md-button--primary }
[Build a first release](getting-started/first-release.md){ .md-button }

!!! warning "Private release contents"

    A complete release contains source identities, expected answers and
    implanted-error details under `private/`. Never place the complete release
    directory in a public bucket. Publish only through a consumer that
    understands the [release contract](release-contract-v2.md).
