# Building assessment releases

`genomepuzzle release` is the supported interface for all new datasets.
Historical generators are quarantined under `genomepuzzle legacy` and do not
produce publishable v2 releases.

## Managed software

Pixi is the only software provider. The Linux cluster environment and every
biological tool are pinned in `pixi.toml` and `pixi.lock`. There is no
repository `bin/` directory and there are no Docker or PATH fallback images.

Install and test the exact environment:

```bash
pixi install --locked
pixi run test
```

## Cluster workflow

Heavy generation and analysis always run through SLURM:

```bash
export GENOMEPUZZLE_ID_SALT="$(openssl rand -hex 32)"

pixi run genomepuzzle release plan \
  --spec releases/my-round/typing-practice.toml \
  --output-dir generated/typing-practice

pixi run genomepuzzle release build \
  --spec releases/my-round/typing-practice.toml \
  --output-dir generated/typing-practice
```

`build` writes a persisted plan and submits a `generate` job followed by an
independent `validate` job with an `afterok` dependency. It never runs a heavy
stage locally when submission fails.

Inspect a running or failed workflow:

```bash
pixi run genomepuzzle release status \
  --plan generated/typing-practice/build/plan.json

pixi run genomepuzzle release logs \
  --plan generated/typing-practice/build/plan.json

pixi run genomepuzzle release resume \
  --plan generated/typing-practice/build/plan.json
```

Each stage records its command, resources, attempts, job ID, node, timestamps,
status and failure. A retry removes only partial release artifacts; it
preserves the build plan, scripts, attempts and logs.

## Release specification

All identities, inputs, implants and expected interpretation live in one
private TOML specification:

```toml
schema_version = "1.0"
release_id = "2026-round-1-assembly-practice"
exercise = "assembly"
mode = "practice"
master_seed = 470
title = "Short-read assembly practice"
description = "Assemble paired reads and identify problematic datasets."
pass_threshold = 0.8
instructions = [
  "Assemble every paired-read dataset.",
  "Return the completed sample_sheet.csv."
]

[inputs]
source_dir = "../../sources/assemblies"

[[samples]]
source_id = "GCA_000000001.1"
identity_key = "assembly-clean-1"
[samples.expected_answers]
species = "Klebsiella pneumoniae"

[[samples]]
source_id = "GCA_000000002.1"
identity_key = "assembly-low-coverage-1"
implant = "LOW_COVERAGE"
[samples.implant_parameters]
read_fraction = 0.15
[samples.expected_answers]
species = "Klebsiella pneumoniae"
```

Public IDs and tool seeds are derived from the release ID, stable identity key,
master seed and private salt. Reordering samples does not change them.
Practice specifications may freeze a `public_id`; challenge specifications
normally omit it.

Paths under `[inputs]` are resolved relative to the specification:

- typing: `source_dir` containing `<source_id>.fasta`;
- assembly and hybrid: `source_dir` containing reference assemblies;
- native outbreak simulation: `base_genome` and `metadata_csv`; or
- pre-simulated outbreak input: `source_dir` and `metadata_csv`.

## Exercise generation

### Genotyping

Typing accepts `NORMAL`, `FRAGMENTED` and `MIXED_CONTIGS`. Kleborate 3.1.3 is
run on the final anonymous FASTA, not the source assembly. Species and `st` are
normalised to the website contract.

```toml
[inputs]
source_dir = "../../sources/typing"

[[samples]]
source_id = "GCA_000000099.1"
implant = "MIXED_CONTIGS"
[samples.implant_parameters]
contaminant_source_id = "GCA_000000100.1"
contamination_fraction = 0.10
```

### Short-read assembly

ART creates paired reads directly from the frozen reference assembly.
Supported implants are `LOW_COVERAGE`, `POOR_QUALITY`, `TRUNCATED` and
`CONTAMINATED`. Every troublesome sample must materialise its requested
implant; silent fallback to `NORMAL` is forbidden.

### Hybrid assembly

ART and Badread produce the short- and long-read tracks. Supported implants
are `LOW_SHORT_COVERAGE`, `LOW_LONG_COVERAGE`, `LONG_READ_QUALITY` and
`CONTAMINATED`. A contaminated hybrid sample receives contamination in both
data modalities.

### Phylogeny and outbreak

The native generator creates shared cluster mutations and private
sample-specific mutations from a frozen base genome, then simulates paired
reads with ART. The metadata CSV requires `Sample`, `Cluster` and `SPECIES`.
Public metadata may include epidemiological fields, but cluster and species
truth remain private.

```toml
[inputs]
base_genome = "../../sources/outbreak/reference.fasta"
metadata_csv = "outbreak_metadata.csv"
```

Existing pre-simulated read cohorts remain supported by using `source_dir`
instead of `base_genome`.

## Acceptance

A release is publishable only when `COMPLETE.json` exists. Sealing fails when:

- a participant file is absent, malformed or has the wrong checksum;
- sample membership differs between manifest, template and answer key;
- an answer is pending or lies outside the submission schema;
- a scored field is missing;
- a source identity remains in a sequence header;
- an implant lacks a passing validation record; or
- any required public/private contract artifact is absent.

Validate and inspect without biological recomputation:

```bash
pixi run genomepuzzle release validate \
  --release-dir generated/typing-practice \
  --require-complete

pixi run genomepuzzle release inspect \
  --release-dir generated/typing-practice
```

The complete digest covers every published public/private artifact. Mutable
SLURM state under `build/` is intentionally excluded.
