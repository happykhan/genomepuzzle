# Building assessment releases

The supported release path is now the `genomepuzzle release` command group.
The old `eqa-test.py` script remains only as historical reference.

## Execution policy

Run computationally heavy generation and biological analysis through SLURM,
including read simulation, assemblies, batch Kleborate analysis, phylogenetic
inference, outbreak generation, full-dataset implants, and QC. Local commands
are limited to specification validation, source indexing, packaging,
manifests, checksums, tests, and small smoke checks.

Every heavy job must record its SLURM job ID, requested resources, tool
versions, seeds, and output directory in private provenance. A failed SLURM
submission must not silently fall back to local execution.

## 1. Write and validate a specification

```toml
schema_version = "1.0"
release_id = "2026-round-1-typing-practice"
exercise = "typing"
mode = "practice"
master_seed = 42

[[samples]]
source_id = "GCF_000000001.1"
identity_key = "typing-clean-1"

[[samples]]
source_id = "GCF_000000002.1"
identity_key = "typing-fragmented-1"
implant = "FRAGMENTED"
[samples.implant_parameters]
fragment_size = 1000
```

Resolve IDs without generating data:

```bash
export GENOMEPUZZLE_ID_SALT="$(openssl rand -hex 32)"
pixi run genomepuzzle release validate-spec \
  --spec releases/2026-round-1/typing-practice.toml \
  --output-json generated/private-typing-map.json
```

Keep the salt, resolved mapping, specifications, and every `private/` output
outside participant-accessible storage.

## 2. Fetch and build genotyping inputs

For NCBI assembly accessions, fetch and checksum every target and contaminant
declared by the specification:

```bash
pixi run genomepuzzle release fetch-assemblies \
  --spec releases/2026-round-1/typing-practice.toml \
  --output-dir generated/typing-sources
```

The cache contains `<source_id>.fasta` plus `sources.json`; a second run reuses
it. Use `--refresh` only when deliberately recalibrating against a newer NCBI
assembly version. Locally supplied assemblies can instead be staged under the
same filename convention.

Supported initial implants are `NORMAL`, `FRAGMENTED`, and `MIXED_CONTIGS`. A
mixed sample must name its private contaminant:

```toml
implant = "MIXED_CONTIGS"
[samples.implant_parameters]
contaminant_source_id = "GCF_000000099.1"
contamination_fraction = 0.10
```

For a single smoke-test sample, the following command can check the packaging
path. Run the complete Kleborate cohort as a generated SLURM job:

```bash
pixi run genomepuzzle release build-typing \
  --spec releases/2026-round-1/typing-practice.toml \
  --source-dir generated/typing-sources \
  --output-dir generated/typing-practice
```

Kleborate 3.1.3 is pinned in `pixi.toml`, so production commands must run
through `pixi run`. The executable override is intended for a separately
pinned wrapper or container, not an unversioned system installation.

`--skip-analysis` is only for preparation and testing. It creates answers
marked `pending_kleborate` and must not be used for a published assessment.

## 3. Package short-read and hybrid releases

The existing simulators still perform the expensive biological generation.
Once implants and reference analysis are complete, prepare a private JSON
object keyed by source ID:

```json
{
  "source-a": {
    "species": "Klebsiella pneumoniae",
    "qc": "failed",
    "diagnosis": "low coverage"
  }
}
```

Final source files use `<source_id>_R1.fastq.gz`,
`<source_id>_R2.fastq.gz`, and, for hybrid data,
`<source_id>_long.fastq.gz`.

```bash
pixi run genomepuzzle release package-reads \
  --spec releases/2026-round-1/hybrid-practice.toml \
  --source-dir final-implanted-reads \
  --expected-answers private-expected-answers.json \
  --output-dir generated/hybrid-practice
```

This step replaces source names with canonical public IDs, checks every
expected file, and creates paired public/private manifests.

## 4. Build outbreak releases

TreeToReads simulation output is frozen before packaging. The source directory
contains paired FASTQs named for tree tips, and the metadata CSV contains a
`Sample` column plus private `Cluster` and `SPECIES` truth.

Supported initial implants are `NORMAL`, `LOW_COVERAGE`, and `CONTAMINATED`.
Cluster truth remains separate from the expected QC inclusion decision.

```bash
pixi run genomepuzzle release build-outbreak \
  --spec releases/2026-round-1/outbreak-practice.toml \
  --source-dir treetoreads-output \
  --metadata outbreak_trees/rooted_outbreak_4_clusters.csv \
  --output-dir generated/outbreak-practice
```

## 5. Release acceptance

A release is ready for the portal only when:

- `COMPLETE` exists;
- `public/dataset_manifest.json` contains no source IDs, implants, or answers;
- every participant file is listed with its size and SHA-256 digest;
- `private/answer_key.json`, `provenance.json`, and
  `implant_manifest.json` are complete;
- expected analyses were run on the final participant files; and
- the intended troublesome samples were manually reviewed.

Publication is handled by `ghrupuzzle/scripts/publish_release.py`. It validates
the package again and uploads public inputs and private truth with separate
credentials.
