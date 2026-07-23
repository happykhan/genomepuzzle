# GenomePuzzle Assessment Dataset Plan

Status: accepted implementation plan
Last updated: 2026-07-23

## Implementation status

- [x] Record the cross-repository assessment plan.
- [x] Add TOML release-spec loading and validation.
- [x] Add stable salted public IDs and row-order-independent sample seeds.
- [x] Add paired public/private manifest writing and leakage checks.
- [x] Add contract tests for identity, checksums, sample counts, and privacy.
- [x] Expose release building and validation through the main CLI.
- [x] Migrate genotyping out of `eqa-test.py`.
- [x] Migrate outbreak packaging and read-level implants out of `eqa-test.py`.
- [x] Integrate final short-read assets with release packages.
- [x] Integrate final hybrid assets with release packages.
- [ ] Add full external-tool integration and calibration runs.

## Purpose

`genomepuzzle` is the generation and truth-management side of the GHRU
Puzzles assessment platform. It must produce reproducible microbial genomics
exercises without exposing the source-to-sample mapping or expected answers to
participants.

The paired [`ghrupuzzle`](https://github.com/ghruproject/ghrupuzzle)
repository is responsible for publishing participant-facing assets,
authentication, submissions, scoring, review, and certificates. This
repository must not contain portal-specific upload or user-management logic.

## Exercise model

The platform has four participant-facing exercises:

1. Kleborate-based genotyping
2. Short-read de novo assembly
3. Hybrid assembly
4. Outbreak investigation, including phylogenetic reconstruction

Every exercise can have a permanently available practice release and a timed
challenge release. Both release types use the same generation and validation
contract. Access and answer-release policy are enforced by `ghrupuzzle`.

## Design principles

- A versioned release specification is the source of truth.
- Generation is deterministic wherever the underlying tools allow it.
- Public identifiers never reveal source accessions or implanted problems.
- Source provenance and expected results remain available in private outputs.
- Each sample receives at most one declared implant.
- Expected results are measured from the final participant-facing files.
- Partial output directories are never treated as completed releases.
- A dataset cannot be published until structural and biological validation
  passes.

## Canonical identifiers

Each generated record has three separate identifiers:

- `release_id`: identifies the assessment round and dataset version.
- `sample_id`: anonymous identifier used in filenames, public manifests,
  participant sheets, submissions, and scoring.
- `source_id`: private accession, SRA run, assembly, or outbreak tree tip.

The private manifest contains the `sample_id` to `source_id` mapping. Public
artifacts contain only `sample_id`.

Public sample identifiers must:

- be independent of input row order;
- be independent of implant type;
- contain no source accession or recognisable source name;
- be collision-checked within and across a release;
- remain stable when a release is resumed; and
- be derived using a release-specific private salt or an explicit ID registry.

Per-sample random seeds must be derived from stable release and source
identities rather than a row index.

## Release specifications

Release specifications should be stored as versioned TOML files, for example:

```text
releases/2026-round-1/
  typing-practice.toml
  typing-challenge.toml
  assembly-practice.toml
  assembly-challenge.toml
  hybrid-practice.toml
  hybrid-challenge.toml
  outbreak-practice.toml
  outbreak-challenge.toml
```

A specification records:

- schema version;
- release ID;
- exercise and practice/challenge mode;
- source accessions, reference genome, or source tree;
- intended sample count;
- public-ID registry or private salt reference;
- master random seed;
- simulator parameters;
- pinned tool and container versions;
- explicit implant assignments and severity;
- expected public columns;
- validation thresholds; and
- the answer-key schema required by the portal scorer.

Sample selection and implant assignment should be explicit in a frozen
challenge specification. A generator may help propose a balanced selection,
but the resolved selection must be recorded before release generation.

## Output contract

Each completed dataset is packaged as:

```text
release/
  public/
    dataset_manifest.json
    sample_sheet.csv
    checksums.sha256
    files/
  private/
    answer_sheet.csv
    implant_manifest.json
    provenance.json
    validation_report.json
```

### Public manifest

The public manifest contains only:

- schema version and release ID;
- exercise and mode;
- release metadata safe for participants;
- anonymous sample IDs;
- relative participant-facing filenames;
- file sizes and SHA-256 checksums; and
- participant sheet metadata.

It must not contain:

- source accessions or source names;
- implant types or severities;
- contaminant identities;
- private answer URLs;
- expected results; or
- private storage locations.

### Private provenance

Private provenance records, per sample:

- public sample ID;
- source identifier and source asset checksum;
- simulator inputs and per-sample seed;
- implant type and exact parameters;
- contaminant source, where applicable;
- generated filenames and checksums;
- pre- and post-implant QC measurements;
- expected analysis results;
- generator Git commit;
- dependency lock hash; and
- external tool and container versions.

The answer sheet must use a versioned, machine-readable schema suitable for
automatic scoring. Accepted aliases or equivalent answer representations
should be explicit rather than embedded in portal code.

## Troublesome samples

Practice and challenge datasets should each contain a small, known number of
realistic troublesome samples. The default target is two per dataset, but the
release specification controls the exact number.

| Exercise | Initial implant set | Expected evidence |
| --- | --- | --- |
| Genotyping | Mixed-contig assembly; fragmented or locus-disrupted assembly | Kleborate calls plus assembly QC |
| Short-read assembly | Low coverage; cross-species contamination | Read QC, assembly statistics, and taxonomy |
| Hybrid assembly | Weak or poor long reads; one-modality contamination or mismatch | Separate short/long QC and cross-modality agreement |
| Outbreak/phylogeny | Low-coverage isolate; contaminated or mixed isolate | QC exclusion decision and preserved cluster truth |

Implants should test biological analysis and QC judgement. Malformed gzip
files, null-byte corruption, and similar transport failures may remain
optional pipeline-resilience cases, but are not the core proficiency
assessment.

## Workflow-specific requirements

### Genotyping

- Move generation out of the deprecated `eqa-test.py` workflow.
- Produce anonymised FASTA inputs through the shared release framework.
- Apply assembly-level implants before calculating expected results.
- Run a pinned Kleborate version against the final FASTA files.
- Record raw Kleborate output and normalised scoring fields privately.
- Treat ambiguous or uncallable loci as explicit expected outcomes.

### Short-read assembly

- Retain one-implant-per-sample planning.
- Replace sequential names such as `sample01` with canonical release IDs.
- Make coverage, quality, contamination, and subsampling parameters explicit
  in the specification.
- Record paired-read counts and QC measurements after implantation.
- Keep participant sheets and private truth derived from the same sample
  context.

### Hybrid assembly

- Use the `long hybrid` workflow rather than the `rapid` compatibility bridge.
- Make Badread and all subsequent sampling deterministic where supported.
- Keep short-read and long-read provenance separate.
- Support implants affecting only one modality, allowing participants to
  detect disagreement between short- and long-read evidence.
- Validate that all three expected FASTQ files exist for every sample before
  completing a release.

### Outbreak and phylogeny

- Move generation out of `eqa-test.py` into a first-class module and CLI.
- Freeze the reference, source tree, cluster truth, epidemiological metadata,
  and simulator configuration in the release specification.
- Preserve a private tree-tip to public-sample mapping.
- Apply read-level implants only after clean outbreak reads have been
  generated.
- Maintain cluster truth separately from the expected QC decision, because a
  troubled isolate may need to be excluded from interpretation.

## Reproducible and restartable generation

The intended interface is:

```bash
genomepuzzle release build \
  --spec releases/2026-round-1/hybrid-practice.toml

genomepuzzle release validate \
  --release-dir generated/2026-round-1/hybrid-practice
```

Generation should:

1. Validate the specification and source table.
2. Resolve and checksum source assets.
3. Generate clean inputs in a staging area.
4. Apply the resolved implant plan.
5. Run expected reference analyses.
6. Build public and private manifests.
7. Run structural and biological validation.
8. Atomically mark or move the release as complete.

A stage may be reused only when its input hashes, parameters, tool versions,
and expected outputs match a recorded stage manifest. Directory existence
alone is never sufficient.

## Validation requirements

Structural validation must confirm:

- all required files exist;
- FASTA and FASTQ files can be parsed;
- paired FASTQ record counts agree;
- sample IDs are unique;
- filenames, sample sheets, and manifests form a one-to-one mapping;
- file sizes and SHA-256 checksums match;
- the configured implant count is correct;
- no sample has multiple declared implants;
- no private identifier or answer field appears publicly; and
- no temporary or source-download files appear in the public package.

Biological validation must confirm:

- normal controls meet expected QC;
- each implant produces the intended measurable signal;
- expected answers were calculated from final output files;
- typing results use the pinned Kleborate workflow;
- assembly and hybrid QC thresholds are calibrated against controls; and
- outbreak cluster truth remains consistent after public-name substitution.

Small fixtures should exercise contracts in normal CI. Full integrations using
external bioinformatics tools should run separately, including through Slurm
where appropriate.

## Interface with `ghrupuzzle`

`ghrupuzzle` consumes only completed release packages. It must not:

- recalculate public sample names;
- infer filenames from answer-sheet rows;
- run simulation tools;
- reconstruct source provenance; or
- publish private files.

The shared boundary is a versioned manifest schema. Changes to that schema
require contract tests in both repositories.

The website publisher receives:

- the public package for participant distribution; and
- the private answer/provenance package through an explicitly private
  administrative path.

## Implementation sequence

1. Define release, public-manifest, private-manifest, and answer-key schemas.
2. Add the stable ID registry and per-sample seed derivation.
3. Add shared release packaging and validation.
4. Modernise genotyping generation.
5. Modernise outbreak and phylogeny generation.
6. Align short-read generation with the shared contract.
7. Align hybrid generation with the shared contract.
8. Add small deterministic fixtures and integration jobs.
9. Generate and calibrate all practice releases.
10. Freeze and generate challenge releases only after the practice pilot.

## Completion criteria

Work in this repository is complete when:

- one command builds each specified dataset;
- rebuilding with the same source assets and environment reproduces the same
  participant files and checksums;
- reordering an input table does not alter public identities;
- every release contains a complete public/private package;
- public artifacts contain no answer or source leakage;
- each troublesome sample is detected by the expected validation evidence;
- all structural and biological checks pass before publishing; and
- `ghrupuzzle` can publish the package without renaming or inferring files.

## Explicit non-goals

- User accounts, submissions, and certificates do not belong in this
  repository.
- The generator does not execute participant-submitted workflows.
- The first release does not attempt to cover every possible sequencing
  failure.
- Challenge selection does not need to be generated dynamically at request
  time.
