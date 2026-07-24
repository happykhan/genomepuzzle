# Troubleshooting

## `sbatch` is unavailable

`release build`, `submit`, `resume`, `status` and `logs` are the audited SLURM
workflow and require a configured cluster. On a Linux laptop, use the direct
exercise command such as `release generate-reads`, then run `release validate`.
GenomePuzzle never silently changes execution mode after a failed submission.

## Public IDs change between runs

Use the same `release_id`, `identity_key` and private salt. For a stable
practice release, declare and freeze `public_id`. Do not use sample row order as
identity.

## A source assembly cannot be found

For typing, assembly and hybrid releases, place it under `[inputs].source_dir`
as `<source_id>.fasta`, `<source_id>.fna`, `<source_id>.fa` or an extensionless
file. Paths are relative to the TOML specification.

## Validation says an implant lacks evidence

The requested implant did not produce the expected measurable result, or the
validation record is absent. Inspect the generation log and
`private/validation_report.json`. Adjust the implant parameters and create a
new intentional release version; do not relabel the sample as troublesome
without evidence.

## Validation finds a source identifier

Check participant filenames, FASTA/FASTQ headers, public metadata and
instructions. Source accessions belong only in private provenance.

## A generation job failed

```bash
pixi run genomepuzzle release status --plan <release>/build/plan.json
pixi run genomepuzzle release logs --plan <release>/build/plan.json
sacct -j <job-id> --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

Correct the underlying input or configuration, then resume. Do not delete the
whole release directory: its build state is the audit record.

## `COMPLETE.json` is missing

The release is not publishable. Check that generation finished, validation
passed and every required public/private artifact exists. Directory existence
alone never means a release is complete.

## Pixi cannot reproduce the environment

Use the checked-in lock file:

```bash
pixi install --locked
```

If a dependency genuinely needs to change, update `pixi.toml`, regenerate
`pixi.lock`, run the full test suite and treat subsequent datasets as produced
by a new software environment.
