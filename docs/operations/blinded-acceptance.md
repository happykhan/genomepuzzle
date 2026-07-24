# Blinded release acceptance

Automated contract validation proves structural consistency. It does not prove
that an independent analyst can recognise faults, obtain the expected results
or complete the exercise with ordinary tools. Production releases therefore
require a blinded acceptance run.

## Roles

Use two logical roles even if one person performs them at different times:

- the release custodian can access the complete public/private bundle;
- the blinded analyst receives only an isolated copy or symlink of `public/`.

Do not rely on personal restraint while private files remain visible in the
same working directory.

## Procedure

1. Freeze a candidate release and run:

   ```bash
   pixi run genomepuzzle release validate \
     --release-dir generated/round/exercise \
     --require-complete
   ```

2. Create a new analyst directory outside the release tree.
3. Expose only public manifests, instructions, schemas and participant files.
4. Create and lock a separate Pixi analysis environment.
5. Record tool versions, commands, SLURM scripts, job IDs and outputs.
6. Run all dataset-scale biological analysis through SLURM. Lightweight file
   counts, CSV preparation, manifest reads and checksum operations may run on
   the login node.
7. Use public evidence to complete every row in the advertised submission
   schema.
8. Treat QC conditionally:
   - `PASS` requires the exercise's analytical answers;
   - `FAIL` requires the categorical failure reason;
   - unavailable analytical values are left blank and are not inferred from
     provenance.
9. For outbreak work, retain the analyst's natural cluster labels. Do not try
   to guess organiser label names.
10. Validate the CSV against the public schema and exact public sample set.
11. Calculate SHA-256, record a timestamp and make the submission read-only.
12. Verify every seal. Only then grant the comparison process access to
    `private/`.
13. Compare using `scoring_policy.json`, including `score_when`, normalizers,
    unordered lists and partition scoring.
14. Write a report containing scientific evidence, discrepancies, operational
    problems and all hashes.

## Acceptance criteria

A release is suitable for publication when:

- all intended PASS/FAIL decisions are independently defensible;
- catastrophic failure reasons are recovered without borderline thresholds;
- passing samples produce plausible analytical results;
- no public filename, header or metadata exposes source identity or truth;
- failed analytical fields are excluded according to the scoring policy;
- outbreak partitions are recoverable without requiring literal labels;
- every mismatch is understood and either corrected or explicitly accepted;
- the workflow is practical for the intended environment.

Do not accept a fault merely because it matches private truth. Two competent
analysts should reach the same categorical decision from the participant data.

## Evidence to retain

Retain, outside participant bundles:

- the locked Pixi environment;
- SLURM scripts and final job IDs;
- tool outputs used to make decisions;
- sealed analyst CSVs and their timestamps;
- the truth-comparison script and result;
- the written validation report.

Large intermediate assemblies may be removed after review if the retained
report, logs, final evidence and checksums are sufficient to reproduce the
decision.

## Laptop calibration

Practice cohorts should remain small enough for sequential laptop analysis.
Challenge-scale hybrid assembly may reasonably take hours when run
sequentially. Validate both that:

- the practice workflow is usable without a cluster; and
- challenge workflows parallelise cleanly through the supplied SLURM model.

The acceptance run may stop expensive analysis of an already categorical FAIL
sample once the evidence is recorded. It must never fabricate assembly
statistics for that sample.
