# 2026 pilot releases

This directory freezes the first production-calibration cohort for all four
exercise types:

| Specification | Exercise | Cohort |
| --- | --- | --- |
| `typing-practice.toml` | Kleborate genotyping | 6 controls, 2 troublesome assemblies |
| `assembly-practice.toml` | Short-read assembly | 6 controls, 2 troublesome read sets |
| `hybrid-practice.toml` | Hybrid assembly | 6 controls, 2 troublesome read sets |
| `outbreak-practice.toml` | Phylogeny and outbreak | 6 controls, 2 troublesome read sets |

The troublesome samples are deliberately anonymous in participant artifacts.
Their source mapping and implant parameters remain in the private release
manifests.

## Build on SLURM

Create and submit a persisted generation and validation workflow from the
repository root:

```bash
pixi run genomepuzzle release build \
  --spec releases/2026-pilot/assembly-practice.toml \
  --output-dir generated/2026-pilot/assembly-practice \
  --partition short
```

Replace the specification and output directory for each exercise. Monitor or
resume the recorded workflow rather than issuing an untracked command:

```bash
pixi run genomepuzzle release status \
  --plan generated/2026-pilot/assembly-practice/build/plan.json

pixi run genomepuzzle release resume \
  --plan generated/2026-pilot/assembly-practice/build/plan.json
```

Generation, assembly, Kleborate and phylogenetic calibration are cluster work.
Do not run them on a login node. Pixi supplies the pinned tools; there is no
repository `bin/` environment.

For small development cohorts on a standalone Linux laptop, use the direct
`release generate-reads` or `release generate-outbreak` commands. They run
sequentially without a SLURM allocation.

## Acceptance

A structurally valid release is not automatically an educationally useful
release. Before publication:

1. require a successful dependent validation job and `COMPLETE.json`;
2. run participant-style calibration from the public files;
3. compare each troublesome sample with a clean control;
4. confirm the problem is detectable but does not merely corrupt the file;
5. record the job IDs, bundle digest, measurements and reviewer decision.

The completed 2026 pilot evidence and the calibration command are documented
in [Pilot calibration](../../docs/operations/pilot-calibration.md).

Do not publish a typing build made with `--skip-analysis`: its private answer
key is marked `pending_kleborate`.
