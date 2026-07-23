# SLURM workflow

GenomePuzzle uses a persisted two-stage workflow:

1. `generate` creates, analyses, packages and seals the release.
2. `validate` independently checks the sealed bundle with an `afterok`
   dependency.

There is no local fallback for heavy biological work.

## Commands

`release plan` resolves the specification and writes scripts without
submission. `release submit` submits an existing plan. `release build` performs
both operations.

```bash
pixi run genomepuzzle release plan \
  --spec releases/round/hybrid-practice.toml \
  --output-dir generated/round/hybrid-practice \
  --partition short

pixi run genomepuzzle release submit \
  --plan generated/round/hybrid-practice/build/plan.json
```

Monitor through GenomePuzzle:

```bash
pixi run genomepuzzle release status --plan <release>/build/plan.json
pixi run genomepuzzle release logs --plan <release>/build/plan.json
```

Use normal cluster tools for additional diagnosis:

```bash
squeue -j <job-id>
sacct -j <job-id> --format=JobID,State,ExitCode,Elapsed,NodeList
```

## Recorded state

The `build/` directory contains:

| Path | Purpose |
| --- | --- |
| `plan.json` | Frozen commands, samples, Git commit and lock digest |
| `scripts/*.sbatch` | Exact submitted job scripts |
| `stages/*.json` | Attempts, job IDs, timestamps, nodes and status |
| `logs/` | SLURM standard output and error |
| `work/` | Intermediate generation material |

Mutable `build/` state is excluded from the final bundle digest.

## Default resources

| Exercise/stage | CPUs | Memory | Time |
| --- | ---: | ---: | ---: |
| Typing generation | 4 | 16 GB | 6 hours |
| Outbreak generation | 4 | 16 GB | 6 hours |
| Assembly/hybrid generation | 8 | 32 GB | 12 hours |
| Contract validation | 1 | 4 GB | 1 hour |

The CLI accepts the cluster partition. Resource values are stored in the plan
and generated script, making the actual request auditable.

## Resume behaviour

```bash
pixi run genomepuzzle release resume --plan <release>/build/plan.json
```

Resume records a new attempt and preserves earlier logs and job IDs. A failed
generation retry clears partial public/private artifacts and the root
completion seal, but retains the workflow audit trail.
