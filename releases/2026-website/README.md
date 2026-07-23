# 2026 website datasets

This directory defines the first production datasets consumed by GHRUPuzzles.
Each exercise has one permanently available practice release and one timed
challenge release.

| Exercise | Practice cohort | Challenge cohort |
| --- | ---: | ---: |
| Genotyping | 4 (3 control, 1 troublesome) | 10 (8 control, 2 troublesome) |
| Short-read assembly | 4 (3 control, 1 troublesome) | 10 (8 control, 2 troublesome) |
| Hybrid assembly | 4 (3 control, 1 troublesome) | 10 (8 control, 2 troublesome) |
| Phylogeny/outbreak | 6 (5 control, 1 troublesome) | 18 (16 control, 2 troublesome) |

Practice IDs are fixed so tutorials and repeat attempts remain comparable.
Challenge IDs are derived from `GENOMEPUZZLE_ID_SALT`; never commit or publish
that salt. The challenge time window is configured by GHRUPuzzles when the
completed release is registered, not in these generator specifications.

Source FASTAs are cached outside generated releases:

```text
cache/assemblies/
```

Build production releases through the normal SLURM workflow. For example:

```bash
export GENOMEPUZZLE_ID_SALT='use-a-private-random-value'
pixi run genomepuzzle release plan \
  --spec releases/2026-website/hybrid-challenge.toml \
  --output-dir generated/2026-website/hybrid-challenge
pixi run genomepuzzle release submit \
  --plan generated/2026-website/hybrid-challenge/build/plan.json
```

Only a release containing `COMPLETE.json` and passing independent validation
is ready for the GHRUPuzzles manifest-driven publisher.
