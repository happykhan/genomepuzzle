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

## Production build record

The 23 July 2026 production build completed through SLURM, passed independent
GenomePuzzle validation and passed a dry run through the GHRUPuzzles
manifest-driven publisher.

| Release | Generate / validate jobs | Bundle SHA-256 |
| --- | --- | --- |
| `typing-practice` | `25214661` / `25214662` | `2f4ad43e874bce44cd7a0999840f580ba9bfa3bed8941409e59c01aef9e87778` |
| `typing-challenge` | `25214659` / `25214660` | `484b7f3a55ea67b1cd1f39b9766c7e5c1e078e788a8e24757a7da99c376c2283` |
| `assembly-practice` | `25214649` / `25214650` | `0ca1d3a5d7dcf5013fbce7ae48c85d502abc20f1cc8e4728fa643e5c5747e087` |
| `assembly-challenge` | `25214647` / `25214648` | `d76700fd2f026cfdb871256a91a5935b40987fe0ce21f8b3551bf71238bd237c` |
| `hybrid-practice` | `25214653` / `25214654` | `84e6056354e5f19a2bcd578587f618c3b480146f1ed54da8f07bb3d17ba286cd` |
| `hybrid-challenge` | `25214651` / `25214652` | `f32d6e6184a24394f6328036f2f4e4f42b9effa0cbd9acec0283d3e55207b9f5` |
| `outbreak-practice` | `25228679` / `25228680` | `7e943d7801569e9bb24ae63e2c07b4b017ed5895983103d6ddf0e7833895a416` |
| `outbreak-challenge` | `25228681` / `25228682` | `92ceb83b567d35c86d375c42a649ad7c80a45da96505af5e9e803b342b5a6240` |
