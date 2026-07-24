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
| `assembly-practice` | `25214649` / `25214650` | `0ca1d3a5d7dcf5013fbce7ae48c85d502abc20f1cc8e4728fa643e5c5747e087` |
| `hybrid-practice` | `25214653` / `25214654` | `84e6056354e5f19a2bcd578587f618c3b480146f1ed54da8f07bb3d17ba286cd` |
| `outbreak-practice` | `25228679` / `25228680` | `7e943d7801569e9bb24ae63e2c07b4b017ed5895983103d6ddf0e7833895a416` |
| `challenge-2-typing` | `25404244` / `25404245` | `795d9d3dd95464bac21387bdb62f1d236f218f25fce9ea28ad871561fd2d52ee` |
| `challenge-2-assembly` | `25404246` / `25404247` | `5344a4c4f1cdfa1b7fc4a29bea14aac488a60f10a5b37e3f20d317357b4d01a1` |
| `challenge-2-hybrid` | `25404248` / `25404249` | `442dc01cb7c3991f3a057905ca621863a96f39346fe0d29b8a57e25bf59a25c4` |
| `challenge-2-outbreak` | `25404250` / `25404251` | `246272a9ef3d0b0e571a9e9b1d6fbdc1bc4abde93c34d1c241511ae734dd6fcc` |

The four original date-prefixed challenge builds were superseded before
publication. Their replacements use the participant-facing `Challenge 2`
name and are recorded after regeneration.

All eight releases were published to R2 on 24 July 2026. Practice participant
artifacts are in the public practice bucket; practice answer material and the
complete Challenge 2 contracts are in the private bucket. A post-upload audit
matched every remote object size and SHA-256 metadata value, confirmed HTTP
200 participant downloads for all practice exercises, and found no Challenge
2 prefix or private-answer path in the public bucket.
