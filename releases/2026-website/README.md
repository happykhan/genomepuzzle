# 2026 website datasets

This directory defines the first production datasets consumed by GHRUPuzzles.
Each exercise has one permanently available practice release and one timed
challenge release.

| Exercise | Practice cohort | Challenge cohort |
| --- | ---: | ---: |
| Genotyping | 5 (3 pass, 2 fail) | 10 (7 pass, 3 fail) |
| Short-read assembly | 5 (3 pass, 2 fail) | 10 (7 pass, 3 fail) |
| Hybrid assembly | 5 (3 pass, 2 fail) | 10 (7 pass, 3 fail) |
| Phylogeny/outbreak | 6 (4 pass, 2 fail) | 18 (14 pass, 4 fail) |

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

## Contract 2.1 production build record

The 24 July 2026 contract 2.1 releases completed through SLURM and passed both
their dependent validation jobs and a separate local validation. The first
seven plans are pinned to Git commit `4330039`; outbreak challenge is pinned
to `371d210`, which adds the regression-tested sample-sheet handling required
by its deliberate `MISSING_R2` fault.

| Release | Generate / validate jobs | Bundle SHA-256 |
| --- | --- | --- |
| `typing-practice` | `25442292` / `25442293` | `c9fcdb57e1648e58d660a49d48399654da363e8c0d53f6ca2a5c4434ce020c7b` |
| `assembly-practice` | `25442295` / `25442296` | `ad3123e7f5e966dbd44cd9a705f2436b326f74fe377b7b6764b390502701717b` |
| `hybrid-practice` | `25442298` / `25442299` | `464de50827d1e86eea482516eb574f499629560cb8e90369efc7df2d9efee0b9` |
| `outbreak-practice` | `25442304` / `25442305` | `c9a1901fa731225922ab21c59d8c076de92e63c65cb890faf94e4212e239b240` |
| `challenge-2-typing` | `25442307` / `25442308` | `8307a36d4b979293ee2ea084b0a68278dc754ebcdd2a3d27fe34b7a38a7d6a12` |
| `challenge-2-assembly` | `25442310` / `25442311` | `e5f4c5d26f50e5bd559d7038ed9f1642ed841015b0fca79432e0ceaebade8578` |
| `challenge-2-hybrid` | `25442313` / `25442314` | `add243400399180bfb31a2ea52b50b9df95bd98f7b72b516aaadd7b79a6bc506` |
| `challenge-2-outbreak` | `25451564` / `25451565` | `ac805d7b6541be07db7dfc7d06651e1ab3fea1d01b15e4ae8f0d171faf84dbcd` |

Manual fault inspection confirmed literal zero-byte and absent roles, exactly
ten read pairs where declared, 40–50% typing contamination, exact 50% read
contamination, and the intended hybrid short/long-read discordance. Public
trees contain no private source, fault or answer keys.

## Superseded production build record

The 23 July 2026 production build completed through SLURM, passed independent
GenomePuzzle validation and passed a dry run through the GHRUPuzzles
manifest-driven publisher. It uses the previous QC contract and is retained
only as provenance; it must not be registered on the website.

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

These eight superseded releases were published to R2 on 24 July 2026. Practice participant
artifacts are in the public practice bucket; practice answer material and the
complete Challenge 2 contracts are in the private bucket. A post-upload audit
matched every remote object size and SHA-256 metadata value, confirmed HTTP
200 participant downloads for all practice exercises, and found no Challenge
2 prefix or private-answer path in the public bucket.

The validated contract 2.1 releases recorded above supersede these objects and
must replace them in R2 before website registration.
