# Pilot calibration

Release validation proves the file and manifest contract. Pilot calibration
separately asks whether the participant-facing data exercise the intended
scientific judgement.

## Method

Calibration must run under SLURM and must read only `public/files`. For
short-read and hybrid exercises, compare one clean control with every
troublesome sample:

```bash
pixi run python -m genomepuzzle.pilot_review \
  --release-dir generated/2026-pilot/assembly-practice-v2 \
  --reference-dir generated/2026-pilot/typing-sources \
  --sample-id Sample_AP001 \
  --sample-id Sample_AP007 \
  --sample-id Sample_AP008 \
  --threads 8 \
  --memory-gb 32
```

The command refuses to run outside a SLURM allocation. It records:

- read counts, bases, lengths, mean quality and GC;
- Mash read screens against the frozen source set;
- participant-style SPAdes assemblies and assembly statistics;
- Kleborate calls from sample-labelled assembly copies; and
- a Mashtree relatedness smoke test for an outbreak release.

Results are written beneath `build/calibration/`, which is intentionally
excluded from the sealed bundle digest. `report.json` records the SLURM job,
Git commit, Pixi lock digest and exact commands. A retry reuses completed
assemblies.

## 2026 pilot evidence

All releases contain eight samples: six unmodified controls and two implanted
problems. The directory suffix is a build attempt, while `release_id` remains
the stable assessment identity.

| Exercise | Release directory | Production job IDs | Bundle SHA-256 |
| --- | --- | --- | --- |
| Typing | `typing-practice-v3` | `25089408` / `25089409` | `bea9eb68e158ea432fe6c0063bdd0c47ffa35449eb3a6dfa54e9dd5865f5c7db` |
| Short-read assembly | `assembly-practice-v2` | `25087440` / `25087441` | `b743c08b0c1d58c1a098e5d6c83e1f6a129be392b0dff0d53104427a4bd19270` |
| Hybrid assembly | `hybrid-practice-v6` | `25131352` / `25145087` / `25145088` | `9c50c5bc6066889b53d56ca9d7e849170959809d864e94eb8bbcd029cd91995e` |
| Outbreak | `outbreak-practice-v2` | `25087054` / `25087055` | `129e17923e81143d9de608c2840999c50510a370ed276ad5b937519422d1b757` |

### Typing

`Sample_TP007` is a fragmented assembly. It contains 5,601 1 kb contigs
instead of the six source contigs. Kleborate reports an N50 warning,
`ST11-3LV`, and no `wzi` call. This is a useful obvious QC failure while the
FASTA itself remains valid.

`Sample_TP008` contains exactly 677,554 contaminant bases, an achieved fraction
of 0.11999994. It remains superficially typeable as `ST101` without a generic
QC warning, but reports the conflicting carbapenemases `KPC-3;OXA-48` and an
altered O locus. This tests review of biologically discordant results rather
than file handling.

Decision: accept both implants.

### Short-read assembly

Calibration job `25112015` compared `Sample_AP001`, `Sample_AP007` and
`Sample_AP008`.

| Sample | Implant | Read pairs | Contigs | Total bases | N50 | Kleborate |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| `Sample_AP001` | control | 759,800 | 196 | 5,631,922 | 190,559 | `ST14`, no QC warning |
| `Sample_AP007` | low coverage | 112,480 | 1,226 | 5,497,638 | 7,875 | `ST11-1LV`, N50 warning |
| `Sample_AP008` | mixture | 905,198 | 455 | 6,404,475 | 90,322 | `ST101`, no QC warning |

The low-coverage sample is approximately 6× and fails assembly QC without
being malformed. The mixed sample has a plausible assembly and no Kleborate
QC warning, while its Mash screen contains two equally strong source signals
and its assembly carries `KPC-3;OXA-48`. It therefore tests whether a learner
looks beyond N50.

Decision: accept both implants.

### Hybrid assembly

The final cohort uses 30× short reads and 10× long reads. Simulation job
`25131352` produced all read tracks in 32 minutes. It exposed a strict-header
validation edge case after generation: Badread's generic `junk_seq` and
`random_seq` labels are anonymous but did not contain the public sample ID.
The policy was corrected, future simulations disable those categories, and
the retained reads were sealed by job `25145087` in 72 seconds. Independent
validation job `25145088` passed all 8 samples and 24 participant files.
Calibration job `25145243` completed in 19 minutes 49 seconds.

| Sample | Implant | Read pairs | Long-read bases | Contigs | Total bases | N50 | Kleborate |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| `Sample_HP001` | control | 569,850 | 57,006,061 | 20 | 5,696,400 | 1,610,515 | `ST14`, no QC warning |
| `Sample_HP007` | low long coverage | 559,725 | 8,581,538 | 130 | 5,566,300 | 287,949 | `ST11`, no QC warning |
| `Sample_HP008` | mixture | 678,598 | 67,615,985 | 211 | 6,426,742 | 4,196,700 | `ST101`, no QC warning |

The 1.5× long-read sample still assembles and types correctly, but loses the
main contiguity benefit of hybrid assembly: its N50 falls from 1.61 Mb to
288 kb. The mixed sample is deliberately deceptive. Its N50 rises to 4.20 Mb,
but its assembly is inflated and fragmented, Mash gives essentially complete
support to both source strains, and Kleborate reports
`KPC-3*;OXA-48`. Reporting N50 alone would therefore produce the wrong QC
decision.

Decision: accept both implants.

### Outbreak

Calibration job `25112016` assembled and typed all eight isolates.

The six controls have 192–208 contigs and N50 values of 183,908–194,740 bp.
`Sample_OP007` has only 85,735 read pairs (approximately 4.5×), 2,784 contigs,
N50 3,008 bp, an N50 warning and an incomplete `ST14-1LV` call. It should be
excluded before inference.

`Sample_OP008` is a cross-cluster read mixture. Its consensus assembly remains
normal-looking at 184 contigs, N50 188,245 bp and `ST14` with no QC warning.
The mixture must therefore be detected from read-level minor alleles rather
than from a broken consensus assembly.

The native truth has 20 cluster-shared SNPs plus five private SNPs per
isolate. Mashtree completes as a packaging and identity smoke test but rounds
most of these distances to zero, so it is not an acceptance test for the
expected clusters. Cluster acceptance requires a reference-based core-SNP
workflow and read-level mixture review.

Decision: accept both implants; do not use Mash-distance topology as the
scored outbreak truth.

## Review rule

A troublesome sample is accepted when:

- all participant files pass structural validation;
- its intended signal is reproduced from the public files;
- at least one routine result can look plausible enough to require judgement;
- the expected QC decision is defensible; and
- the sample ID remains traceable through every analysis output.

Reject or recalibrate an implant that only truncates a file, makes every tool
crash, is indistinguishable from controls, or is identified only by a leaked
source label.

## Generator findings

Production proving changed the generator itself:

- ART and Badread copy reference contig names into read headers. GenomePuzzle
  now gives each simulator a temporary reference containing only public,
  sample-scoped contig names.
- Generated reads are validated and copied during sealing rather than
  decompressed and recompressed solely to replace their headers.
- Laptop defaults are 30× short reads and 10× long reads. A direct
  `release generate-reads` command runs sequentially without SLURM; the
  production workflow uses bounded sample parallelism when CPUs are allocated.
- Badread junk and random reads are disabled for new cohorts. The privacy
  validator still recognises their generic labels as anonymous so an
  interrupted older run can be recovered safely.
