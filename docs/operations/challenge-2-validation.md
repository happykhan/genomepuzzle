# Challenge 2 blinded validation

Challenge 2 was independently analysed from participant-visible files before
its private truth was opened. This page preserves the evidence that the
scientific exercises are solvable and the deliberate faults are categorical.

## Outcome

| Check | Result |
| --- | ---: |
| `qc_status` | 48/48 correct |
| `failure_reason` | 48/48 correct |
| Species on passing samples | 35/35 correct |
| Assembly scored fields | 27/27 correct |
| Hybrid scored fields | 27/27 correct |
| Outbreak non-partition fields | 50/50 correct |
| Outbreak passing-sample relationships | 91/91 correct |
| Typing scored fields | 75/76 correct |

The only miss was representational rather than biological: the analyst
submitted a blank `wzi` where Kleborate and the answer key used the literal
`-`. GHRU Puzzles should normalise these equivalent unavailable values or
offer a controlled form value.

All 13 deliberately troublesome challenge samples received the exact intended
participant-facing failure reason.

## Blinding

The analyst workspace exposed only each release's `public/` directory. The
four result CSVs were made read-only and SHA-256 sealed before any answer key,
implant manifest, scoring policy or private validation report was read.

| Exercise | Sealed submission SHA-256 |
| --- | --- |
| Typing | `6c00aaf7bd70ad71aeff1e13fe9df44903ee37c960ea0ae49102e52f598b57d5` |
| Assembly | `6dcbbe0b539f65171c86382c682fe253385dec8aa85b0004f909e2b02ecf4299` |
| Hybrid | `68faee4dda3fb88f882dd75059544ca7a1264109af8c99a6281aea2cfcb86792` |
| Outbreak | `8f382fd481311bc28a4985fddc008ff9c2f353868b42cc87187e4c3919620b83` |

## Environment and compute

The isolated Pixi environment contained Kleborate 3.2.4, SeqKit 2.13.0,
SPAdes 4.0.0, Unicycler 0.5.1, Snippy 4.6.0, IQ-TREE 3.1.2, Mash 2.3,
SKA 0.5.1, fastp 0.24.0 and QUAST 5.3.0.

Biological work ran through SLURM:

| Work | Successful jobs |
| --- | --- |
| Typing statistics and Kleborate | `25472733` |
| Short-read assembly and Kleborate | `25472734`, `25472735` |
| Hybrid assembly and Kleborate | `25472736`, `25472737` |
| Outbreak assembly and Kleborate | `25472738`, `25472739` |
| Outbreak Snippy array and core tree | `25478467`, `25478468` |

Early orchestration attempts failed before useful biological computation
because commands lacked `pixi run` or analyst-workspace symlinks were wrong.
The final jobs above are the evidence-bearing runs.

One already-failed contaminated hybrid task was stopped after 39 minutes when
all valid assemblies had completed and public evidence was conclusive.
Conditional scoring correctly avoids demanding expensive assembly statistics
for a sample already rejected on QC.

## Observed signals

### Typing

- a literal zero-byte assembly;
- 5,748 contigs with N50 1,000 for extreme fragmentation;
- an approximately 10.27 Mb mixed assembly for contamination.

Seven passing assemblies gave coherent Kleborate species, ST, K/O antigen and
carbapenemase calls.

### Short-read assembly

- a zero-byte R2;
- 41 retained read pairs;
- an approximately 10.19 Mb contaminated assembly with a Kleborate total-size
  warning.

Seven passing read sets produced plausible SPAdes assemblies and the expected
species.

### Hybrid assembly

- a missing long-read role;
- a discordant short/long Mash distance of 0.239495 with only 86/10,000 shared
  hashes, versus approximately 0.074–0.077 and 1,552–1,671 shared hashes for
  coherent samples;
- a contamination sample with 928,320 read pairs and 6,320 long reads, well
  above the ordinary sample range.

The discordant sample still produced a superficially plausible short-read-led
assembly. This usefully requires participants to compare modalities rather
than treating assembly completion as proof of coherence.

### Outbreak

- a zero-byte R1;
- a missing R2;
- 83 retained read pairs;
- an approximately 10.22 Mb contaminated assembly.

The 14 passing samples produced a 129-site core SNP alignment. Within inferred
clusters, distances were 9–10 SNPs; between clusters they were 49–50 SNPs.
The three-way partition matched truth exactly.

## Conclusion

The Challenge 2 fault vocabulary is fair and unambiguous. Variable retained
counts of 41 and 83 demonstrate that the assessment does not depend on a
recognisable fixed ten-read signature. Future production releases should
repeat the blinded process in
[Blinded release acceptance](blinded-acceptance.md).
