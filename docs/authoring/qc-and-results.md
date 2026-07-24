# QC failures and participant results

GenomePuzzle uses one canonical participant-results vocabulary across every
exercise. A column has the same name, meaning and value format wherever it
appears, but each exercise includes only results relevant to its own task.
Assembly and outbreak do not silently become additional genotyping exercises.
CSV column names are lower-case `snake_case`; the website may render
human-friendly labels such as **ST** and **K locus**.

## QC decision

Every exercise begins with:

```csv
sample_id,qc_status,failure_reason
```

`qc_status` is exactly `PASS` or `FAIL`. `failure_reason` is `NONE` for a
passing sample and one exercise-supported categorical reason for a failing
sample. There is no `REVIEW` state.

The assessment deliberately uses catastrophic failures, not threshold
near-misses. Two competent analysts should reach the same QC decision without
debating whether a dataset is just usable.

## Public reasons and private faults

The participant reports a stable `failure_reason`. The private implant
manifest retains the exact `fault_type` and achieved evidence:

```yaml
qc_status: FAIL
failure_reason: TOO_FEW_READS
fault_type: TRUNCATED_R2_TO_10_READS
validation:
  observed_r1_reads: 100000
  observed_r2_reads: 10
```

Several precise implants may map to one reportable reason. This keeps the
participant vocabulary small while preserving reproducibility.

There is deliberately no `MATE_COUNT_MISMATCH` reason. Read trimming can
legitimately create singletons, making ordinary mate-count differences a poor
proficiency question.

## Failure vocabulary

| Reason | Exercises | Required interpretation |
| --- | --- | --- |
| `NONE` | All | Sample passes QC |
| `EMPTY_FILE` | All | A required supplied sequence file contains no records |
| `MISSING_MATE` | Assembly, hybrid, outbreak | One paired-read role is absent |
| `TOO_FEW_READS` | Assembly, hybrid, outbreak | An absurdly small short-read dataset, normally about ten pairs or a severely truncated mate |
| `LOW_COVERAGE` | Assembly, hybrid, outbreak | Measured short-read coverage is catastrophically low, normally no more than approximately 1× |
| `MISSING_LONG_READS` | Hybrid | The long-read role is absent |
| `TOO_FEW_LONG_READS` | Hybrid | Only an absurdly small number of long reads is supplied |
| `CONTAMINATED` | All | A clearly different species contributes 30–90% of sequence, normally about 50% |
| `WRONG_ORGANISM` | All | The material is effectively another organism rather than a mixture |
| `DISCORDANT_READ_SETS` | Hybrid | Short and long reads represent different organisms |
| `EXTREME_FRAGMENTATION` | Typing | An assembly is so fragmented that it is categorically unsuitable |

At effectively 100% replacement, use `WRONG_ORGANISM`, not `CONTAMINATED`.
The contaminant must be a clearly different species and the achieved mixture
must be validated from the final participant material.

## Conditional scoring

QC fields are always scored. Downstream analytical fields are scored only for
samples whose expected `qc_status` is `PASS`.

| Expected QC | Scored fields |
| --- | --- |
| `PASS` | `qc_status`, `failure_reason`, and configured analytical fields |
| `FAIL` | `qc_status` and `failure_reason` |

For outbreak exercises, failed samples are excluded from partition scoring.
Their biological cluster truth remains in private provenance, but a cluster
answer cannot rescue or harm an otherwise correct failure decision.

The conditions are emitted in `private/scoring_policy.json`; the website must
consume them rather than duplicating exercise-specific logic.

## Canonical field formats

| Field | Canonical format |
| --- | --- |
| `species` | Full canonical taxon call; scoring may accept a species-level form for a subspecies call |
| `st` | Number without an `ST` prefix; importers may accept and normalise the prefix |
| `bla_carb` | Semicolon-separated, order-independent gene calls |
| `contig_count` | Integer |
| `total_length` | Integer bases |
| `n50` | Integer bases |
| `longest_contig` | Integer bases |
| `notes` | Optional free text; never automatically scored |

Unknown analytical values are empty. Do not mix `NA`, `N/A`, `ND`, `unknown`
and `-` in participant answers.

Assembly statistics are collected consistently but are not exact-string
scored: different valid assemblers need range/tolerance scoring or manual
review. `assembler` and `notes` are provenance supplied by the participant and
are never automatically scored.

## Exercise columns

Genotyping:

```csv
sample_id,qc_status,failure_reason,species,st,k_locus,capsule_type,wzi,o_locus,o_type,bla_carb,notes
```

Short-read assembly:

```csv
sample_id,qc_status,failure_reason,species,assembler,contig_count,total_length,n50,longest_contig,notes
```

Hybrid assembly:

```csv
sample_id,qc_status,failure_reason,species,assembler,contig_count,total_length,n50,longest_contig,notes
```

Outbreak investigation:

```csv
sample_id,qc_status,failure_reason,species,cluster,notes
```
