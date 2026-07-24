# Hybrid assembly

Hybrid releases combine ART short reads with Badread long reads. The separate
data modalities allow an assessment to test whether participants detect weak
long-read data or disagreement between evidence sources.

## Supported implants

| Implant | Important parameters | Intended signal |
| --- | --- | --- |
| `LOW_SHORT_COVERAGE` | `read_fraction` | No more than approximately 1× Illumina depth |
| `TRUNCATE_TO_READ_PAIRS` | `retained_read_pairs` | An explicit 1–100 short-read pairs |
| `MISSING_LONG_READS` | — | Long-read role is absent |
| `ZERO_BYTE_LONG_READS` | — | Long-read role is a literal zero-byte file |
| `TEN_LONG_READS` | — | Exactly ten long reads |
| `CONTAMINATED` | `contaminant_source_id`, `contamination_fraction` | A 30–90% different-species mixture in both modalities |
| `WRONG_ORGANISM` | `replacement_source_id` | Both modalities are another organism |
| `DISCORDANT_READ_SETS` | `contaminant_source_id` | Short and long reads are different organisms |

Normal long-read quantity defaults to `10x` and can be overridden with
`long_quantity`.

```toml
[[samples]]
source_id = "GCA_000000003.1"
implant = "CONTAMINATED"

[samples.implant_parameters]
contaminant_source_id = "GCA_000000004.1"
contamination_fraction = 0.50
long_quantity = "20x"
```

A contaminated hybrid sample is altered in both its short- and long-read
tracks. This avoids a misleading private implant label that does not match the
material participants actually receive.

## Outputs

Normal samples must have all three roles before sealing:

- `read_1`;
- `read_2`; and
- `long_reads`.

GenomePuzzle rejects incomplete triplets unless the exact private
`fault_type` declares the missing role. Every fault requires passing
materialisation evidence.
