# Hybrid assembly

Hybrid releases combine ART short reads with Badread long reads. The separate
data modalities allow an assessment to test whether participants detect weak
long-read data or disagreement between evidence sources.

## Supported implants

| Implant | Important parameters | Intended signal |
| --- | --- | --- |
| `LOW_SHORT_COVERAGE` | `read_fraction` | Weak Illumina support |
| `LOW_LONG_COVERAGE` | `read_fraction` | Weak long-read support |
| `LONG_READ_QUALITY` | `min_quality`, `max_quality` | Degraded long-read evidence |
| `CONTAMINATED` | `contaminant_source_id`, `contamination_fraction` | Mixture in both modalities |

Normal long-read quantity defaults to `10x` and can be overridden with
`long_quantity`.

```toml
[[samples]]
source_id = "GCA_000000003.1"
implant = "CONTAMINATED"

[samples.implant_parameters]
contaminant_source_id = "GCA_000000004.1"
contamination_fraction = 0.20
long_quantity = "20x"
```

A contaminated hybrid sample is altered in both its short- and long-read
tracks. This avoids a misleading private implant label that does not match the
material participants actually receive.

## Outputs

Every sample must have all three roles before sealing:

- `read_1`;
- `read_2`; and
- `long_reads`.

GenomePuzzle rejects incomplete triplets and any troublesome sample without
passing implant evidence.
