# Genotyping

Typing releases provide anonymous FASTA assemblies. GenomePuzzle runs pinned
Kleborate against the final files and stores raw output and normalised answers
privately.

## Source layout

```text
sources/typing/
└── <source_id>.fasta
```

`.fna`, `.fa` and extensionless FASTA files are also recognised.

## Supported implants

| Implant | Important parameters | Intended signal |
| --- | --- | --- |
| `NORMAL` or `NONE` | — | Clean assembly |
| `FRAGMENTED` | `fragment_size` | Categorically extreme fragmentation |
| `MIXED_CONTIGS` | `contaminant_source_id`, `contamination_fraction` | A 30–50% different-species mixture |
| `ZERO_BYTE_ASSEMBLY` | — | A literal zero-byte FASTA |
| `WRONG_ORGANISM` | `replacement_source_id` | Complete replacement with another organism |

```toml
[[samples]]
source_id = "GCA_000000099.1"
implant = "MIXED_CONTIGS"

[samples.implant_parameters]
contaminant_source_id = "GCA_000000100.1"
contamination_fraction = 0.45
```

The contaminant assembly must exist in the same source directory. Source names
are removed from public headers.

## Participant answers

The standard template includes species, sequence type, capsule and O-antigen
typing, `wzi`, and carbapenemase calls. `st` is the canonical name; internal
analyser labels such as `kleborate_st` are not exposed.

Before publishing, inspect difficult cases to ensure the expected Kleborate
result represents a useful proficiency question rather than an accidental
generator failure.
