# genomepuzzle

`genomepuzzle` generates microbial genomics training datasets for short-read assembly, hybrid assembly, contamination exercises, typing exercises, and outbreak reconstruction.

The codebase is organized around the same core pattern for exercise generation:

1. build a clean base sample or dataset
2. choose an explicit implant plan
3. apply at most one implant per sample
4. emit public files plus private answer/manifest outputs

## Installation

Use Pixi.

```bash
pixi install
pixi run genomepuzzle --help
```

The Pixi environment declares:
- Python runtime
- `typer` and `rich` for the CLI
- `art`
- `seqtk`
- `sra-tools`
- `ncbi-datasets-cli`
- `pigz`

If a tool is already on `PATH`, `genomepuzzle` will use it. If not, it falls back to the bundled `bin/` directory where possible.

For authenticated NCBI downloads, set `NCBI_DATASETS_API_KEY` in the environment. The project no longer ships a hardcoded API key.

## CLI

```bash
pixi run genomepuzzle --help
pixi run genomepuzzle short --help
pixi run genomepuzzle long --help
pixi run genomepuzzle rapid --help
```

Top-level commands remain available for convenience, and grouped subcommands expose the cleaner project structure:
- `short simulate`
- `short errors`
- `short contamination`
- `long hybrid`
- `rapid`

## Main workflows

### Clean short-read dataset

```bash
pixi run genomepuzzle short simulate \
  --num-samples 10 \
  --samplelist datasets/samplelist.csv \
  --species "K. pneumoniae" \
  --output-dir output_dataset
```

### Short-read assembly exercise with implants

```bash
pixi run genomepuzzle short errors \
  --sample-sheet output_dataset/sample_sheet.csv \
  --contamination-list datasets/samplelist.csv \
  --error-proportion 0.6 \
  --output-dir output_final
```

Current short-read implant types:
- `CONTAMINATED`
- `LOW_COVERAGE`
- `POOR_QUALITY`
- `TRUNCATED`
- `CORRUPT`

Each sample receives at most one implant.

### Hybrid assembly exercise

```bash
pixi run genomepuzzle long hybrid \
  --samplelist datasets/rapid_data.csv \
  --mode challenge \
  --output-dir hybrid_dataset
```

Current hybrid implant types:
- `CONTAMINATED`
- `LOW_SHORT_COVERAGE`
- `LOW_LONG_COVERAGE`
- `LONG_READ_QUALITY`

Each sample receives at most one implant.

## Outputs

Modern exercise generators emit:
- `sample_sheet.csv`
- `answer_sheet.csv`
- `implant_manifest.csv`
- `implant_manifest.json`

The public `sample_sheet.csv` hides the answer key. The manifest files carry the implanted truth for internal use.

## CSV Schemas

Input and output table schemas are documented in [docs/csv_schemas.md](docs/csv_schemas.md).

## Project shape

The most important modules are:
- `genomepuzzle/main.py`: Typer/Rich CLI
- `genomepuzzle/short_read.py`: short-read dataset implant planning and application
- `genomepuzzle/hybrid.py`: hybrid dataset implant planning and application
- `genomepuzzle/create_error.py`: FASTQ mutation helpers
- `genomepuzzle/simulate_reads.py`: clean short-read generation
- `genomepuzzle/eqa-test.py`: deprecated legacy umbrella workflow

## Validation

The project is set up to run under Pixi:

```bash
pixi run lint
pixi run test
```

## What still needs work

- add real integration tests that exercise external bioinformatics tools, not just Python-level wiring
- keep shrinking legacy compatibility code such as `eqa-test.py`
- document expected biological assumptions for each implant type
