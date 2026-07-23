# Build a first release

This walkthrough uses a typing practice release. The same planning, submission
and validation commands apply to every exercise.

## 1. Prepare source assemblies

Create a private source directory containing one FASTA per `source_id`:

```text
sources/typing/
├── GCA_000000001.1.fasta
├── GCA_000000002.1.fasta
└── GCA_000000003.1.fasta
```

Sequence headers may contain source information because GenomePuzzle rewrites
participant-facing headers. The source files themselves must never be
published.

## 2. Copy and edit an example

```bash
mkdir -p releases/my-first-round
cp releases/examples/typing-practice.toml \
  releases/my-first-round/typing-practice.toml
```

Update `release_id`, `[inputs].source_dir`, the samples and any implants. Check
the specification without running biological tools:

```bash
pixi run genomepuzzle release validate-spec \
  --spec releases/my-first-round/typing-practice.toml
```

## 3. Review the plan

```bash
pixi run genomepuzzle release plan \
  --spec releases/my-first-round/typing-practice.toml \
  --output-dir generated/my-first-round/typing-practice
```

Review:

- `build/plan.json` for resolved identities, commands and resources;
- `build/scripts/*.sbatch` for the exact cluster jobs; and
- `build/stages/*.json` for initial stage state.

Planning is lightweight and does not submit work.

## 4. Submit generation

The convenience command creates the plan and submits the workflow:

```bash
pixi run genomepuzzle release build \
  --spec releases/my-first-round/typing-practice.toml \
  --output-dir generated/my-first-round/typing-practice
```

Generation runs first. Independent validation is submitted with an `afterok`
dependency and therefore cannot seal a failed generation.

```bash
pixi run genomepuzzle release status \
  --plan generated/my-first-round/typing-practice/build/plan.json

pixi run genomepuzzle release logs \
  --plan generated/my-first-round/typing-practice/build/plan.json
```

## 5. Verify the sealed release

```bash
pixi run genomepuzzle release validate \
  --release-dir generated/my-first-round/typing-practice \
  --require-complete

pixi run genomepuzzle release inspect \
  --release-dir generated/my-first-round/typing-practice
```

A publishable release has `COMPLETE.json` and reports `status: complete`.
Inspect the private answer key and validation evidence before handing the
bundle to a publisher.

## If a job fails

Fix the source or specification, then use the recorded plan:

```bash
pixi run genomepuzzle release resume \
  --plan generated/my-first-round/typing-practice/build/plan.json
```

Resume preserves previous job IDs, attempts, scripts and logs. It removes only
incomplete release artifacts before retrying the failed stage.
