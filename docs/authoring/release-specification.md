# Release specifications

A versioned TOML file is the source of truth for one practice or challenge
dataset. Keep specifications under `releases/<round>/` and review them like
code.

## Top-level fields

| Field | Required | Meaning |
| --- | --- | --- |
| `schema_version` | Yes | Specification schema; currently `"1.0"` |
| `release_id` | Yes | Stable identifier containing letters, numbers, `.`, `_` or `-` |
| `exercise` | Yes | `typing`, `assembly`, `hybrid` or `outbreak` |
| `mode` | Yes | `practice` or `challenge` |
| `master_seed` | Recommended | Non-negative seed; defaults to `42` |
| `id_salt_env` | No | Environment variable holding the private salt |
| `title` | No | Participant-facing title |
| `description` | No | Participant-facing summary |
| `instructions` | No | Ordered participant instructions |
| `pass_threshold` | No | Score required to pass; defaults to `0.8` |

`release_id` identifies the whole bundle. `source_id` identifies private input
material. `sample_id` is the anonymous public identity derived from the
release, an identity key and a private salt—or supplied explicitly for a
practice dataset.

## Inputs

Paths in `[inputs]` are resolved relative to the specification file.

```toml
[inputs]
source_dir = "../../sources/typing"
```

Typing, assembly and hybrid releases use `source_dir`. Native outbreak
simulation uses `base_genome` plus `metadata_csv`; a pre-simulated outbreak
uses `source_dir` plus `metadata_csv`.

## Samples

Each `[[samples]]` table declares exactly one private source:

```toml
[[samples]]
source_id = "GCA_000000002.1"
identity_key = "stable-purpose-specific-key"
implant = "LOW_COVERAGE"

[samples.implant_parameters]
read_fraction = 0.15

[samples.expected_answers]
species = "Klebsiella pneumoniae"
```

`identity_key` defaults to `source_id`. It makes identity and random-seed
derivation stable when input order changes. It must be unique within the
release.

`public_id` may freeze an identifier for a practice dataset:

```toml
public_id = "Sample_A002"
```

Avoid explicit IDs for blind challenges; use a private salt instead.

Expected answers in the specification are useful when the correct
interpretation is known in advance. An exercise analyser may add or normalise
fields from the final participant-facing material. Expected results must
describe those final files, not pristine sources.

## Reproducibility rules

- Freeze source assets before release generation.
- Keep `release_id`, `identity_key` and `master_seed` stable when resuming.
- Assign at most one implant to each sample.
- Record implant parameters explicitly rather than relying on changing
  defaults.
- Never reorder or edit a submitted plan in place; update the specification and
  create a deliberate new release version when its scientific content changes.
- Treat `pixi.lock` and the generator Git commit as part of release provenance.

See the ready-to-edit files under
[`releases/examples`](https://github.com/happykhan/genomepuzzle/tree/main/releases/examples).
