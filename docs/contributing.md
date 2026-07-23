# Contributing

GenomePuzzle treats reproducibility, privacy and auditable failure as product
features.

## Development checks

```bash
pixi install --locked
pixi run lint
pixi run test
pixi run docs-build
```

The login node may run unit tests, specification validation, packaging and
small smoke tests. Read simulation, assembly, Kleborate batches, phylogenetic
inference, dataset-wide QC and other heavy biological work must use the
repository's SLURM workflow.

## Change expectations

- Add tests for contract, privacy or resume behaviour when those areas change.
- Keep public and private data models explicit.
- Preserve deterministic seeds and row-order-independent identity.
- Pin new runtime software in Pixi; do not add binaries to the repository.
- Update this manual when CLI behaviour, supported implants or the release
  contract changes.
- Do not add a silent local fallback for failed SLURM submission.

## Documentation

Documentation uses MkDocs Material. Add pages under `docs/`, include them in
`mkdocs.yml`, then run the strict build. Pushes to `main` publish the generated
site through GitHub Pages; pull requests build it without deploying.

The repository's design record remains in
[GenomePuzzle Assessment Dataset Plan](assessment-platform-plan.md).
