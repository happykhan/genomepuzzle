# Installation

GenomePuzzle supports Linux laptops and SLURM clusters. The
biological tools and Python application are installed from the checked-in Pixi
lock file.

## Prerequisites

- Git
- [Pixi](https://pixi.sh/)
- sufficient storage for source assemblies, simulated reads and release bundles
- optionally, working `sbatch`, `squeue` and `sacct` commands for production
  cluster workflows

Clone the repository and install the locked environment:

```bash
git clone https://github.com/happykhan/genomepuzzle.git
cd genomepuzzle
pixi install --locked
pixi run genomepuzzle --help
```

Run the lightweight verification suite on the login node:

```bash
pixi run lint
pixi run test
```

The environment includes ART, Badread, Kleborate, SPAdes, Flye, minimap2,
Mashtree, IQ-TREE and the NCBI data utilities. GenomePuzzle resolves these
tools only from the active Pixi environment.

## Private identity salt

Challenge releases normally derive anonymous public IDs using an HMAC salt.
Generate and export a private value before planning a release:

```bash
export GENOMEPUZZLE_ID_SALT="$(openssl rand -hex 32)"
```

Keep the salt in your secret manager. Reusing the same release specification
with a different salt changes its public sample IDs. Practice specifications
may instead freeze explicit `public_id` values.

!!! danger

    Do not commit a production identity salt, source-to-sample mapping or
    generated `private/` directory.

## Documentation

Preview this manual locally:

```bash
pixi run docs-serve
```

Build it with strict link and configuration checks:

```bash
pixi run docs-build
```
