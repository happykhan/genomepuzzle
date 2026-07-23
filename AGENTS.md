# Repository execution policy

All computationally heavy biological data generation and analysis must run
through SLURM. This includes read simulation, short-read and hybrid assembly,
Kleborate batches, phylogenetic inference, outbreak generation, contamination
or error implantation over full datasets, and dataset-wide QC.

The login/local process may perform only lightweight work: release-spec
validation, source indexing, manifest and checksum creation, packaging,
unit tests, small smoke tests, and deployment preparation.

Use the repository's SLURM helpers and generated `sbatch` scripts. Record the
SLURM job ID, command/tool versions, resource request, seed, and output path in
the private release provenance. Do not silently fall back to running a heavy
job locally when SLURM submission fails.
