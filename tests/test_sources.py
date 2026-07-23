import json
import zipfile
from pathlib import Path

from genomepuzzle.release import load_release_spec
from genomepuzzle.sources import (
    fetch_assembly_sources,
    required_assembly_accessions,
)


def _write_spec(path: Path) -> None:
    path.write_text(
        """
schema_version = "1.0"
release_id = "typing-source-test"
exercise = "typing"
mode = "practice"
master_seed = 1

[[samples]]
source_id = "GCA_000001"
public_id = "Sample_one"

[[samples]]
source_id = "GCA_000002"
public_id = "Sample_two"
implant = "MIXED_CONTIGS"
[samples.implant_parameters]
contaminant_source_id = "GCA_000003"
contamination_fraction = 0.1
""".strip()
        + "\n",
        encoding="utf-8",
    )


def test_required_accessions_include_contaminants(tmp_path):
    spec_path = tmp_path / "release.toml"
    _write_spec(spec_path)
    assert required_assembly_accessions(load_release_spec(spec_path)) == (
        "GCA_000001",
        "GCA_000002",
        "GCA_000003",
    )


def test_fetch_sources_stages_stable_files_and_reuses_cache(tmp_path):
    spec_path = tmp_path / "release.toml"
    _write_spec(spec_path)
    calls = []

    def fake_runner(command, check):
        calls.append(command)
        package = Path(command[command.index("--filename") + 1])
        with zipfile.ZipFile(package, "w") as archive:
            for accession, base in (
                ("GCA_000001.1", "A"),
                ("GCA_000002.2", "C"),
                ("GCA_000003.1", "G"),
            ):
                archive.writestr(
                    "ncbi_dataset/data/{0}/{0}_genomic.fna".format(accession),
                    ">source\n" + base * 100 + "\n",
                )

    destination = tmp_path / "sources"
    manifest_path = fetch_assembly_sources(
        load_release_spec(spec_path),
        destination,
        runner=fake_runner,
    )
    fetch_assembly_sources(
        load_release_spec(spec_path),
        destination,
        runner=fake_runner,
    )

    assert len(calls) == 1
    assert (destination / "GCA_000002.fasta").read_text().startswith(">source\nC")
    manifest = json.loads(manifest_path.read_text())
    assert [row["source_id"] for row in manifest["sources"]] == [
        "GCA_000001",
        "GCA_000002",
        "GCA_000003",
    ]
    assert all(len(row["sha256"]) == 64 for row in manifest["sources"])
