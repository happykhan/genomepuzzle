import gzip
import json

from genomepuzzle.outbreak_generation import generate_outbreak_release
from genomepuzzle.release import load_release_spec


def _fake_short(reference, r1, r2, **_):
    public_id = reference.read_text(encoding="utf-8").splitlines()[0][1:]
    for path, mate in ((r1, 1), (r2, 2)):
        with gzip.open(path, "wt") as handle:
            for index in range(20):
                handle.write(
                    "@{0}-read-{2}/{1}\nACGT\n+\nIIII\n".format(
                        public_id, mate, index
                    )
                )


def test_native_outbreak_generation_records_cluster_mutations(tmp_path, monkeypatch):
    reference = tmp_path / "reference.fasta"
    reference.write_text(">reference\n" + "ACGT" * 1000 + "\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "Sample,Cluster,SPECIES,Host\n"
        "tip-1,A,Klebsiella pneumoniae,Human\n"
        "tip-2,A,Klebsiella pneumoniae,Human\n",
        encoding="utf-8",
    )
    spec_path = tmp_path / "outbreak.toml"
    spec_path.write_text(
        """
release_id = "native-outbreak"
exercise = "outbreak"
mode = "practice"
master_seed = 9

[[samples]]
source_id = "tip-1"
public_id = "Sample_1"

[[samples]]
source_id = "tip-2"
public_id = "Sample_2"
implant = "LOW_COVERAGE"
[samples.implant_parameters]
read_fraction = 0.5
""".strip()
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "genomepuzzle.outbreak_generation._simulate_short_reads", _fake_short
    )

    output = tmp_path / "release"
    generate_outbreak_release(
        load_release_spec(spec_path),
        reference,
        metadata,
        output,
        shared_cluster_snps=10,
        private_snps=2,
    )
    simulation = json.loads((output / "build/outbreak_simulation.json").read_text())

    assert (output / "COMPLETE.json").is_file()
    assert simulation["samples"][0]["shared_mutations"] == 10
    assert simulation["samples"][0]["private_mutations"] == 2
    assert simulation["samples"][1]["total_mutations"] == 12
