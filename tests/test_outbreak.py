import csv
import gzip
import json
from pathlib import Path

from genomepuzzle.outbreak import build_outbreak_release
from genomepuzzle.release import load_release_spec


def _write_fastq(path: Path, prefix: str, count: int = 20) -> None:
    with gzip.open(path, "wt") as handle:
        for index in range(count):
            handle.write("@{0}-{1}\nACGT\n+\nIIII\n".format(prefix, index))


def test_build_outbreak_release_keeps_cluster_truth_private(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    for sample in ("tip1", "tip2", "outsider"):
        _write_fastq(source_dir / "{0}_R1.fastq.gz".format(sample), sample + "-1")
        _write_fastq(source_dir / "{0}_R2.fastq.gz".format(sample), sample + "-2")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "Sample,Cluster,Host,Location,SPECIES\n"
        "tip1,1,Human,Oxford,Klebsiella pneumoniae\n"
        "tip2,1,Human,London,Klebsiella pneumoniae\n",
        encoding="utf-8",
    )
    spec_path = tmp_path / "outbreak.toml"
    spec_path.write_text(
        """
schema_version = "1.0"
release_id = "outbreak-2026"
exercise = "outbreak"
mode = "challenge"

[[samples]]
source_id = "tip1"
public_id = "Sample_normal"

[[samples]]
source_id = "tip2"
public_id = "Sample_low"
implant = "LOW_COVERAGE"
[samples.implant_parameters]
read_fraction = 0.5
source_coverage = 1.0

[[samples]]
source_id = "tip1"
identity_key = "tip1-contaminated"
public_id = "Sample_dirty"
implant = "CONTAMINATED"
[samples.implant_parameters]
contaminant_source_id = "outsider"
contamination_fraction = 0.5
""".strip()
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "release"
    manifests = build_outbreak_release(
        load_release_spec(spec_path), source_dir, metadata, output
    )
    public_text = Path(manifests["public_manifest"]).read_text()
    public = json.loads(public_text)
    private = json.loads(Path(manifests["private_manifest"]).read_text())
    assert "Cluster" not in public_text
    assert "SPECIES" not in public_text
    assert "outsider" not in public_text
    assert public["samples"][0]["metadata"] == {
        "Host": "Human",
        "Location": "Oxford",
    }
    assert private["samples"][1]["expected_answers"]["qc_status"] == "FAIL"
    assert (
        private["samples"][1]["expected_answers"]["failure_reason"]
        == "LOW_COVERAGE"
    )
    assert private["samples"][2]["expected_answers"]["cluster"] == "1"

    with open(output / "public/sample_sheet.csv", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["sample_id"] == "Sample_normal"
    assert "Cluster" not in rows[0]
    assert (output / "COMPLETE.json").is_file()
