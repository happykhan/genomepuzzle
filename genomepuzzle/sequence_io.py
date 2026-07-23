"""Streaming, deterministic sequence-file handling for participant assets."""

from __future__ import annotations

import gzip
import io
import os
import re
import tempfile
from pathlib import Path
from typing import Iterator, TextIO


def _open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        if "w" in mode:
            raw = open(path, "wb")
            compressed = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
            return io.TextIOWrapper(compressed, encoding="utf-8", newline="")
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def fastq_records(path: str | os.PathLike[str]) -> Iterator[tuple[str, str, str, str]]:
    """Yield validated FASTQ records without loading a dataset into memory."""

    source = Path(path)
    with _open_text(source, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                raise ValueError("truncated FASTQ record in {0}".format(source))
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError("invalid FASTQ structure in {0}".format(source))
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise ValueError("FASTQ sequence/quality length mismatch in {0}".format(source))
            yield header, sequence, separator, quality


def _pair_key(header: str) -> str:
    key = header[1:].strip().split()[0]
    key = re.sub(r"([/_-])R?[12]([/_-])", r"\1\2", key, count=1)
    for suffix in ("/1", "/2", "_R1", "_R2", "_1", "_2", "-1", "-2"):
        if key.endswith(suffix):
            return key[: -len(suffix)]
    return key


def validate_paired_fastq(
    r1: str | os.PathLike[str],
    r2: str | os.PathLike[str],
    sample_id: str | None = None,
) -> int:
    """Validate paired FASTQs, optionally requiring a public ID in every header."""

    count = 0
    iterator_r1 = fastq_records(r1)
    iterator_r2 = fastq_records(r2)
    while True:
        record_r1 = next(iterator_r1, None)
        record_r2 = next(iterator_r2, None)
        if record_r1 is None and record_r2 is None:
            break
        if record_r1 is None or record_r2 is None:
            raise ValueError("paired FASTQs contain different record counts")
        if _pair_key(record_r1[0]) != _pair_key(record_r2[0]):
            raise ValueError("paired FASTQ identifiers differ")
        if sample_id and (
            sample_id not in record_r1[0] or sample_id not in record_r2[0]
        ):
            raise ValueError("non-anonymous paired FASTQ header")
        count += 1
    if count == 0:
        raise ValueError("FASTQ pair contains no reads")
    return count


def validate_single_fastq(
    path: str | os.PathLike[str], sample_id: str | None = None
) -> int:
    """Validate one FASTQ, optionally requiring a public ID in every header."""

    count = 0
    for count, (header, _, _, _) in enumerate(fastq_records(path), start=1):
        if sample_id and sample_id not in header:
            raise ValueError("non-anonymous FASTQ header")
    if count == 0:
        raise ValueError("FASTQ contains no reads")
    return count


def anonymize_paired_fastq(
    source_r1: str | os.PathLike[str],
    source_r2: str | os.PathLike[str],
    output_r1: str | os.PathLike[str],
    output_r2: str | os.PathLike[str],
    sample_id: str,
) -> int:
    """Reheader a FASTQ pair while validating pairing and record structure."""

    target_r1 = Path(output_r1)
    target_r2 = Path(output_r2)
    target_r1.parent.mkdir(parents=True, exist_ok=True)
    target_r2.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    iterator_r1 = fastq_records(source_r1)
    iterator_r2 = fastq_records(source_r2)
    with _open_text(target_r1, "wt") as handle_r1, _open_text(target_r2, "wt") as handle_r2:
        while True:
            record_r1 = next(iterator_r1, None)
            record_r2 = next(iterator_r2, None)
            if record_r1 is None and record_r2 is None:
                break
            if record_r1 is None or record_r2 is None:
                raise ValueError("paired FASTQs contain different record counts")
            if _pair_key(record_r1[0]) != _pair_key(record_r2[0]):
                raise ValueError("paired FASTQ identifiers differ")
            count += 1
            identifier = "{0}_read_{1:09d}".format(sample_id, count)
            handle_r1.write("@{0}/1\n{1}+\n{2}".format(identifier, record_r1[1], record_r1[3]))
            handle_r2.write("@{0}/2\n{1}+\n{2}".format(identifier, record_r2[1], record_r2[3]))
    if count == 0:
        raise ValueError("FASTQ pair contains no reads")
    return count


def anonymize_paired_fastq_in_place(
    r1: str | os.PathLike[str], r2: str | os.PathLike[str], sample_id: str
) -> int:
    r1_path = Path(r1)
    r2_path = Path(r2)
    temp_dir = Path(tempfile.mkdtemp(prefix=".reheader-", dir=str(r1_path.parent)))
    try:
        temp_r1 = temp_dir / r1_path.name
        temp_r2 = temp_dir / r2_path.name
        count = anonymize_paired_fastq(r1_path, r2_path, temp_r1, temp_r2, sample_id)
        os.replace(temp_r1, r1_path)
        os.replace(temp_r2, r2_path)
        return count
    finally:
        try:
            temp_dir.rmdir()
        except OSError:
            pass


def anonymize_single_fastq_in_place(
    path: str | os.PathLike[str], sample_id: str, role: str = "long"
) -> int:
    source = Path(path)
    temp = source.parent / ".reheader-{0}".format(source.name)
    count = 0
    try:
        with _open_text(temp, "wt") as handle:
            for count, (_, sequence, _, quality) in enumerate(
                fastq_records(source), start=1
            ):
                handle.write(
                    "@{0}_{1}_{2:09d}\n{3}+\n{4}".format(
                        sample_id, role, count, sequence, quality
                    )
                )
        if count == 0:
            raise ValueError("FASTQ contains no reads")
        os.replace(temp, source)
        return count
    finally:
        if temp.exists():
            temp.unlink()
