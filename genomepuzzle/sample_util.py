import gzip
import logging
import os
import shutil
import subprocess

from genomepuzzle.assembly_stats import calculate_assembly_stats
from genomepuzzle.runtime import require_tool, run_command


def _move_output(source, destination):
    if not os.path.exists(source):
        raise RuntimeError("Expected output missing: {0}".format(source))
    shutil.move(source, destination)


def _find_first_tool(candidates):
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    raise RuntimeError(
        "Required Pixi executable not found; expected one of: {0}".format(
            ", ".join(candidates)
        )
    )


def run_spades(r1, r2, output_dir, output_prefix):
    local_r1 = os.path.join(output_dir, os.path.basename(r1))
    local_r2 = os.path.join(output_dir, os.path.basename(r2))
    if os.path.abspath(r1) != os.path.abspath(local_r1):
        shutil.copy2(r1, local_r1)
    if os.path.abspath(r2) != os.path.abspath(local_r2):
        shutil.copy2(r2, local_r2)

    native_spades = _find_first_tool(["spades.py", "spades"])
    command = [
        native_spades,
        "-1",
        local_r1,
        "-2",
        local_r2,
        "-o",
        os.path.join(output_dir, "{0}_spades".format(output_prefix)),
    ]
    run_command(command)

    output_contigs = os.path.join(output_dir, "{0}_spades".format(output_prefix), "scaffolds.fasta")
    record = {}
    if os.path.exists(output_contigs):
        stats = calculate_assembly_stats(output_contigs)
        record["spades_total_bases"] = stats["Total assembly size"]
        record["spades_number_of_contigs"] = stats["Number of contigs"]
        record["spades_n50"] = stats["N50"]
        record["spades_gc_content"] = stats["GC content (%)"]
    else:
        logging.error("no spades output")
    destination = os.path.join(output_dir, "{0}_spades.fasta".format(output_prefix))
    _move_output(output_contigs, destination)
    shutil.rmtree(os.path.join(output_dir, "{0}_spades".format(output_prefix)))
    record["spades_assembly"] = os.path.basename(destination)
    return record


def run_badread(output_prefix, assembly_path, output_dir):
    assembly_file = os.path.basename(assembly_path)
    local_assembly = os.path.join(output_dir, assembly_file)
    shutil.copy2(assembly_path, local_assembly)

    command = [
        require_tool("badread"),
        "simulate",
        "--reference",
        local_assembly,
        "--quantity",
        "20x",
    ]
    output_path = os.path.join(output_dir, "{0}_long.fastq.gz".format(output_prefix))
    with open(output_path, "wb") as raw_output:
        with subprocess.Popen(command, stdout=subprocess.PIPE) as process:
            with gzip.GzipFile(fileobj=raw_output, mode="wb") as gz_output:
                shutil.copyfileobj(process.stdout, gz_output)
            return_code = process.wait()
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, command)
    os.remove(local_assembly)
    return {"long_reads_path": output_path}


def run_flye(output_prefix, genomelen, long_read, output_dir):
    command = [
        require_tool("flye"),
        "--nano-raw",
        long_read,
        "-o",
        os.path.join(output_dir, "{0}_flye".format(output_prefix)),
        "-g",
        str(genomelen),
    ]
    run_command(command)
    output_assembly = os.path.join(output_dir, "{0}_longassembly.fasta".format(output_prefix))
    scaffolds = os.path.join(output_dir, "{0}_flye".format(output_prefix), "scaffolds.fasta")
    _move_output(scaffolds, output_assembly)
    shutil.rmtree(os.path.join(output_dir, "{0}_flye".format(output_prefix)))
    stats = calculate_assembly_stats(output_assembly)
    record = {}
    record["flye_total_bases"] = stats["Total assembly size"]
    record["flye_number_of_contigs"] = stats["Number of contigs"]
    record["flye_n50"] = stats["N50"]
    record["flye_gc_content"] = stats["GC content (%)"]
    record["flye_assembly"] = os.path.basename(output_assembly)
    return record
