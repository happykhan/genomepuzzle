import gzip
import logging
import os
import shutil
import subprocess

from genomepuzzle.assembly_stats import calculate_assembly_stats
from genomepuzzle.runtime import require_tool, resolve_tool, run_command


DOCKER_PLATFORM = "linux/x86_64"


def _docker_base(abs_output_dir, image):
    require_tool("docker")
    return [
        "docker",
        "run",
        "--rm",
        "--platform",
        DOCKER_PLATFORM,
        "-v",
        "{0}:/data".format(abs_output_dir),
        image,
    ]


def _move_output(source, destination):
    if not os.path.exists(source):
        raise RuntimeError("Expected output missing: {0}".format(source))
    shutil.move(source, destination)


def _find_first_tool(candidates):
    for candidate in candidates:
        resolved = resolve_tool(candidate)
        if shutil.which(resolved) or os.path.exists(resolved):
            return resolved
    return None


def run_spades(r1, r2, output_dir, output_prefix):
    abs_output_dir = os.path.abspath(output_dir)
    local_r1 = os.path.join(output_dir, os.path.basename(r1))
    local_r2 = os.path.join(output_dir, os.path.basename(r2))
    if os.path.abspath(r1) != os.path.abspath(local_r1):
        shutil.copy2(r1, local_r1)
    if os.path.abspath(r2) != os.path.abspath(local_r2):
        shutil.copy2(r2, local_r2)

    native_spades = _find_first_tool(["spades.py", "spades"])
    if native_spades:
        command = [
            native_spades,
            "-1",
            local_r1,
            "-2",
            local_r2,
            "-o",
            os.path.join(output_dir, "{0}_spades".format(output_prefix)),
        ]
    else:
        command = _docker_base(abs_output_dir, "quay.io/biocontainers/spades:3.11.0--py36_0")
        command.extend(
            [
                "spades.py",
                "-1",
                "/data/{0}".format(os.path.basename(local_r1)),
                "-2",
                "/data/{0}".format(os.path.basename(local_r2)),
                "-o",
                "/data/{0}_spades".format(output_prefix),
            ]
        )
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
    abs_output_dir = os.path.abspath(output_dir)
    assembly_file = os.path.basename(assembly_path)
    local_assembly = os.path.join(output_dir, assembly_file)
    shutil.copy2(assembly_path, local_assembly)

    native_badread = _find_first_tool(["badread"])
    if native_badread:
        command = [
            native_badread,
            "simulate",
            "--reference",
            local_assembly,
            "--quantity",
            "20x",
        ]
    else:
        command = _docker_base(abs_output_dir, "quay.io/biocontainers/badread:0.4.1--pyhdfd78af_0")
        command.extend(
            [
                "badread",
                "simulate",
                "--reference",
                "/data/{0}".format(assembly_file),
                "--quantity",
                "20x",
            ]
        )
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
    abs_output_dir = os.path.abspath(output_dir)
    long_read_file = os.path.basename(long_read)
    native_flye = _find_first_tool(["flye"])
    if native_flye:
        command = [
            native_flye,
            "--nano-raw",
            long_read,
            "-o",
            os.path.join(output_dir, "{0}_flye".format(output_prefix)),
            "-g",
            str(genomelen),
        ]
    else:
        command = _docker_base(abs_output_dir, "quay.io/biocontainers/flye:2.3.6--py27ha92aebf_3")
        command.extend(
            [
                "flye",
                "--nano-raw",
                "/data/{0}".format(long_read_file),
                "-o",
                "/data/{0}_flye".format(output_prefix),
                "-g",
                str(genomelen),
            ]
        )
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
