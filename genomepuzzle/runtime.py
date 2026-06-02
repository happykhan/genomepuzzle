"""
Runtime helpers for external tool resolution.
"""

import os
import shutil
import subprocess
import gzip


def project_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def bundled_bin_dir():
    return os.path.join(project_root(), "bin")


def resolve_tool(name):
    """
    Resolve an executable name from PATH first, then from a configured or bundled bin dir.
    """
    path_candidate = shutil.which(name)
    if path_candidate:
        return path_candidate

    configured_bin = os.environ.get("GENOMEPUZZLE_BIN_DIR")
    if configured_bin:
        candidate = os.path.join(configured_bin, name)
        if os.path.exists(candidate):
            return candidate

    bundled = os.path.join(bundled_bin_dir(), name)
    if os.path.exists(bundled):
        return bundled

    return name


def run_command(args, **kwargs):
    """
    Thin checked subprocess wrapper.
    """
    return subprocess.run(args, check=True, **kwargs)


def open_maybe_gzip(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    if "b" in mode:
        return open(path, mode)
    return open(path, mode, encoding="utf-8")


def require_tool(name):
    candidate = resolve_tool(name)
    if not shutil.which(candidate) and not os.path.exists(candidate):
        raise RuntimeError("Required executable not found: {0}".format(name))
    return candidate


def gzip_file(path, threads=None):
    compressor = resolve_tool("pigz") if threads else resolve_tool("gzip")
    command = [compressor, "-f"]
    if threads:
        command.extend(["-p", str(threads)])
    command.append(path)
    return run_command(command)


def unzip_archive(archive_path, output_dir):
    return run_command(
        ["unzip", "-o", archive_path, "-d", output_dir],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def seqtk_sample(input_fastq, output_fastq, seed, amount):
    seqtk = require_tool("seqtk")
    with open_maybe_gzip(input_fastq, "rb") as input_handle:
        with gzip.open(output_fastq, "wb") as output_handle:
            process = subprocess.Popen(
                [seqtk, "sample", "-s", str(seed), "-", str(amount)],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )
            shutil.copyfileobj(input_handle, process.stdin)
            process.stdin.close()
            shutil.copyfileobj(process.stdout, output_handle)
            process.stdout.close()
            return_code = process.wait()
            if return_code != 0:
                raise subprocess.CalledProcessError(
                    return_code,
                    [seqtk, "sample", "-s", str(seed), "-", str(amount)],
                )


def seqtk_shift_quality(input_fastq, output_fastq, decrement):
    seqtk = require_tool("seqtk")
    with gzip.open(output_fastq, "wb") as output_handle:
        process = subprocess.Popen(
            [seqtk, "seq", "-Q{0}".format(decrement), input_fastq],
            stdout=subprocess.PIPE,
        )
        shutil.copyfileobj(process.stdout, output_handle)
        process.stdout.close()
        return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(
                return_code,
                [seqtk, "seq", "-Q{0}".format(decrement), input_fastq],
            )


def concatenate_files(input_paths, output_path):
    with open(output_path, "wb") as output_handle:
        for input_path in input_paths:
            with open(input_path, "rb") as input_handle:
                shutil.copyfileobj(input_handle, output_handle)
