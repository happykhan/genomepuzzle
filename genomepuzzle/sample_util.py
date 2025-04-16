
import os 
import logging 
from genomepuzzle.assembly_stats import calculate_assembly_stats

@staticmethod
def run_spades(r1, r2, output_dir, output_prefix):
    abs_output_dir = os.path.abspath(output_dir)
    # copy r1 to the output dir if needed
    if not r1.startswith(output_dir):
        os.system(f"cp {r1} {output_dir}")
    # copy r2 to the output dir if needed
    if not r2.startswith(output_dir):
        os.system(f"cp {r2} {output_dir}")
    os.system(
        f"docker run --platform linux/x86_64 -v {abs_output_dir}:/data "
        f"quay.io/biocontainers/spades:3.11.0--py36_0 spades.py "
        f"-1  /data/{os.path.basename(r1)} "
        f"-2 /data/{os.path.basename(r2)} "
        f"-o /data/{output_prefix}_spades"
    )
    output_contigs = f"{output_dir}/{output_prefix}_spades/scaffolds.fasta"
    record = {}
    if os.path.exists(output_contigs):
        stats = calculate_assembly_stats(output_contigs)
        record["spades_total_bases"] = stats["Total assembly size"]
        record["spades_number_of_contigs"] = stats["Number of contigs"]
        record["spades_n50"] = stats["N50"]
        record["spades_gc_content"] = stats["GC content (%)"]
    else:
        logging.error("no spades output")
    # copy the assembly to the output dir
    os.system(f"mv {output_contigs} {output_dir}/{output_prefix}_spades.fasta")
    # delete spades output dir
    os.system(f"rm -r {output_dir}/{output_prefix}_spades")
    record["spades_assembly"] = f"{output_prefix}_scaffolds.fasta"
    return record


def run_badread(output_prefix, assembly_path, output_dir):
    abs_output_dir = os.path.abspath(output_dir)
    assembly_file = os.path.basename(assembly_path)
    # copy the assembly to the output dir
    os.system(f"cp {assembly_path} {output_dir}")
    os.system(
        f"docker run --platform linux/x86_64 -v {abs_output_dir}:/data "
        f"quay.io/biocontainers/badread:0.4.1--pyhdfd78af_0 badread simulate "
        f"--reference /data/{assembly_file} --quantity 20x | gzip > {abs_output_dir}/{output_prefix}_long.fastq.gz"
    )
    record = {}
    record['long_reads_path'] = f"{output_dir}/{output_prefix}_long.fastq.gz"
    os.remove(f"{output_dir}/{assembly_file}")
    return record

def run_flye(output_prefix, genomelen, long_read, output_dir):
    abs_output_dir = os.path.abspath(output_dir)
    long_read_file = os.path.basename(long_read)
    os.system(
        f"docker run --platform linux/x86_64 -v {abs_output_dir}:/data "
        f"quay.io/biocontainers/flye:2.3.6--py27ha92aebf_3 flye --nano-raw "
        f"/data/{long_read_file} -o /data/{output_prefix}_flye -g {genomelen}"
    )
    # copy the assembly to the output dir
    output_assembly = f"{output_dir}/{output_prefix}_longassembly.fasta"
    os.system(f"mv {output_dir}/{output_prefix}_flye/scaffolds.fasta {output_assembly}")
    os.system(f"rm -r {output_dir}/{output_prefix}_flye")
    stats = calculate_assembly_stats(output_assembly)
    record = {}
    record["flye_total_bases"] = stats["Total assembly size"]
    record["flye_number_of_contigs"] = stats["Number of contigs"]
    record["flye_n50"] = stats["N50"]
    record["flye_gc_content"] = stats["GC content (%)"]
    record["flye_assembly"] = os.path.basename(output_assembly)
    return record
