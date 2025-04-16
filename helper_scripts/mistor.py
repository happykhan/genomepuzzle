#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def logmsg(message):
    sys.stderr.write(f"{os.path.basename(__file__)}: {message}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Annotate a genome quickly using blast and an MLST scheme."
    )
    parser.add_argument("assembly", help="Assembly fasta file")
    parser.add_argument(
        "--mlstdir",
        required=True,
        help="Location of MLST files. All files must have .fna extension.",
    )
    parser.add_argument("--numcpus", type=int, default=1, help="Number of CPUs to use")
    parser.add_argument("--tempdir", help="Location of temporary files (optional)")
    args = parser.parse_args()

    tempdir = args.tempdir or tempfile.mkdtemp(
        prefix=f"{os.path.basename(__file__)}.XXXXXX"
    )
    mlstdir = args.mlstdir
    numcpus = args.numcpus
    asm = args.assembly

    mlst_locus = [
        os.path.join(mlstdir, f) for f in os.listdir(mlstdir) if f.endswith(".fna")
    ]
    num_loci = len(mlst_locus)
    if num_loci == 0:
        sys.exit(f"ERROR: no loci were found in {mlstdir}")
    elif num_loci < 7:
        logmsg(f"WARNING: there are only {num_loci} loci in this scheme.")

    blastdb = os.path.join(tempdir, "assembly.fna")
    shutil.copy(asm, blastdb)
    subprocess.run(["makeblastdb", "-dbtype", "nucl", "-in", blastdb], check=True)

    num_loci_per_thread = (num_loci // numcpus) + 1
    with ThreadPoolExecutor(max_workers=numcpus) as executor:
        futures = []
        for i in range(numcpus):
            thread_locus = mlst_locus[
                i * num_loci_per_thread : (i + 1) * num_loci_per_thread
            ]
            futures.append(
                executor.submit(annotation_worker, blastdb, thread_locus, tempdir)
            )
            logmsg(f"{len(thread_locus)} loci being compared in thread {i + 1}")

        features = []
        for future in futures:
            features.extend(future.result())

    seq_records = SeqIO.to_dict(SeqIO.parse(asm, "fasta"))
    for feature in sorted(features, key=lambda f: f.location.start):
        seq_id = feature.qualifiers["seq_id"][0]
        seq_records[seq_id].features.append(feature)

    logmsg(f"Done adding all {len(features)} features.")
    logmsg("Printing to stdout")
    SeqIO.write(seq_records.values(), sys.stdout, "genbank")


def annotation_worker(blastdb, loci, tempdir):
    features = []
    blast_result = os.path.join(tempdir, f"blast{os.getpid()}.tsv")
    with open(blast_result, "w") as blast_result_fh:
        for locus in loci:
            locusname = os.path.basename(locus).replace(".fna", "")
            result = subprocess.run(
                [
                    "blastn",
                    "-query",
                    locus,
                    "-db",
                    blastdb,
                    "-num_threads",
                    "1",
                    "-outfmt",
                    "6",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            best_hit = sorted(
                result.stdout.splitlines(),
                key=lambda x: float(x.split("\t")[11]),
                reverse=True,
            )[0]
            blast_result_fh.write(best_hit + "\n")

    with open(blast_result) as bls_fh:
        for line in bls_fh:
            fields = line.strip().split("\t")
            (
                allele,
                contig,
                identity,
                aln_length,
                mismatch,
                gap,
                qstart,
                qend,
                sstart,
                send,
                evalue,
                score,
            ) = fields
            locusname = allele

            strand = 1
            feat_start, feat_end = int(sstart), int(send)
            if feat_start > feat_end:
                strand = -1
                feat_start, feat_end = feat_end, feat_start

            feature = SeqFeature(
                location=FeatureLocation(feat_start - 1, feat_end, strand),
                type="gene",
                qualifiers={
                    "source": os.path.basename(__file__),
                    "locus": locusname,
                    "gene": allele,
                    "note": f"Evidence:{os.path.basename(__file__)}",
                    "seq_id": contig,
                },
            )
            features.append(feature)

    return sorted(features, key=lambda f: f.location.start)


if __name__ == "__main__":
    main()
