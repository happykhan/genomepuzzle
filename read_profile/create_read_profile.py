import argparse
import csv
import subprocess
import logging
import os
import json
import gzip

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def fetch_assembly(accessions, output_dir):
    if "ncbi_dataset.zip" not in os.listdir(output_dir):
        accessions = " ".join(accessions)
        command = f"./bin/datasets download genome accession {accessions}"
        subprocess.run(command, shell=True, check=True)
        # Move the downloaded file to the output directory
        subprocess.run(f"mv ncbi_dataset.zip {output_dir}", shell=True, check=True)
        logging.info(f"Downloaded assembly files to {output_dir}")
    else:
        logging.info(f"Assembly files already downloaded in {output_dir}")
    # Unzip the downloaded file
    try:
        subprocess.run(f"unzip -o {os.path.join(output_dir, 'ncbi_dataset.zip')} -d {output_dir}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to unzip assembly files: {e}")
    logging.info(f"Unzipped assembly files in {output_dir}")
    

# take csv of assemblies as input
parser = argparse.ArgumentParser(description='Process some assemblies.')
parser.add_argument('--csv_file', type=str, help='Path to the CSV file of assemblies', default='kleb.csv')
parser.add_argument('--output_dir', type=str, help='Path to the output directory', default='output')
parser.add_argument('--tiles', action='store_true', help='Error with tiles [broken do not use]')
args = parser.parse_args()

# open csv file
samples = []
with open(args.csv_file, mode='r', encoding='utf-8-sig') as file:
    csv_reader = csv.DictReader(file)
    for row in csv_reader:
        if row.get('ASSEMBLY') and row.get('SHORT_READS'):
            samples.append({'id' : row['ID'], 'assembly': row['ASSEMBLY'], 'reads': row['SHORT_READS']})
    
logging.info(f"Found {len(samples)} samples in the CSV file")
# Create output directory if it doesn't exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
    logging.info(f"Created output directory at {args.output_dir}")
else:
    logging.info(f"Output directory already exists at {args.output_dir}")

# fetch assembly file with ncbi datasets e.g. ./datasets download genome accession GCA_902156115
assembly_accesions = [sample['assembly'] for sample in samples]
fetch_assembly(assembly_accesions, args.output_dir)
# for each fasta file in output_dir / ncbi_dataset / data  ; get path, and add path to samples list 
fasta_dir = os.path.join(args.output_dir, 'ncbi_dataset', 'data')
for root, dirs, files in os.walk(fasta_dir):
    for file in files:
        if file.endswith('.fna'):
            fasta_path = os.path.join(root, file)
            for sample in samples:
                if sample['assembly'] in fasta_path:
                    sample['fasta_path'] = fasta_path

for sample in samples:
    read_accession = sample['reads']
    read_r1_path = os.path.join(args.output_dir, f"{sample['id']}_1.fastq.gz")
    read_r2_path = os.path.join(args.output_dir, f"{sample['id']}_2.fastq.gz")
    unnammed_r1 = os.path.join(args.output_dir, f"un.{sample['id']}_1.fastq.gz")
    unnammed_r2 = os.path.join(args.output_dir, f"un.{sample['id']}_2.fastq.gz")    
    reseq_stats_file = os.path.join(args.output_dir, f"{read_accession}.bam.reseq")
    reseq_prob_file = os.path.join(args.output_dir, f"{read_accession}.bam.reseq.ipf")
    sorted_bam = os.path.join(args.output_dir, f"{sample['id']}.bam")
    if not os.path.exists(reseq_stats_file) and not os.path.exists(reseq_prob_file):
        if not os.path.exists(read_r1_path):
            command = f"ffq {read_accession} --ftp" # TODO: Add hanlding to use AWS instead.
            try:
                result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                logging.error(f"Command '{e.cmd}' failed with exit status {e.returncode}. Retrying...")
                try:
                    result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Command '{e.cmd}' failed again with exit status {e.returncode}. Skipping this sample...")
                    continue
            json_output = result.stdout
            # Assuming `json_output` is the JSON string containing the data
            data = json.loads(json_output)
            urls = [entry['url'] for entry in data]
            for url in urls:
                if "_1.fastq.gz" in url:
                    subprocess.run(f"wget -O {unnammed_r1} {url}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                elif "_2.fastq.gz" in url:
                    subprocess.run(f"wget -O {unnammed_r2} {url}", shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            logging.info(f"Downloaded reads for {read_accession}")
            if args.tiles:
                # TODO: This process seems to be broken I am giving the wrong input. 
                command = f"docker run --platform linux/amd64 -v {os.path.abspath(args.output_dir)}:/output quay.io/biocontainers/reseq:1.1--py38h7ce28ed_4 reseq-prepare-names.py /output/un.{sample['id']}_1.fastq.gz /output/un.{sample['id']}_2.fastq.gz"
                logging.info(f"Running command: {command}")
                with gzip.open(read_r1_path, 'wb') as f:
                    subprocess.run(command, shell=True, check=True, stdout=f)
                command = f"docker run --platform linux/amd64 -v {os.path.abspath(args.output_dir)}:/output quay.io/biocontainers/reseq:1.1--py38h7ce28ed_4 reseq-prepare-names.py /output/un.{sample['id']}_2.fastq.gz /output/un.{sample['id']}_1.fastq.gz"
                logging.info(f"Running command: {command}")
                with gzip.open(read_r2_path, 'wb') as f:
                    subprocess.run(command, shell=True, check=True, stdout=f)                
            else:
                # move files to the correct name
                os.rename(unnammed_r1, read_r1_path)
                os.rename(unnammed_r2, read_r2_path)
        else:
            logging.info(f"Reads already downloaded for {read_accession}")
        sample['read_r1'] = read_r1_path
        sample['read_r2'] = read_r2_path
        # use reseq to create read profile
        profile_output_dir = os.path.join(args.output_dir)    
        if os.path.exists(sorted_bam):
            logging.info(f"Reads already mapped to assembly for sample {sample['id']}")
        else:
            # Run minimap2 to map reads to assembly
            command = f"bin/minimap2 -ax sr {sample['fasta_path']} {sample['read_r1']} {sample['read_r2']} | samtools sort -m 10G -@ 4 -T _tmp -o {sorted_bam} -"
            logging.info(f"Running command: {command}")
            subprocess.run(command, shell=True, check=True)
            logging.info(f"Mapped reads to assembly and sorted BAM for sample {sample['id']}")

            # Index BAM file
            command = f"samtools index {sorted_bam}"
            subprocess.run(command, shell=True, check=True)
            logging.info(f"Indexed BAM file for sample {sample['id']}")

        # Create profile 
        if args.tiles:
            command = f"docker run --platform linux/amd64 -v {os.path.abspath(profile_output_dir)}:/{profile_output_dir} quay.io/biocontainers/reseq:1.1--py38h7ce28ed_4 reseq illuminaPE --tiles -j 0 -r /{sample['fasta_path']} -b /{sorted_bam} --stopAfterEstimation"
        else:
            command = f"docker run --platform linux/amd64 -v {os.path.abspath(profile_output_dir)}:/{profile_output_dir} quay.io/biocontainers/reseq:1.1--py38h7ce28ed_4 reseq illuminaPE -j 0 -r /{sample['fasta_path']} -b /{sorted_bam} --stopAfterEstimation"
        logging.info(f"Running command: {command}")
        try:
            subprocess.run(command, shell=True, check=True)
            # rename output files
            command = f"mv {sorted_bam}.reseq {reseq_stats_file}"
            logging.info(f"Running command: {command}")
            subprocess.run(command, shell=True, check=True)
            command = f"mv {sorted_bam}.reseq.ipf {reseq_prob_file}"
            logging.info(f"Running command: {command}")
            subprocess.run(command, shell=True, check=True)
            logging.info(f"Completed read profile estimation for sample {sample['id']}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}")
            continue
    # Check other big files are deleted 
    # Remove the big files
    if os.path.exists(read_r1_path):
        os.remove(read_r1_path)
    if os.path.exists(read_r2_path):
        os.remove(read_r2_path)
    if os.path.exists(sorted_bam):
        os.remove(sorted_bam)
    if os.path.exists(f"{sorted_bam}.bai"):
        os.remove(f"{sorted_bam}.bai")
    if os.path.exists(unnammed_r1):
        os.remove(unnammed_r1)
    if os.path.exists(unnammed_r2):
        os.remove(unnammed_r2)
    else:
        logging.info(f"Read profile files already created for sample {sample['id']}, Skipping...")


# Once completed delete ncbi_dataset.zip and ncbi_dataset folder
subprocess.run(f"rm -rf {os.path.join(args.output_dir, 'ncbi_dataset.zip')}", shell=True, check=True)
subprocess.run(f"rm -rf {os.path.join(args.output_dir, 'ncbi_dataset')}", shell=True, check=True)
subprocess.run(f"rm  {os.path.join(args.output_dir, 'README.md')}", shell=True, check=True)

logging.info("Cleaned up downloaded assembly files")








