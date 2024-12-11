import logging 
import os
import csv 
import subprocess
import random
from genomepuzzle.util import check_input_table
from genomepuzzle.util import check_operating_system

def cleanup_output_dir(output_dir):
    logging.info(f"Cleaning up output directory {output_dir}")
    os.remove(os.path.join(output_dir, 'ncbi_dataset.zip'))
    ncbi_dataset_dir = os.path.join(output_dir, 'ncbi_dataset')
    if os.path.exists(ncbi_dataset_dir):
        subprocess.run(f"rm -rf {ncbi_dataset_dir}", shell=True, check=True)
        logging.info(f"Removed directory {ncbi_dataset_dir}")
    # remove md5sum.txt if it exists
    md5sum_file = os.path.join(output_dir, 'md5sum.txt')
    if os.path.exists(md5sum_file):
        os.remove(md5sum_file)
        logging.info(f"Removed file {md5sum_file}")
    # remove README.md if it exists
    readme_file = os.path.join(output_dir, 'README.md')
    if os.path.exists(readme_file):
        os.remove(readme_file)
        logging.info(f"Removed file {readme_file}")
    
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
    

def run_art(sample, output_dir, reference_genome, output_r1, output_r2):

        art_r1 = os.path.join(output_dir, f"{sample['public_name']}_R1.fq")
        art_r2 = os.path.join(output_dir, f"{sample['public_name']}_R2.fq")
        command = f'bin/art_illumina -ss {sample['platform']} -i {reference_genome} -l {sample['read_length']} -f {sample['coverage']} -o {os.path.join(output_dir, f"{sample['public_name']}_R")} -p -m {sample['fragment_length']} -s {sample['standard_deviation']} --rndSeed {sample['random_seed']} -na'
        # rename the output files to the desired names
        

        subprocess.run(command, shell=True, check=True)
        logging.info(f"Running command: {command}")
        # gzip output_r1
        logging.info(f"gzipping {art_r1}...")
        subprocess.run(f"gzip -f {art_r1}", shell=True, check=True)
        # gzip output_r2
        logging.info(f"gzipping {art_r2}...")
        subprocess.run(f"gzip -f {art_r2}", shell=True, check=True)
        # rename the files
        os.rename(art_r1 + '.gz', output_r1)
        os.rename(art_r2 + '.gz', output_r2)
        
        
def download_reads(sample, output_dir, output_r1, output_r2):
        fastqdump_output_r1 = os.path.join(output_dir, f"{sample['SHORT_READS']}_1.fastq")
        fastqdump_output_r2 = os.path.join(output_dir, f"{sample['SHORT_READS']}_2.fastq")                    
        logging.info(f"Fetching reads for {sample['SAMPLE_NAME']}")
        command = f"bin/fasterq-dump --outdir {output_dir} --split-files {sample['SHORT_READS']}"
        subprocess.run(command, shell=True, check=True)
        logging.info(f"Downloaded reads for {sample['SAMPLE_NAME']}")

        logging.info(f"gzipping {fastqdump_output_r1}...")
        subprocess.run(f"gzip -f {fastqdump_output_r1}", shell=True, check=True)
        logging.info(f"gzipping {fastqdump_output_r2}...")
        subprocess.run(f"gzip -f {fastqdump_output_r2}", shell=True, check=True)
        # rename the files
        os.rename(fastqdump_output_r1 + '.gz', output_r1)
        os.rename(fastqdump_output_r2 + '.gz', output_r2)

def fetch_reads(all_sample_list, output_dir):
    for sample in all_sample_list:
        if sample.get('USE_ORIGINAL_READS').upper() == 'TRUE':
            if sample.get('SHORT_READS'):
                output_r1 = os.path.join(output_dir, f"{sample['SAMPLE_NAME']}_R1.fastq.gz")
                output_r2 = os.path.join(output_dir, f"{sample['SAMPLE_NAME']}_R2.fastq.gz")                
                if  not os.path.exists(output_r1) or not os.path.exists(output_r2):
                    download_reads(sample, output_dir, output_r1, output_r2)
                else:
                    logging.info(f"Reads already downloaded for {sample['SAMPLE_NAME']}")
                sample['r1'] = output_r1
                sample['r2'] = output_r2
                sample['coverage'] = -1 
                sample['read_length'] = -1
                sample['platform'] = 'Unknown'
                sample['fragment_length'] = -1
                sample['standard_deviation'] = -1
                sample['random_seed'] = -1
                sample['QC'] = 'PASSED'
                sample['ERROR'] = 'NONE'
                sample['Notes'] = 'None'             
            else:
                logging.error(f"Sample {sample['SAMPLE_NAME']} does not have a read record.")
                raise ValueError(f"Sample {sample['SAMPLE_NAME']} does not have a read record.")
    return all_sample_list

def simulate_reads(num_samples, sample_list, species, output_dir, random_seed=42):
    check_operating_system()
    # import sample list as dict 
    all_sample_list = check_input_table(sample_list)
    random.seed(random_seed)
    # create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory at {output_dir}")

    # Create a reduced sample list 
    species_sample_list = [x for x in all_sample_list if x.get('SPECIES') == species]
    print(f"Generating {num_samples} samples for {species} in {output_dir}")
    if num_samples < len(species_sample_list):
        species_sample_list = random.sample(species_sample_list, num_samples)
    elif num_samples > len(species_sample_list):
        for i in range(num_samples - len(species_sample_list)):
            species_sample_list.append(random.choice(species_sample_list))    
    
    species_sample_list = fetch_reads(species_sample_list, output_dir)
    # download the assembly genomes for use from the sample list (where sequence reads are not available)
    assembly_list = [x['ASSEMBLY'] for x in species_sample_list if not x.get('r1') ]
    if assembly_list:
        fetch_assembly(assembly_list, output_dir)
        # generate samples
        for sample in species_sample_list:
            sample['public_name'] = sample['SAMPLE_NAME']
            reference_genome = [ os.path.join(output_dir, 'ncbi_dataset/data/', x) for x in os.listdir(os.path.join(output_dir, 'ncbi_dataset/data/')) if x.startswith(sample['ASSEMBLY'])][0]
            reference_genome = [ os.path.join(reference_genome, x) for x in os.listdir(reference_genome) if x.endswith('.fna')][0]
            output_r1 = os.path.join(output_dir, f"{sample['public_name']}_R1.fastq.gz")
            output_r2 = os.path.join(output_dir, f"{sample['public_name']}_R2.fastq.gz")
            
            sample['r1'] = output_r1
            sample['r2'] = output_r2
            # coverage is a random number between 40 and 60
            sample['coverage'] = random.randint(40, 60)
            sample['read_length'] = 150
            sample['platform'] = 'HS25'
            sample['fragment_length'] = 200
            sample['standard_deviation'] = 10
            sample['random_seed'] = 42
            sample['QC'] = 'PASSED'
            sample['ERROR'] = 'NONE'
            sample['Notes'] = 'None' 
            if not os.path.exists(output_r1) or not os.path.exists(output_r2):
                # use ART 
                run_art(sample, output_dir, reference_genome, output_r1, output_r2)
            else:
                logging.info(f"Reads already generated for {sample['SAMPLE_NAME']}")
    # write sample sheet
    sample_sheet = os.path.join(output_dir, 'sample_sheet.csv')
    with open(sample_sheet, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=species_sample_list[0].keys())
        writer.writeheader()
        for sample in species_sample_list:
            sample['r1'] = os.path.basename(sample['r1'])
            sample['r2'] = os.path.basename(sample['r2'])
            writer.writerow(sample)
    # clean up output directory
    cleanup_output_dir(output_dir)
    logging.info(f"Finished generating {num_samples} samples for {species} in {output_dir}")

