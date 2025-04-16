import os
import subprocess


def create_config_file(output_dir, config_file_name, example_tree, base_genome):
    # pick a branch in the tree and rename to ref
    # copy tree file to output dir
    tree_output_path = os.path.join(output_dir, os.path.basename(example_tree))
    with open(example_tree, "r", encoding="utf-8") as tree_file:
        tree_content = tree_file.read()
    with open(tree_output_path, "w", encoding="utf-8") as tree_output_file:
        tree_output_file.write(tree_content)

    # copy base genome to output dir
    base_genome_output_path = os.path.join(output_dir, os.path.basename(base_genome))
    with open(base_genome, "r", encoding="utf-8") as genome_file:
        genome_content = genome_file.read()
    with open(base_genome_output_path, "w", encoding="utf-8") as genome_output_file:
        genome_output_file.write(genome_content)

    # calculate genome length
    genome_length = 0
    with open(base_genome, "r", encoding="utf-8") as genome_file:
        for line in genome_file:
            if not line.startswith(">"):
                genome_length += sum(1 for char in line if char in "ATCG")
    print(f"Genome length: {genome_length}")
    invariate_sites = max(genome_length // 50, 10000)
    print(f"Number of invariate sites: {invariate_sites}")
    config_content = f"""
#REQUIRED PARAMETERS
treefile_path = /data/{os.path.basename(example_tree)} #Must be newick or Nexus format, and include branch lengths
number_of_variable_sites = {invariate_sites} # The number of variable sites to simulate
base_genome_name = Outgroup #Should be the label of a tip in your tree
base_genome_path = /data/{os.path.basename(base_genome)}
output_dir = /data/simulated

#parameters of evolutionary model (comma separated), in order ac, ag, at, cg, ct, gc (gc = 1)
rate_matrix = 1,1,1,1,1,1

#parameters for read simulation
coverage = 20 #either an integer or a file name of a comma delimited file with tip names and coverage

#Optional evolutionary model parameters
gamma_shape = 5 #default is no rate variation across sites

#parameters for clustering of variable site locations (OPTIONAL)
mutation_clustering = ON
percent_clustered = 0.25 #The percentage of variable sites whose distance to another site is drawn from the clustering distribution
exponential_mean = 125 #Minimum allowed value = 2
    """
    # create config file in output dir
    config_path = os.path.join(output_dir, config_file_name)
    with open(config_path, "w", encoding="utf-8") as config_file:
        config_file.write(config_content)
    return config_file_name

def run_docker(config_file_name , output_dir):
    current_dir = os.path.abspath(output_dir)
    docker_image = "snacktavish/treetoreads"
    config_file = f"/data/{config_file_name}"
    
    command = [
        "docker", "run",
        "-v", f"{current_dir}:/data",
        docker_image,
        config_file
    ]
    print("Running docker command")
    print(" ".join(command))
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    
    if result.returncode == 0:
        print("Docker run successful")
        print(result.stdout)
    else:
        print("Docker run failed")
        print(result.stderr)

def run_mashtree(output_dir):
    fasta_files_dir = os.path.join(output_dir, 'simulated', "fasta_files")
    os.makedirs(fasta_files_dir, exist_ok=True)
    
    file_of_files_path = os.path.join(fasta_files_dir, "file_of_files.txt")
    with open(file_of_files_path, "w", encoding="utf-8") as file_of_files:
        for fasta_file in os.listdir(fasta_files_dir):
            if fasta_file.endswith(".fasta"):
                file_of_files.write(f"/data/simulated/fasta_files/{fasta_file}\n")

    command = [
        "docker", "run",
        "-v", f"{os.path.abspath(output_dir)}:/data",
        "staphb/mashtree",
        "mashtree",
        "--mindepth", "0",
        "--numcpus", "4",
        "--outtree", "/data/mashtree.nwk",
        "--file-of-files", "/data/simulated/fasta_files/file_of_files.txt"
    ]
    print("Running mashtree command")
    print(" ".join(command))
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    
    if result.returncode == 0:
        print("Mashtree run successful")
        print(result.stdout)
    else:
        print("Mashtree run failed")
        print(result.stderr)



if __name__ == "__main__":
    output_dir = "small_run"    
    os.makedirs(output_dir, exist_ok=True)
    config_path = create_config_file(output_dir, "small_run.cfg", "small_run/small_tree.nwk", "small_run/GCF_000292505.1_longassembly.fasta")
    run_docker(config_path, output_dir)
    run_mashtree(output_dir)