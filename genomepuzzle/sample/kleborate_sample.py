from genomepuzzle.sample.basic_sample import BasicSample
import os 
import logging
import subprocess
import csv

class KleborateSample(BasicSample):

    def __init__(
        self,
        input_dict: dict,
        output_dir: str,
        random_seed: int = 42,
        bin_dir: str = "bin",
        work_dir: str = "work",        
    ):
        super().__init__(input_dict, output_dir, random_seed=random_seed, bin_dir=bin_dir, work_dir=work_dir)
        self.analysis_results = None

    def to_dict(self):
        dic = super().to_dict()
        dic["analysis_results"] = self.analysis_results
        return dic

    def to_short_dict(self, public_name=False):
        if public_name:
            sample_name = self.public_name
            fasta = f"{self.public_name}.fasta"
        else:
            sample_name = self.sample_name   
            fasta = os.path.basename(self.assembly_file)
        short_dic = {
            'sample': sample_name,
            'st': '',
            'k_locus': '',
            'capsule_type': '',
            'wzi': '',
            'o_locus': '',
            'o_type': '',
            'bla_carb': '',
            'fasta': fasta,
        }
        return short_dic

    def run_analysis(self):
        self.fetch_assembly()
        kleboutput = os.path.join(self.work_dir, f"{self.sample_name}_kleborate_output", "klebsiella_pneumo_complex_output.txt")
        if not os.path.exists(kleboutput):
            logging.info("Running kleborate for %s", self.sample_name)
            # Run kleborate with docker
            command = (
                f"docker run --rm -v {os.path.abspath(self.work_dir)}:/data "
                f"quay.io/biocontainers/kleborate:3.1.3--pyhdfd78af_0 kleborate -a /data/{os.path.basename(self.assembly_file)} "
                f"-o /data/{self.sample_name}_kleborate_output -p kpsc"
            )
            subprocess.run(command, shell=True, check=True)
            logging.info("Kleborate analysis complete. Output saved to %s", kleboutput)
        self.analysis_results = kleboutput

    def get_results(self):
        kleboutput = os.path.join(self.work_dir, f"{self.sample_name}_kleborate_output", "klebsiella_pneumo_complex_output.txt")
        if not self.assembly_file:
            self.fetch_assembly()
        if os.path.exists(kleboutput):
            self.analysis_results = kleboutput
        if not self.analysis_results:
            self.run_analysis()
        # Parse kleborate output
        kleborate_results = {}
        with open(self.analysis_results, "r") as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                kleborate_results = row
                break  # Assuming we only need the first row
        return {
            'kleborate_st': kleborate_results['klebsiella_pneumo_complex__mlst__ST'],
            'k_locus': kleborate_results['klebsiella_pneumo_complex__kaptive__K_locus'],
            'capsule_type': kleborate_results['klebsiella_pneumo_complex__kaptive__K_type'],
            'wzi': kleborate_results['klebsiella_pneumo_complex__wzi__wzi'],
            'o_locus': kleborate_results['klebsiella_pneumo_complex__kaptive__O_locus'],
            'o_type': kleborate_results['klebsiella_pneumo_complex__kaptive__O_type'],
            'bla_carb': kleborate_results['klebsiella_pneumo_complex__amr__Bla_Carb_acquired'],
            'fasta': os.path.basename(self.assembly_file), 
        }
