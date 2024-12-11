import csv 
import logging
import platform

def check_input_table(sample_list_file):
    required_columns = {'SAMPLE_NAME', 'SHORT_READS', 'ASSEMBLY', 'SPECIES', 'ST', 'AMR' , 'USE_ORIGINAL_READS'}
    with open(sample_list_file, 'r') as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames)
        if not required_columns.issubset(columns):
            missing_columns = required_columns - columns
            raise ValueError(f"Input file is missing required columns: {', '.join(sorted(missing_columns))}")
    logging.info(f"Input file {sample_list_file} has all required columns.")
    # import sample list as dict 
    all_sample_list = []
    with open(sample_list_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_sample_list.append(row)
    for sample in all_sample_list:
        if not sample.get('ASSEMBLY') and not sample.get('READS'):
            logging.error(f"Sample {sample['SAMPLE_NAME']} does not have an assembly or read record.")
            raise ValueError(f"Sample {sample['SAMPLE_NAME']} does not have an assembly or read record.")
    logging.info(f"Input file {sample_list_file} has valid rows.")
    return all_sample_list


def check_operating_system():
    os_name = platform.system()
    if os_name == 'Darwin':
        logging.info("Operating system is macOS.")
    else:
        logging.error(f"Operating system {os_name} is not supported. Only macOS is supported.")
        raise EnvironmentError(f"Operating system {os_name} is not supported. Only macOS is supported.")
    return os_name
