import csv
import logging
import platform


def check_input_table(sample_list_file):
    """
    Validates the input sample list file and returns its contents as a list of dictionaries.

    The function performs the following checks:
    1. Ensures the input file contains all required columns: "SAMPLE_NAME", "SHORT_READS", 
       "ASSEMBLY", "SPECIES", "ST", "AMR", and "USE_ORIGINAL_READS".
    2. Ensures each sample in the file has either an "ASSEMBLY" or "READS" record.

    Args:
        sample_list_file (str): Path to the input sample list file in CSV format.

    Returns:
        list: A list of dictionaries, where each dictionary represents a row in the input file.

    Raises:
        ValueError: If the input file is missing required columns or if any sample does not have 
                    an "ASSEMBLY" or "READS" record.

    Logs:
        Info: Logs when the input file has all required columns and valid rows.
        Error: Logs if any sample does not have an "ASSEMBLY" or "READS" record.
    """
    required_columns = {
        "SAMPLE_NAME",
        "SHORT_READS",
        "ASSEMBLY",
        "SPECIES",
        "ST",
        "AMR",
        "USE_ORIGINAL_READS",
    }
    with open(sample_list_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        columns = set(reader.fieldnames)
        if not required_columns.issubset(columns):
            missing_columns = required_columns - columns
            raise ValueError(
                f"Input file is missing required columns: {', '.join(sorted(missing_columns))}"
            )
    logging.info("Input file %s has all required columns.", sample_list_file)
    # import sample list as dict
    all_sample_list = []
    with open(sample_list_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_sample_list.append(row)
    for sample in all_sample_list:
        if not sample.get("ASSEMBLY") and not sample.get("READS"):
            logging.error(
                "Sample %s does not have an assembly or read record.",
                sample["SAMPLE_NAME"],
            )
            raise ValueError(
                f"Sample {sample['SAMPLE_NAME']} does not have an assembly or read record."
            )
    logging.info("Input file %s has valid rows.", sample_list_file)
    return all_sample_list


def check_operating_system():
    """
    Checks the operating system of the current environment.

    This function determines the operating system using the platform module.
    If the operating system is macOS (Darwin), it logs an informational message.
    If the operating system is not macOS, it logs an error message and raises
    an EnvironmentError.

    Returns:
        str: The name of the operating system.

    Raises:
        EnvironmentError: If the operating system is not macOS.
    """
    os_name = platform.system()
    if os_name == "Darwin":
        logging.info("Operating system is macOS.")
    else:
        logging.error(
            "Operating system %s is not supported. Only macOS is supported.", os_name
        )
        raise EnvironmentError(
            f"Operating system {os_name} is not supported. Only macOS is supported."
        )
    return os_name
