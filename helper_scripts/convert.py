"""
This module provides functionality to convert a CSV file containing genome data
into a tidy format suitable for further analysis.

The input CSV file is expected to have the following columns:
- \ufeffID: Sample identifier
- SHORT_READS: Short reads information
- ASSEMBLY: Assembly information
- Species: Species information
- ST: Sequence type

The output CSV file will have the following columns:
- SAMPLE_NAME: Sample identifier
- SHORT_READS: Short reads information
- ASSEMBLY: Assembly information
- SPECIES: Species information (with " sensu stricto" removed)
- ST: Sequence type (with '~' removed)
- AMR: Antimicrobial resistance information (set to "Unknown")
- USE_ORIGINAL_READS: Flag indicating whether to use original reads (set to "FALSE")

Usage:
    Run this module to read the input CSV file, process the data, and write the
    output to a new CSV file in a tidy format.
"""
import csv

with open("/Users/nfareed/code/genomepuzzle/kleb.csv", "r", encoding="utf-8") as file:
    reader = csv.DictReader(file)
    lines = list(reader)

# Define the output table structure
output_headers = [
    "SAMPLE_NAME",
    "SHORT_READS",
    "ASSEMBLY",
    "SPECIES",
    "ST",
    "AMR",
    "USE_ORIGINAL_READS",
]

with open(
    "/Users/nfareed/code/genomepuzzle/kleb-tidy.csv", "w", encoding="utf-8", newline=""
) as file:
    writer = csv.writer(file)
    writer.writerow(output_headers)
    for row in lines:
        sample_name = row["\ufeffID"]
        short_reads = row["SHORT_READS"]
        assembly = row["ASSEMBLY"]
        species = row["Species"].replace(" sensu stricto", "")
        st = row["ST"].replace('~', '')
        AMR = "Unknown"  # Set as Unknown as per example
        USE_ORIGINAL_READS = "FALSE"
        writer.writerow(
            [sample_name, short_reads, assembly, species, st, AMR, USE_ORIGINAL_READS]
        )
