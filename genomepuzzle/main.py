"""
This module provides command-line interface for genome sample generation and error introduction.
"""

import argparse
import logging
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from genomepuzzle.create_error import introduce_errors
from genomepuzzle.simulate_reads import simulate_reads
from genomepuzzle.rapid import rapid
from genomepuzzle.contamination import contamination_menu

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def parse_arguments():
    """
    Parses command-line arguments for genome sample generation and error introduction.

    This function sets up an argument parser with two subcommands:
    1. 'simulate': Generates genome samples based on provided parameters.
    2. 'errors': Introduces errors into existing genome samples.

    Returns:
        argparse.Namespace: Parsed command-line arguments.

    Subcommands:
        simulate:
            --num_samples (int): Number of samples to generate (default: 3).
            --samplelist (str): Path to the list of samples to use (default: 'samplelist.csv').
            --species (str): Species to use for generating samples (default: 'K. pneumoniae').
            --output_dir (str): Directory to save the generated samples (default: 'output_dataset').
            --random_seed (int): Random seed for reproducibility (default: 42).

        errors:
            --samplelist (str): Path to the list of samples to use (default: 'output_dataset/sample_sheet.csv').
            --error_proportion (float): Proportion of bad samples to introduce (default: 0.6).
            --random_seed (int): Random seed for reproducibility (default: 42).
            --contamination_list (str): Path to the list of contaminants to use (default: 'samplelist.csv').
            --output_dir (str): Directory to save the generated samples with errors (default: 'output_final').
    """
    parser = argparse.ArgumentParser(description="Generate genome samples.")
    subparsers = parser.add_subparsers(dest="command")

    # Subparser for simulate reads
    simulate_parser = subparsers.add_parser("simulate", help="Simulate reads")
    simulate_parser.add_argument(
        "--num_samples", type=int, help="Number of samples to generate", default=10
    )
    simulate_parser.add_argument(
        "--samplelist",
        type=str,
        help="List of samples to use",
        default="samplelist.csv",
    )
    simulate_parser.add_argument(
        "--species",
        type=str,
        help="Species to use for generating samples",
        default="K. pneumoniae",
    )
    simulate_parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the generated samples",
        default="output_dataset",
    )
    simulate_parser.add_argument(
        "--random_seed", type=int, help="Random seed for reproducibility", default=42
    )
    simulate_parser.set_defaults(
        func=lambda args: simulate_reads(
            args.num_samples,
            args.samplelist,
            args.species,
            args.output_dir,
            args.random_seed,
        )
    )

    # Subparser for introduce errors
    error_parser = subparsers.add_parser("errors", help="Introduce errors in reads")
    error_parser.add_argument(
        "--samplelist",
        type=str,
        help="List of samples to use",
        default="output_dataset/sample_sheet.csv",
    )
    error_parser.add_argument(
        "--error_proportion", type=float, help="Proportion of bad samples", default=0.6
    )
    error_parser.add_argument(
        "--random_seed", type=int, help="Random seed for reproducibility", default=42
    )
    error_parser.add_argument(
        "--contamination_list",
        type=str,
        help="List of contaminants to use",
        default="samplelist.csv",
    )
    error_parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the generated samples",
        default="output_final",
    )
    error_parser.set_defaults(
        func=lambda args: introduce_errors(
            args.samplelist,
            args.error_proportion,
            args.contamination_list,
            args.output_dir,
            args.random_seed,
        )
    )
    # Subparser for rapid test (generate small datasets)
    rapid_parser = subparsers.add_parser("rapid", help="Generate small genome assembly datasets")
    rapid_parser.add_argument(
        "--samplelist",
        type=str,
        help="List of samples to use",
        default="rapid_data.csv",
    )
    rapid_parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the generated samples",
        default="rapid_dataset",
    )
    rapid_parser.set_defaults(
        func=lambda args: rapid(
            args.output_dir,
            args.samplelist,
        )
    )
    contaimination_parser = subparsers.add_parser("contamination", help="Generate contamination samples")
    contaimination_parser.add_argument(
        "--num_samples", type=int, help="Number of samples to generate", default=10
    )
    contaimination_parser.add_argument(
        "--samplelist",
        type=str,
        help="List of samples to use",
        default="samplelist.csv",
    )
    contaimination_parser.add_argument(
        "--species",
        type=str,
        help="Species to use for generating samples",
        default="K. pneumoniae",
    )
    contaimination_parser.add_argument(
        "--type",
        type=str,
        choices=["Species", "ST"],
        help="Contamination type to use for generating samples (Species or ST)",
        default="Species",
    )
    contaimination_parser.add_argument(
        "--random_seed", type=int, help="Random seed for reproducibility", default=42
    )    
    contaimination_parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the generated samples",
        default="contamination_dataset",
    )
    contaimination_parser.add_argument(
        "--assemble",
        action="store_true",
        help="Flag to indicate if the samples should be assembled",
    )
    contaimination_parser.set_defaults(
        func=lambda args: contamination_menu(
            args.num_samples,
            args.samplelist,
            args.species,
            args.type,
            args.output_dir,
            args.assmble,
            args.random_seed,
        )
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    if args.command:
        args.func(args)
    else:
        logging.error("No command given")
        sys.exit(1)
