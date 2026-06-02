"""
Typer-based command-line interface for genomepuzzle.
"""

import logging
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genomepuzzle.contamination import contamination_menu
from genomepuzzle.create_error import introduce_errors
from genomepuzzle.hybrid import create_hybrid_dataset
from genomepuzzle.long_qc import (
    compare_qc_to_manifest,
    summarize_hybrid_dataset,
    write_qc_outputs,
    write_report_outputs,
)
from genomepuzzle.rapid import rapid
from genomepuzzle.simulate_reads import simulate_reads
from genomepuzzle.slurm import build_hybrid_sbatch_script, submit_sbatch_script

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

try:
    import typer
    from rich.console import Console
    from rich.table import Table
except ImportError:  # pragma: no cover - exercised outside the pixi env
    typer = None
    Console = None
    Table = None


console = Console() if Console else None


def _require_cli_dependencies():
    if typer is None:
        raise RuntimeError(
            "Typer/Rich are not installed. Create the managed environment with "
            "`pixi install` and run commands via `pixi run genomepuzzle ...`."
        )


def _print_run_summary(command_name, rows):
    if not console or not Table:
        return
    table = Table(title="genomepuzzle")
    table.add_column("Command", style="bold cyan")
    table.add_column("Value", overflow="fold")
    for key, value in rows:
        table.add_row(key, str(value))
    console.print(table)


if typer is not None:
    app = typer.Typer(
        help="Generate microbial genomics training datasets with explicit implants.",
        no_args_is_help=True,
        add_completion=False,
        rich_markup_mode="rich",
    )
    short_app = typer.Typer(
        rich_markup_mode="rich",
        help="Short-read dataset generation and implant workflows.",
    )
    long_app = typer.Typer(
        rich_markup_mode="rich",
        help="Long-read (ONT/PacBio) QC and hybrid assembly workflows.",
    )
    legacy_app = typer.Typer(
        rich_markup_mode="rich",
        help="Legacy compatibility commands that should be phased out.",
    )
    app.add_typer(short_app, name="short")
    app.add_typer(long_app, name="long")
    app.add_typer(legacy_app, name="legacy")

    @app.command("simulate")
    @short_app.command("simulate")
    def simulate_command(
        num_samples: int = typer.Option(10, help="Number of samples to generate."),
        samplelist: str = typer.Option(
            "samplelist.csv", help="CSV file describing candidate source samples."
        ),
        species: str = typer.Option(
            "K. pneumoniae", help="Species to sample from the input table."
        ),
        output_dir: str = typer.Option(
            "output_dataset", help="Directory to write the generated reads."
        ),
        random_seed: int = typer.Option(42, help="Random seed for reproducibility."),
    ):
        _print_run_summary(
            "simulate",
            [
                ("samples", num_samples),
                ("samplelist", samplelist),
                ("species", species),
                ("output", output_dir),
                ("seed", random_seed),
            ],
        )
        simulate_reads(num_samples, samplelist, species, output_dir, random_seed)

    @app.command("errors")
    @short_app.command("errors")
    def errors_command(
        sample_sheet: str = typer.Option(
            "output_dataset/sample_sheet.csv",
            "--sample-sheet",
            help="Clean short-read sample sheet to mutate.",
        ),
        error_proportion: float = typer.Option(
            0.6, help="Fraction of samples to receive exactly one implant."
        ),
        contamination_list: str = typer.Option(
            "samplelist.csv",
            help="CSV file to draw contamination assemblies from.",
        ),
        output_dir: str = typer.Option(
            "output_final", help="Directory to write the implanted dataset."
        ),
        random_seed: int = typer.Option(42, help="Random seed for reproducibility."),
    ):
        _print_run_summary(
            "errors",
            [
                ("sample_sheet", sample_sheet),
                ("error_proportion", error_proportion),
                ("contamination_list", contamination_list),
                ("output", output_dir),
                ("seed", random_seed),
            ],
        )
        introduce_errors(
            sample_sheet, error_proportion, contamination_list, output_dir, random_seed
        )

    @app.command("rapid")
    def rapid_command(
        samplelist: str = typer.Option(
            "rapid_data.csv", help="CSV file describing assembly accessions."
        ),
        output_dir: str = typer.Option(
            "rapid_dataset", help="Directory to write the generated dataset."
        ),
    ):
        _print_run_summary(
            "rapid", [("samplelist", samplelist), ("output", output_dir)]
        )
        rapid(output_dir, samplelist)

    @app.command("hybrid")
    @long_app.command("hybrid")
    def hybrid_command(
        samplelist: str = typer.Option(
            "rapid_data.csv", help="CSV file describing source assemblies."
        ),
        contamination_list: str = typer.Option(
            None, help="Optional CSV file of assemblies to use as contaminants."
        ),
        mode: str = typer.Option(
            "challenge",
            help="Preset implant profile: practice, challenge, or none.",
        ),
        output_dir: str = typer.Option(
            "hybrid_dataset", help="Directory to write the generated hybrid dataset."
        ),
        random_seed: int = typer.Option(42, help="Random seed for reproducibility."),
    ):
        _print_run_summary(
            "hybrid",
            [
                ("samplelist", samplelist),
                ("contamination_list", contamination_list or "<self>"),
                ("mode", mode),
                ("output", output_dir),
                ("seed", random_seed),
            ],
        )
        create_hybrid_dataset(
            output_dir=output_dir,
            samplelist=samplelist,
            contamination_list=contamination_list,
            mode=mode,
            random_seed=random_seed,
        )

    @long_app.command("qc")
    def long_qc_command(
        sample_sheet: str = typer.Option(
            ...,
            "--sample-sheet",
            help="Hybrid sample sheet to summarize.",
        ),
        dataset_dir: str = typer.Option(
            None,
            help="Directory containing the FASTQ files. Defaults to the sample sheet directory.",
        ),
        output_csv: str = typer.Option(
            "long_qc_summary.csv",
            help="Path to write the QC summary CSV.",
        ),
        output_json: str = typer.Option(
            None,
            help="Optional path to also write the QC summary as JSON.",
        ),
    ):
        _print_run_summary(
            "long qc",
            [
                ("sample_sheet", sample_sheet),
                ("dataset_dir", dataset_dir or os.path.dirname(sample_sheet) or "."),
                ("output_csv", output_csv),
                ("output_json", output_json or "<none>"),
            ],
        )
        sample_rows = summarize_hybrid_dataset(sample_sheet, dataset_dir=dataset_dir)
        write_qc_outputs(sample_rows, output_csv=output_csv, output_json=output_json)

    @long_app.command("report")
    def long_report_command(
        sample_sheet: str = typer.Option(..., "--sample-sheet", help="Hybrid sample sheet."),
        manifest: str = typer.Option(
            ...,
            "--manifest",
            help="Hybrid implant_manifest.csv to compare against.",
        ),
        dataset_dir: str = typer.Option(
            None,
            help="Directory containing FASTQ files. Defaults to the sample sheet directory.",
        ),
        output_csv: str = typer.Option(
            "hybrid_qc_report.csv",
            help="Path to write the implant-vs-QC report CSV.",
        ),
        output_json: str = typer.Option(
            None,
            help="Optional path to also write the report as JSON.",
        ),
    ):
        _print_run_summary(
            "long report",
            [
                ("sample_sheet", sample_sheet),
                ("manifest", manifest),
                ("dataset_dir", dataset_dir or os.path.dirname(sample_sheet) or "."),
                ("output_csv", output_csv),
                ("output_json", output_json or "<none>"),
            ],
        )
        qc_rows = summarize_hybrid_dataset(sample_sheet, dataset_dir=dataset_dir)
        report_rows = compare_qc_to_manifest(qc_rows, manifest)
        write_report_outputs(report_rows, output_csv=output_csv, output_json=output_json)

    @long_app.command("hybrid-slurm")
    def long_hybrid_slurm_command(
        samplelist: str = typer.Option(
            "datasets/rapid_data.csv", help="CSV file describing source assemblies."
        ),
        output_dir: str = typer.Option(
            "hybrid_dataset", help="Directory to write the generated hybrid dataset."
        ),
        mode: str = typer.Option("challenge", help="Preset implant profile."),
        contamination_list: str = typer.Option(
            None, help="Optional CSV file of assemblies to use as contaminants."
        ),
        random_seed: int = typer.Option(42, help="Random seed for reproducibility."),
        partition: str = typer.Option("short", help="Slurm partition."),
        cpus_per_task: int = typer.Option(8, help="Slurm CPUs per task."),
        mem_gb: int = typer.Option(32, help="Slurm memory in GB."),
        time_limit: str = typer.Option("12:00:00", help="Slurm time limit."),
        script_path: str = typer.Option(
            None, help="Optional path to write the sbatch script."
        ),
        submit: bool = typer.Option(
            True, help="Submit the job immediately after writing the script."
        ),
    ):
        repo_dir = os.getcwd()
        output_dir_abs = os.path.abspath(output_dir)
        script_path = script_path or os.path.join(output_dir_abs, "run_hybrid.sbatch")
        os.makedirs(output_dir_abs, exist_ok=True)
        script_text = build_hybrid_sbatch_script(
            samplelist=os.path.abspath(samplelist),
            output_dir=output_dir_abs,
            mode=mode,
            contamination_list=os.path.abspath(contamination_list)
            if contamination_list
            else None,
            random_seed=random_seed,
            partition=partition,
            cpus_per_task=cpus_per_task,
            mem_gb=mem_gb,
            time_limit=time_limit,
            repo_dir=repo_dir,
        )
        with open(script_path, "w", encoding="utf-8") as handle:
            handle.write(script_text)
        _print_run_summary(
            "long hybrid-slurm",
            [
                ("samplelist", samplelist),
                ("output_dir", output_dir_abs),
                ("mode", mode),
                ("partition", partition),
                ("cpus", cpus_per_task),
                ("mem_gb", mem_gb),
                ("time_limit", time_limit),
                ("script_path", script_path),
            ],
        )
        if submit:
            job_id = submit_sbatch_script(script_path)
            if console:
                console.print("Submitted Slurm job: [bold]{0}[/bold]".format(job_id))

    @app.command("contamination")
    @short_app.command("contamination")
    def contamination_command(
        num_samples: int = typer.Option(10, help="Number of samples to generate."),
        samplelist: str = typer.Option(
            "samplelist.csv", help="CSV file describing source samples."
        ),
        species: str = typer.Option(
            "K. pneumoniae", help="Primary species to contaminate."
        ),
        contamination_type: str = typer.Option(
            "Species",
            "--type",
            help="Contamination mode: Species or ST.",
        ),
        output_dir: str = typer.Option(
            "contamination_dataset", help="Directory to write the dataset."
        ),
        assemble: bool = typer.Option(
            False, help="Also assemble the contaminated reads."
        ),
        random_seed: int = typer.Option(42, help="Random seed for reproducibility."),
    ):
        _print_run_summary(
            "contamination",
            [
                ("samples", num_samples),
                ("samplelist", samplelist),
                ("species", species),
                ("type", contamination_type),
                ("output", output_dir),
                ("assemble", assemble),
                ("seed", random_seed),
            ],
        )
        contamination_menu(
            num_samples,
            samplelist,
            species,
            contamination_type,
            output_dir,
            assemble,
            random_seed,
        )

    @legacy_app.command("eqa")
    def legacy_eqa_note():
        """
        Legacy umbrella dataset generation is deprecated.
        """
        raise typer.BadParameter(
            "The old eqa-test.py entrypoint is being retired. Use "
            "`short simulate`, `short errors`, `long hybrid`, `rapid`, and "
            "`short contamination` directly."
        )

else:
    app = None


def main():
    _require_cli_dependencies()
    app()


if __name__ == "__main__":
    main()
