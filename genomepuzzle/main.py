"""
Typer-based command-line interface for genomepuzzle.
"""

import json
import logging
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genomepuzzle.contamination import contamination_menu
from genomepuzzle.combined_pack import build_combined_pack, validate_combined_pack
from genomepuzzle.contract import inspect_release, validate_release_bundle
from genomepuzzle.create_error import introduce_errors
from genomepuzzle.hybrid import create_hybrid_dataset
from genomepuzzle.long_qc import (
    compare_qc_to_manifest,
    summarize_hybrid_dataset,
    write_qc_outputs,
    write_report_outputs,
)
from genomepuzzle.outbreak import build_outbreak_release
from genomepuzzle.outbreak_generation import generate_outbreak_release
from genomepuzzle.rapid import rapid
from genomepuzzle.read_generation import generate_read_release
from genomepuzzle.reads_release import load_expected_answers, package_read_release
from genomepuzzle.release import load_release_spec, resolve_release_samples
from genomepuzzle.simulate_reads import simulate_reads
from genomepuzzle.slurm import build_hybrid_sbatch_script, submit_sbatch_script
from genomepuzzle.sources import fetch_assembly_sources
from genomepuzzle.typing import build_typing_release, run_kleborate
from genomepuzzle.workflow import (
    build_release_plan,
    run_stage,
    submit_plan,
    workflow_status,
    write_release_plan,
)

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
    long_app = typer.Typer(
        rich_markup_mode="rich",
        help="Read-only QC and reporting for long-read datasets.",
    )
    legacy_app = typer.Typer(
        rich_markup_mode="rich",
        help="Legacy compatibility commands that should be phased out.",
    )
    release_app = typer.Typer(
        rich_markup_mode="rich",
        help="Plan, build, resume, validate and inspect assessment releases.",
    )
    app.add_typer(long_app, name="long")
    app.add_typer(release_app, name="release")
    app.add_typer(legacy_app, name="legacy")

    @release_app.command("validate-spec")
    def validate_release_spec_command(
        spec: str = typer.Option(
            ...,
            "--spec",
            help="Versioned TOML release specification.",
        ),
        id_salt: str = typer.Option(
            None,
            "--id-salt",
            help=(
                "Private ID salt. If omitted, the environment variable named "
                "by id_salt_env in the specification is used."
            ),
        ),
        output_json: str = typer.Option(
            None,
            "--output-json",
            help="Optional private path to write resolved source-to-public mappings.",
        ),
    ):
        release_spec = load_release_spec(spec)
        try:
            resolved = resolve_release_samples(release_spec, id_salt=id_salt)
        except ValueError as exc:
            raise typer.BadParameter(str(exc))
        implant_counts = {}
        for sample in resolved:
            implant_counts[sample.implant] = implant_counts.get(sample.implant, 0) + 1
        _print_run_summary(
            "release validate-spec",
            [
                ("release_id", release_spec.release_id),
                ("exercise", release_spec.exercise),
                ("mode", release_spec.mode),
                ("samples", len(resolved)),
                (
                    "implants",
                    ", ".join(
                        "{0}={1}".format(key, implant_counts[key])
                        for key in sorted(implant_counts)
                    ),
                ),
            ],
        )
        if output_json:
            output_path = os.path.abspath(output_json)
            output_parent = os.path.dirname(output_path)
            if output_parent:
                os.makedirs(output_parent, exist_ok=True)
            with open(output_path, "w", encoding="utf-8") as handle:
                json.dump(
                    {
                        "schema_version": release_spec.schema_version,
                        "release_id": release_spec.release_id,
                        "exercise": release_spec.exercise,
                        "mode": release_spec.mode,
                        "samples": [
                            {
                                "source_id": sample.source_id,
                                "identity_key": sample.identity_key,
                                "sample_id": sample.sample_id,
                                "random_seed": sample.random_seed,
                                "implant": sample.implant,
                                "implant_parameters": dict(
                                    sample.implant_parameters
                                ),
                            }
                            for sample in resolved
                        ],
                    },
                    handle,
                    indent=2,
                    sort_keys=True,
                )
                handle.write("\n")
            if console:
                console.print(
                    "Wrote private resolved mapping: [bold]{0}[/bold]".format(
                        output_path
                    )
                )

    @release_app.command("plan")
    def plan_release_command(
        spec: str = typer.Option(..., "--spec", help="Versioned release TOML."),
        output_dir: str = typer.Option(
            ..., "--output-dir", help="Release directory and workflow state."
        ),
        id_salt: str = typer.Option(
            None,
            "--id-salt",
            help="Private ID salt; otherwise use the specification environment variable.",
        ),
        partition: str = typer.Option("short", "--partition", help="SLURM partition."),
    ):
        plan = build_release_plan(
            spec,
            output_dir,
            id_salt=id_salt,
            partition=partition,
        )
        plan_path = write_release_plan(plan)
        _print_run_summary(
            "release plan",
            [
                ("release_id", plan.release_id),
                ("samples", len(plan.samples)),
                ("stages", ", ".join(stage.name for stage in plan.stages)),
                ("plan", plan_path),
            ],
        )

    @release_app.command("build")
    def build_release_command(
        spec: str = typer.Option(..., "--spec", help="Versioned release TOML."),
        output_dir: str = typer.Option(..., "--output-dir", help="Release directory."),
        id_salt: str = typer.Option(
            None,
            "--id-salt",
            help="Private ID salt; otherwise use the specification environment variable.",
        ),
        partition: str = typer.Option("short", "--partition", help="SLURM partition."),
    ):
        plan = build_release_plan(
            spec,
            output_dir,
            id_salt=id_salt,
            partition=partition,
        )
        plan_path = write_release_plan(plan)
        jobs = submit_plan(plan_path)
        _print_run_summary(
            "release build",
            [
                ("release_id", plan.release_id),
                ("plan", plan_path),
                ("submitted", ", ".join("{0}={1}".format(*item) for item in jobs.items())),
            ],
        )

    @release_app.command("submit")
    def submit_release_command(
        plan: str = typer.Option(..., "--plan", help="Workflow build/plan.json."),
    ):
        jobs = submit_plan(plan)
        if not jobs:
            raise typer.BadParameter("no stages are eligible for submission")
        _print_run_summary("release submit", list(jobs.items()))

    @release_app.command("resume")
    def resume_release_command(
        plan: str = typer.Option(..., "--plan", help="Workflow build/plan.json."),
    ):
        jobs = submit_plan(plan, retry=True)
        if not jobs:
            raise typer.BadParameter("no failed or planned stages are eligible")
        _print_run_summary("release resume", list(jobs.items()))

    @release_app.command("status")
    def release_status_command(
        plan: str = typer.Option(..., "--plan", help="Workflow build/plan.json."),
    ):
        states = workflow_status(plan)
        _print_run_summary(
            "release status",
            [
                (
                    state["stage"],
                    "{0}{1}".format(
                        state["status"],
                        " ({0})".format(state["job_id"]) if state.get("job_id") else "",
                    ),
                )
                for state in states
            ],
        )

    @release_app.command("logs")
    def release_logs_command(
        plan: str = typer.Option(..., "--plan", help="Workflow build/plan.json."),
    ):
        log_dir = os.path.join(os.path.dirname(os.path.abspath(plan)), "logs")
        files = (
            sorted(os.path.join(log_dir, name) for name in os.listdir(log_dir))
            if os.path.isdir(log_dir)
            else []
        )
        _print_run_summary(
            "release logs",
            [("log_dir", log_dir), ("files", "\n".join(files) or "<none>")],
        )

    @release_app.command("run-stage", hidden=True)
    def run_release_stage_command(
        plan: str = typer.Option(..., "--plan"),
        stage: str = typer.Option(..., "--stage"),
    ):
        run_stage(plan, stage)

    @release_app.command("validate")
    def validate_release_command(
        release_dir: str = typer.Option(..., "--release-dir"),
        require_complete: bool = typer.Option(False, "--require-complete"),
    ):
        report = validate_release_bundle(
            release_dir, require_complete=require_complete
        )
        _print_run_summary(
            "release validate",
            [
                ("release_id", report["release_id"]),
                ("status", report["status"]),
                ("samples", report["sample_count"]),
                ("participant_files", report["participant_file_count"]),
            ],
        )

    @release_app.command("inspect")
    def inspect_release_command(
        release_dir: str = typer.Option(..., "--release-dir"),
    ):
        details = inspect_release(release_dir)
        _print_run_summary("release inspect", list(details.items()))

    @release_app.command("build-combined-pack")
    def build_combined_pack_command(
        spec: str = typer.Option(..., "--spec", help="Combined pack TOML."),
        short_read_release: str = typer.Option(
            ...,
            "--short-read-release",
            help="Completed short-read assembly release directory.",
        ),
        long_read_release: str = typer.Option(
            ...,
            "--long-read-release",
            help="Completed hybrid assembly release directory.",
        ),
        output_dir: str = typer.Option(
            ..., "--output-dir", help="New combined pack directory."
        ),
    ):
        output = build_combined_pack(
            spec,
            short_read_release=short_read_release,
            long_read_release=long_read_release,
            output_dir=output_dir,
        )
        report = validate_combined_pack(output)
        _print_run_summary(
            "release build-combined-pack",
            [
                ("pack_id", report["pack_id"]),
                ("status", report["status"]),
                ("short-read samples", report["samples"]["short-read"]),
                ("long-read samples", report["samples"]["long-read"]),
                ("output", output),
            ],
        )

    @release_app.command("validate-combined-pack")
    def validate_combined_pack_command(
        pack_dir: str = typer.Option(..., "--pack-dir"),
    ):
        report = validate_combined_pack(pack_dir)
        _print_run_summary(
            "release validate-combined-pack",
            [
                ("pack_id", report["pack_id"]),
                ("status", report["status"]),
                ("short-read samples", report["samples"]["short-read"]),
                ("long-read samples", report["samples"]["long-read"]),
            ],
        )

    @release_app.command("build-typing")
    def build_typing_release_command(
        spec: str = typer.Option(..., "--spec", help="Typing release TOML."),
        source_dir: str = typer.Option(
            ..., "--source-dir", help="Directory containing source FASTA files."
        ),
        output_dir: str = typer.Option(
            ..., "--output-dir", help="New directory for the release package."
        ),
        id_salt: str = typer.Option(
            None,
            "--id-salt",
            help="Private ID salt; otherwise use the specification environment variable.",
        ),
        skip_analysis: bool = typer.Option(
            False,
            "--skip-analysis",
            help="Prepare files with pending answers instead of running Kleborate.",
        ),
        kleborate_executable: str = typer.Option(
            "kleborate",
            "--kleborate-executable",
            help="Pinned Kleborate executable or wrapper.",
        ),
    ):
        release_spec = load_release_spec(spec)
        analyser = None
        if not skip_analysis:
            analyser = lambda path: run_kleborate(path, kleborate_executable)
        manifests = build_typing_release(
            release_spec,
            source_dir,
            output_dir,
            id_salt=id_salt,
            analyser=analyser,
        )
        _print_run_summary(
            "release build-typing",
            [
                ("release_id", release_spec.release_id),
                ("samples", len(release_spec.samples)),
                ("output", os.path.abspath(output_dir)),
                ("analysis", "pending" if skip_analysis else "Kleborate"),
                ("public_manifest", manifests["public_manifest"]),
            ],
        )

    @release_app.command("fetch-assemblies")
    def fetch_release_assemblies_command(
        spec: str = typer.Option(
            ..., "--spec", help="Release TOML containing NCBI assembly accessions."
        ),
        output_dir: str = typer.Option(
            ..., "--output-dir", help="Assembly source cache directory."
        ),
        datasets_executable: str = typer.Option(
            "datasets",
            "--datasets-executable",
            help="NCBI datasets executable from the managed environment.",
        ),
        refresh: bool = typer.Option(
            False,
            "--refresh",
            help="Download and atomically replace already-cached source FASTAs.",
        ),
    ):
        release_spec = load_release_spec(spec)
        manifest = fetch_assembly_sources(
            release_spec,
            output_dir,
            datasets_executable=datasets_executable,
            refresh=refresh,
        )
        _print_run_summary(
            "release fetch-assemblies",
            [
                ("release_id", release_spec.release_id),
                ("sources", len(release_spec.samples)),
                ("output", os.path.abspath(output_dir)),
                ("source_manifest", manifest),
            ],
        )

    @release_app.command("build-outbreak")
    def build_outbreak_release_command(
        spec: str = typer.Option(..., "--spec", help="Outbreak release TOML."),
        source_dir: str = typer.Option(
            ..., "--source-dir", help="Directory containing TreeToReads FASTQ files."
        ),
        metadata_csv: str = typer.Option(
            ...,
            "--metadata",
            help="Frozen outbreak metadata and private cluster truth CSV.",
        ),
        output_dir: str = typer.Option(
            ..., "--output-dir", help="New directory for the release package."
        ),
        id_salt: str = typer.Option(
            None,
            "--id-salt",
            help="Private ID salt; otherwise use the specification environment variable.",
        ),
    ):
        release_spec = load_release_spec(spec)
        manifests = build_outbreak_release(
            release_spec,
            source_dir,
            metadata_csv,
            output_dir,
            id_salt=id_salt,
        )
        _print_run_summary(
            "release build-outbreak",
            [
                ("release_id", release_spec.release_id),
                ("samples", len(release_spec.samples)),
                ("output", os.path.abspath(output_dir)),
                ("public_manifest", manifests["public_manifest"]),
            ],
        )

    @release_app.command("generate-outbreak")
    def generate_outbreak_release_command(
        spec: str = typer.Option(..., "--spec"),
        metadata_csv: str = typer.Option(..., "--metadata"),
        output_dir: str = typer.Option(..., "--output-dir"),
        base_genome: str = typer.Option(None, "--base-genome"),
        fault_genome: str = typer.Option(None, "--fault-genome"),
        source_dir: str = typer.Option(None, "--source-dir"),
        id_salt: str = typer.Option(None, "--id-salt"),
    ):
        release_spec = load_release_spec(spec)
        if base_genome:
            manifests = generate_outbreak_release(
                release_spec,
                base_genome,
                metadata_csv,
                output_dir,
                id_salt=id_salt,
                fault_genome=fault_genome,
            )
        elif source_dir:
            manifests = build_outbreak_release(
                release_spec,
                source_dir,
                metadata_csv,
                output_dir,
                id_salt=id_salt,
            )
        else:
            raise typer.BadParameter("--base-genome or --source-dir is required")
        _print_run_summary(
            "release generate-outbreak",
            [
                ("release_id", release_spec.release_id),
                ("samples", len(release_spec.samples)),
                ("public_manifest", manifests["public_manifest"]),
            ],
        )

    @release_app.command("package-reads")
    def package_read_release_command(
        spec: str = typer.Option(
            ..., "--spec", help="Assembly or hybrid release TOML."
        ),
        source_dir: str = typer.Option(
            ...,
            "--source-dir",
            help="Directory containing final implanted FASTQ files.",
        ),
        expected_answers: str = typer.Option(
            None,
            "--expected-answers",
            help=(
                "Private JSON answer object keyed by source ID. Optional when every "
                "[[samples]] has a [samples.expected_answers] table."
            ),
        ),
        implant_validations: str = typer.Option(
            None,
            "--implant-validations",
            help="Private JSON validation object keyed by source ID.",
        ),
        output_dir: str = typer.Option(
            ..., "--output-dir", help="New directory for the release package."
        ),
        id_salt: str = typer.Option(
            None,
            "--id-salt",
            help="Private ID salt; otherwise use the specification environment variable.",
        ),
        preanonymized: bool = typer.Option(
            False,
            "--preanonymized",
            help="Validate and copy FASTQs whose headers already contain public sample IDs.",
        ),
    ):
        release_spec = load_release_spec(spec)
        manifests = package_read_release(
            release_spec,
            source_dir,
            load_expected_answers(expected_answers) if expected_answers else None,
            output_dir,
            id_salt=id_salt,
            implant_validations=(
                load_expected_answers(implant_validations)
                if implant_validations
                else None
            ),
            preanonymized=preanonymized,
        )
        _print_run_summary(
            "release package-reads",
            [
                ("release_id", release_spec.release_id),
                ("exercise", release_spec.exercise),
                ("samples", len(release_spec.samples)),
                ("output", os.path.abspath(output_dir)),
                ("public_manifest", manifests["public_manifest"]),
            ],
        )

    @release_app.command("generate-reads")
    def generate_read_release_command(
        spec: str = typer.Option(..., "--spec"),
        source_dir: str = typer.Option(..., "--source-dir"),
        output_dir: str = typer.Option(..., "--output-dir"),
        id_salt: str = typer.Option(None, "--id-salt"),
    ):
        release_spec = load_release_spec(spec)
        manifests = generate_read_release(
            release_spec,
            source_dir,
            output_dir,
            id_salt=id_salt,
        )
        _print_run_summary(
            "release generate-reads",
            [
                ("release_id", release_spec.release_id),
                ("exercise", release_spec.exercise),
                ("samples", len(release_spec.samples)),
                ("public_manifest", manifests["public_manifest"]),
            ],
        )

    @legacy_app.command("simulate")
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

    @legacy_app.command("errors")
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

    @legacy_app.command("rapid")
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

    @legacy_app.command("hybrid")
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

    @legacy_app.command("hybrid-slurm")
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

    @legacy_app.command("contamination")
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
            "`genomepuzzle release build --spec ... --output-dir ...`."
        )

else:
    app = None


def main():
    _require_cli_dependencies()
    app()


if __name__ == "__main__":
    main()
