from genomepuzzle.slurm import build_hybrid_sbatch_script


def test_build_hybrid_sbatch_script():
    script = build_hybrid_sbatch_script(
        samplelist="/tmp/rapid.csv",
        output_dir="/tmp/hybrid-output",
        mode="practice",
        contamination_list="/tmp/contaminants.csv",
        random_seed=7,
        partition="short",
        cpus_per_task=4,
        mem_gb=16,
        time_limit="02:00:00",
        repo_dir="/repo",
    )
    assert "#SBATCH -p short" in script
    assert "#SBATCH -c 4" in script
    assert "--samplelist" in script
    assert "--mode practice" in script
    assert "--contamination-list" in script
    assert "pixi" in script
