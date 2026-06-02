from genomepuzzle.short_read import (
    LOW_COVERAGE_FRACTION,
    POOR_QUALITY_MAX,
    POOR_QUALITY_MIN,
    SHORT_READ_ERROR_TYPES,
    TRUNCATED_READ_LENGTH,
    build_short_read_error_plan,
    make_sample_context,
)


def test_build_short_read_error_plan_respects_one_implant_per_sample():
    plan = build_short_read_error_plan(10, error_proportion=0.5, random_seed=7)
    assert len(plan) == 10
    assert sum(1 for item in plan if item != "NORMAL") == 5
    assert all(item in SHORT_READ_ERROR_TYPES for item in plan)


def test_build_short_read_error_plan_zero_errors():
    plan = build_short_read_error_plan(4, error_proportion=0, random_seed=11)
    assert plan == ["NORMAL", "NORMAL", "NORMAL", "NORMAL"]


def test_make_sample_context_assigns_public_names():
    record = {
        "SPECIES": "Klebsiella pneumoniae",
        "r1": "/tmp/input_R1.fastq.gz",
        "r2": "/tmp/input_R2.fastq.gz",
    }
    context = make_sample_context(record, index=2, random_seed=42, source_dir=".")
    assert context.public_name == "sample03"
    assert context.implant.error_type == "NORMAL"


def test_short_read_implant_constants_are_obvious():
    assert LOW_COVERAGE_FRACTION == 0.15
    assert POOR_QUALITY_MIN == 5
    assert POOR_QUALITY_MAX == 14
    assert TRUNCATED_READ_LENGTH == 35
