from genomepuzzle.hybrid import HYBRID_ERROR_TYPES, build_hybrid_error_plan, stable_public_name


def test_build_hybrid_error_plan_practice():
    plan = build_hybrid_error_plan(6, mode="practice", random_seed=42)
    assert len(plan) == 6
    assert set(plan).issubset(set(HYBRID_ERROR_TYPES))
    assert "CONTAMINATED" in plan
    assert "LOW_SHORT_COVERAGE" in plan
    assert "LOW_LONG_COVERAGE" in plan
    assert "LONG_READ_QUALITY" in plan


def test_build_hybrid_error_plan_none():
    plan = build_hybrid_error_plan(4, mode="none", random_seed=42)
    assert plan == ["NORMAL", "NORMAL", "NORMAL", "NORMAL"]


def test_build_hybrid_error_plan_challenge_distribution():
    plan = build_hybrid_error_plan(10, mode="challenge", random_seed=42)
    assert len(plan) == 10
    assert plan.count("CONTAMINATED") == 2
    assert plan.count("LOW_SHORT_COVERAGE") == 2
    assert plan.count("LOW_LONG_COVERAGE") == 2
    assert plan.count("LONG_READ_QUALITY") == 2
    assert plan.count("NORMAL") == 2


def test_stable_public_name_reproducible():
    assert stable_public_name("GCF_000001.1", 42) == stable_public_name("GCF_000001.1", 42)
    assert stable_public_name("GCF_000001.1", 42) != stable_public_name("GCF_000001.1", 43)
