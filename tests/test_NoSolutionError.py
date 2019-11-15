import pytest
from dnachisel import (
    DnaOptimizationProblem,
    AvoidChanges,
    AvoidPattern,
    EnforceGCContent,
    NoSolutionError,
)


def test_no_solution_error_frozen_region():
    problem = DnaOptimizationProblem(
        sequence="AAAAATCGTCTCTTTT",
        constraints=[AvoidChanges(), AvoidPattern("BsmBI_site")],
        logger=None,
    )
    with pytest.raises(NoSolutionError) as err:
        problem.resolve_constraints()
    assert "region that cannot be mutated" in str(err.value)


def test_no_solution_error_random_search():
    problem = DnaOptimizationProblem(
        sequence="TTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        constraints=[
            AvoidChanges(location=(0, 10)),
            EnforceGCContent(mini=0.8),
        ],
        logger=None,
    )
    with pytest.raises(NoSolutionError) as err:
        problem.resolve_constraints()
    assert "Random search did not" in str(err.value)


def test_no_solution_error_exhaustive_search():
    problem = DnaOptimizationProblem(
        sequence="TTTTTTT",
        constraints=[
            AvoidChanges(location=(0, 4)),
            EnforceGCContent(mini=0.8),
        ],
        logger=None,
    )
    with pytest.raises(NoSolutionError) as err:
        problem.resolve_constraints()
    assert "Exhaustive search failed" in str(err.value)
