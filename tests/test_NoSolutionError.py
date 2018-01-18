import pytest
from dnachisel import (DnaOptimizationProblem, AvoidChanges, AvoidPattern,
                       NoSolutionError)

def test_no_solution_error():
    problem = DnaOptimizationProblem(
        sequence="AAAAATCGTCTCTTTT",
        constraints=[AvoidChanges(), AvoidPattern(enzyme='BsmBI')]
    )
    with pytest.raises(NoSolutionError):
        problem.resolve_constraints()
