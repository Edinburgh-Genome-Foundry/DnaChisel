from dnachisel import DnaOptimizationProblem, AvoidStopCodons
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally

def test_AvoidStopCodons():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence="ATTGCCATCTAA",
        constraints=[AvoidStopCodons()]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
