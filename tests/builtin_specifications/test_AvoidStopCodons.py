from dnachisel import DnaOptimizationProblem, AvoidStopCodons, translate
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally


def test_AvoidStopCodons():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence="".join(["ATT", "TAG", "GCC", "TGA", "ATC", "TAA"]),
        constraints=[AvoidStopCodons()],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert "*" not in translate(problem.sequence)
