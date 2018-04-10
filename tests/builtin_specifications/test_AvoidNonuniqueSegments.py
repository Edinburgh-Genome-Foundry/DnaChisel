from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidNonuniqueSegments)
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally
def test_AvoidNonuniqueSegments_as_constraint():
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000, seed=123)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidNonuniqueSegments(8, location=(0, 1000))])
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()

def test_AvoidNonuniqueSegments_as_objective():
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000, seed=123)
    specification = AvoidNonuniqueSegments(8, location=(0, 1000))
    problem = DnaOptimizationProblem(sequence=sequence,
                                     objectives=[specification])
    problem.optimize()
    assert specification.evaluate(problem).passes
