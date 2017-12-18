from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidNonuniqueSegments)
import numpy
numpy.random.seed(123)

# Note: we are not providing a location for AvoidChanges: it applies globally

def test_AvoidNonuniqueSegments():

    sequence = random_dna_sequence(1000, seed=123)
    specification = AvoidNonuniqueSegments(8, location=(0, 1000))
    problem = DnaOptimizationProblem(sequence=sequence,
                                     objectives=[specification])
    problem.optimize()
    assert specification.evaluate(problem).passes
