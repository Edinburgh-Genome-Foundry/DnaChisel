"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidChanges, sequences_differences,
                       EnforceGCContent)
import numpy
numpy.random.seed(123)

# Note: we are not providing a location for AvoidChanges: it applies globally

def test_basics():
    results = []
    for boost in (0, 0.1, 1, 10.0):
        sequence = random_dna_sequence(1000, seed=123)
        problem = DnaOptimizationProblem(
            sequence=sequence,
            objectives=[
                EnforceGCContent(gc_min=0.45, gc_max=0.55, gc_window=80),
                AvoidChanges(boost=boost).as_passive_objective()])

        problem.optimize(max_random_iters=1000, progress_bars=3)
        differences = sequences_differences(problem.sequence,
                                            problem.sequence_before)
        results.append(differences)
    assert results[0] > 40
    assert(results[0] > results[1] > results[2] > results[3])
    assert results[-1] == 0
