"""Example of use of the AvoidPAttern specification"""

from dnachisel import DnaOptimizationProblem, random_dna_sequence, AvoidPattern
import numpy

def test_avoid_pattern_basics():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern(enzyme="BsaI")]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()

def test_avoid_pattern_overlapping_locations():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence="AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        constraints=[AvoidPattern("NAN")]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert "A" not in problem.sequence[1:-1]
