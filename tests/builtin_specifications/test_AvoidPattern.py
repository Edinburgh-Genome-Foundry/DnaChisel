"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidPattern, RepeatedKmerPattern, AvoidChanges)
import numpy

def test_avoid_pattern_basics():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern("BsaI_site")]
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

def test_avoid_repeated_small_kmers():
    problem = DnaOptimizationProblem(
        sequence="AGAAGAAGAAGAAGAAGATTTTTTTTTTTTTGGAGGAGGAGGACCCCCCCCCCCCGAGG",
        constraints=[AvoidPattern(RepeatedKmerPattern(3, 3))]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()

def test_pattern_and_reverse():
    bsmbi = "CGTCTC"
    bsmbi_rev = "GAGACG"
    sequence = 10 * bsmbi + 25 * bsmbi_rev + 15 * bsmbi + 15 * bsmbi_rev
    problem = DnaOptimizationProblem(
        sequence,
        constraints=[AvoidPattern('BsmBI_site')],
        objectives=[AvoidChanges()]
    )
    problem.resolve_constraints()
    problem.optimize()
    assert sum(problem.sequence_edits_as_array()) < 70