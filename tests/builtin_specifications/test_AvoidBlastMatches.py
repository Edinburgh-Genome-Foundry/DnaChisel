"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

from dnachisel import (AvoidBlastMatches, random_dna_sequence,
                       DnaOptimizationProblem)
import numpy
numpy.random.seed(123)


def test_avoid_blast_matches():
    avoided_seqs = [random_dna_sequence(20, seed=i) for i in range(10)]
    constraint = AvoidBlastMatches(sequences=avoided_seqs, min_align_length=8)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(4500, seed=12),
        constraints=[constraint]
    )
    assert not problem.all_constraints_pass()
    cst_eval = constraint.evaluate(problem)
    assert len(cst_eval.locations) == 2
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
