"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

import os
from dnachisel import (AvoidBlastMatches, random_dna_sequence,
                       DnaOptimizationProblem, load_record)

sequence_path = os.path.join("tests", "data", "example_sequence.gbk")
sequence = str(load_record(sequence_path).seq.upper())

def test_avoid_blast_matches():
    avoided_seqs = ["GTCCTCATGCGAAAGCTACGATCGCCAACCCTGT",
                    "ACCCACCTCGTTACGTCCACGGCACGAGGAATGATCTCGAGTTGCTTT"]
    constraint = AvoidBlastMatches(sequences=avoided_seqs, min_align_length=8)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[constraint]
    )
    assert not problem.all_constraints_pass()
    cst_eval = constraint.evaluate(problem)
    assert len(cst_eval.locations) == 10
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
