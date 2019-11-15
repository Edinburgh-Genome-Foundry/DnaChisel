"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

import os
from dnachisel import (
    AvoidBlastMatches,
    random_dna_sequence,
    DnaOptimizationProblem,
    load_record,
)
from genome_collector import GenomeCollection

sequence_path = os.path.join("tests", "data", "example_sequence.gbk")
sequence = str(load_record(sequence_path).seq.upper())


def test_avoid_blast_matches_with_list():
    avoided_seqs = [
        "GTCCTCATGCGAAAGCTACGATCGCCAACCCTGT",
        "ACCCACCTCGTTACGTCCACGGCACGAGGAATGATCTCGAGTTGCTTT",
    ]
    constraint = AvoidBlastMatches(
        sequences=avoided_seqs, min_align_length=8, word_size=7
    )
    problem = DnaOptimizationProblem(
        sequence=sequence, constraints=[constraint],
        logger=None
    )
    assert not problem.all_constraints_pass()
    cst_eval = constraint.evaluate(problem)
    assert len(cst_eval.locations) == 10
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_avoid_phage_blast_matches():
    PHAGE_TAXID = "697289"
    collection = GenomeCollection()
    blastdb = collection.get_taxid_blastdb_path(PHAGE_TAXID, db_type="nucl")
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(30, seed=123),
        constraints=[
            AvoidBlastMatches(
                blast_db=blastdb, min_align_length=10, word_size=7
            )
        ],
        logger=None
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
