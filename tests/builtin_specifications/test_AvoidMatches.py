import os
from dnachisel import (
    AvoidMatches,
    random_dna_sequence,
    DnaOptimizationProblem,
    load_record,
)
from genome_collector import GenomeCollection

sequence_path = os.path.join("tests", "data", "example_sequence.gbk")
sequence = str(load_record(sequence_path).seq.upper())


def test_avoid_matches_with_list():
    pattern_1 = "CGTCTC"
    pattern_2 = "TGCACA"
    sequence = 10 * "A" + pattern_1 + 20 * "A" + pattern_2 + 10 * "A"
    avoided_seqs = [
        10 * "G" + pattern_1 + 10 * "G",
        10 * "G" + pattern_2 + 10 * "G",
    ]
    constraint = AvoidMatches(sequences=avoided_seqs, match_length=6)
    problem = DnaOptimizationProblem(
        sequence=sequence, constraints=[constraint], logger=None
    )
    cst_eval = constraint.evaluate(problem)
    assert len(cst_eval.locations) == 2
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    constraint.remove_temp_directory()


def test_avoid_matches_with_phage():
    PHAGE_TAXID = "697289"
    collection = GenomeCollection()
    index = collection.get_taxid_bowtie_index_path(PHAGE_TAXID, version="1")
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(30, seed=123),
        constraints=[AvoidMatches(bowtie_index=index, match_length=10)],
        logger=None,
    )
    all_breaches = problem.constraints_evaluations().all_locations()
    assert len(all_breaches) == 5
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
