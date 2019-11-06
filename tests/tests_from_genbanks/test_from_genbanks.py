import os
import numpy
from dnachisel import (
    CircularDnaOptimizationProblem,
    DnaOptimizationProblem,
    random_dna_sequence,
    sequence_to_biopython_record,
    annotate_record,
)


def test_circular_example():
    """This example has a BsmBI cross origin site (location -3 -- 3)"""
    path = os.path.join(
        "tests", "tests_from_genbanks", "genbanks", "circular_example_1.gb"
    )
    problem = CircularDnaOptimizationProblem.from_record(path)
    evals = problem.constraints_evaluations()
    assert str(evals.evaluations[0].locations[0]) == "-3-3(+)"
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_cuba_example_1():
    path = os.path.join(
        "tests", "tests_from_genbanks", "genbanks", "cuba_example_1.gbk"
    )
    problem = DnaOptimizationProblem.from_record(path)
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert problem.objective_scores_sum() < -100
    problem.optimize()
    assert problem.objective_scores_sum() > -0.1


def test_all_shorthands():
    """This test compiles all shorthands as a check that nothing is broken."""
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000)
    record = sequence_to_biopython_record(sequence)
    annotate_record(record, (100, 900), label="@no(CATG)")
    annotate_record(record, (100, 900), label="@gc(40-60%)")
    annotate_record(record, (100, 900), label="@insert(AarI_site)")
    annotate_record(record, (650, 752), label="@cds")
    annotate_record(record, (100, 200), label="@keep")
    annotate_record(record, (250, 273), label="@primer")
    annotate_record(record, (250, 280), label="@change")
    annotate_record(record, (943, 950), label="@sequence(AKGNTKT)")
    annotate_record(record, (955, 958), label="@sequence(ATT|ATC|GGG)")
    problem = DnaOptimizationProblem.from_record(record)
    assert len(problem.constraints) == 13  # AllowPrimer counts for 4 specs.
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
