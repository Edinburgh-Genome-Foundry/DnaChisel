import os
import numpy
from dnachisel import (
    CircularDnaOptimizationProblem,
    DnaOptimizationProblem,
    random_dna_sequence,
    sequence_to_biopython_record,
    Specification,
    annotate_record,
)


def test_circular_example():
    """This example has a BsmBI cross origin site (location -3 -- 3)."""
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


def test_rca_example():
    """Test a Genbank with ~harmonize_rca feature."""
    path = os.path.join("tests", "tests_from_genbanks", "genbanks", "rca_example.gb")
    problem = DnaOptimizationProblem.from_record(path)
    assert str(problem.objectives) == "[HarmonizeRCA[0-105(+)](e_coli -> h_sapiens)]"
    assert problem.objectives[0].original_species == "e_coli"
    assert problem.objectives[0].species == "h_sapiens"
    problem.optimize()


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


def test_record_with_multispec_feature():
    sequence = random_dna_sequence(100)
    record = sequence_to_biopython_record(sequence)
    label = "@gc(40-60%/20bp) & @no(BsaI_site) & @keep"
    annotate_record(record, label=label)
    problem = DnaOptimizationProblem.from_record(record)
    assert len(problem.constraints) == 3
    c1, c2, c3 = problem.constraints
    assert c1.mini == 0.4
    assert c2.pattern.name == "BsaI"


def test_feature_to_spec():
    sequence = random_dna_sequence(100)
    record = sequence_to_biopython_record(sequence)
    label = "@gc(40-60%/20bp) & @no(BsaI_site) & @keep"
    annotate_record(record, label=label)
    feature = record.features[0]
    specs = Specification.list_from_biopython_feature(feature)
    assert len(specs) == 3
