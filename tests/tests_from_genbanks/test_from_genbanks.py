from dnachisel import (
    CircularDnaOptimizationProblem,
    DnaOptimizationProblem
)
import os


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
