from dnachisel import DnaOptimizationProblem, load_record
import os


def test_cuba_example_1():
    path = os.path.join('tests', 'tests_from_genbanks', 'genbanks',
                        'cuba_example_1.gbk')
    record = load_record(path)
    problem = DnaOptimizationProblem.from_record(record)
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert problem.objective_scores_sum() < -100
    problem.optimize()
    assert problem.objective_scores_sum() > - 0.1
