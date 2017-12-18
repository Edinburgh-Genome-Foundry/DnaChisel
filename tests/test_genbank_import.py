from dnachisel import DnaOptimizationProblem, load_record
import os

example_sequence_path = os.path.join('tests', 'data', 'example_sequence.gbk')

def test_genbank_import():
    record = load_record(example_sequence_path)
    problem = DnaOptimizationProblem.from_record(record)
    assert len(problem.constraints) == 5
    assert len(problem.objectives) == 3
