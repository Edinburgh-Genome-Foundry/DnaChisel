import os
import matplotlib
matplotlib.use("Agg")

from dnachisel.reports import optimization_with_report
from dnachisel import (random_dna_sequence, DnaOptimizationProblem,
                       AvoidPattern, AvoidChanges)


def test_optimization_with_report(tmpdir):
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern(enzyme='BsmBI')]
    )

    target = os.path.join(str(tmpdir), 'with_solution')
    os.mkdir(target)
    assert os.listdir(target) == []
    success, message, data = optimization_with_report(target, problem)
    assert success
    assert os.listdir(target) != []

def test_optimization_with_report_no_solution(tmpdir):
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern(enzyme='BsmBI'), AvoidChanges()]
    )
    target = os.path.join(str(tmpdir), 'no_solution')
    os.mkdir(target)
    assert os.listdir(target) == []
    success, message, data = optimization_with_report(target, problem)
    assert not success
    assert os.listdir(target) != []
