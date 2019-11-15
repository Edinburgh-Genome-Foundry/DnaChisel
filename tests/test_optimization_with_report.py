import os
import matplotlib

matplotlib.use("Agg")

from dnachisel import (
    random_dna_sequence,
    DnaOptimizationProblem,
    AvoidPattern,
    AvoidChanges,
)


def test_optimize_with_report(tmpdir):
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern("BsmBI_site")],
        logger=None,
    )

    target = os.path.join(str(tmpdir), "with_solution")
    os.mkdir(target)
    assert os.listdir(target) == []
    success, message, data = problem.optimize_with_report(target)
    assert success
    assert os.listdir(target) != []


def test_optimize_with_report_no_solution(tmpdir):
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern("BsmBI_site"), AvoidChanges()],
        logger=None,
    )
    target = os.path.join(str(tmpdir), "no_solution")
    os.mkdir(target)
    assert os.listdir(target) == []
    success, message, data = problem.optimize_with_report(target)
    assert not success
    assert os.listdir(target) != []
