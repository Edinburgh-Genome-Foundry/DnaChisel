"""Example of use of the AvoidPAttern specification"""

from dnachisel import (
    DnaOptimizationProblem,
    random_dna_sequence,
    AvoidPattern,
    EnforceTerminalGCContent,
)
import numpy


def test_basics():
    numpy.random.seed(123)
    probas = {"A": 0.2, "T": 0.2, "G": 0.3, "C": 0.3}
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, probas=probas, seed=123),
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceTerminalGCContent(mini=0.2, maxi=0.4, window_size=50),
        ],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
