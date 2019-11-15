"""Example of use of the AvoidPAttern specification"""

from dnachisel import (
    DnaOptimizationProblem,
    random_dna_sequence,
    AvoidPattern,
    EnforceGCContent,
)
import numpy


def test_EnforceGCContents():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceGCContent(mini=0.3, maxi=0.7, window=50),
        ],
        objectives=[EnforceGCContent(target=0.4)],
        logger=None,
    )

    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_parameters_from_string():
    for pattern, expected in [
        ("35%", (None, None, 0.35, None)),
        ("35%/20bp", (None, None, 0.35, 20)),
        ("5-55%", (0.05, 0.55, None, None)),
        ("5-55%/400bp", (0.05, 0.55, None, 400)),
    ]:
        mini, maxi, target, w = EnforceGCContent.string_to_parameters(pattern)
        assert (mini, maxi, target, w) == expected
