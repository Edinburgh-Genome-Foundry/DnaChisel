"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidPattern, EnforceGCContent)
import numpy

def test_EnforceGCContents():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[
            AvoidPattern(enzyme="BsaI"),
            EnforceGCContent(mini=0.3, maxi=0.7, window=50)
        ],
        objectives=[EnforceGCContent(target=0.4)]
    )

    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
