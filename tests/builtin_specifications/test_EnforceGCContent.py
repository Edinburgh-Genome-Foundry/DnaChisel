"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidPattern, EnforceGCContent)

def test_basics():

    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[
            AvoidPattern(enzyme="BsaI"),
            EnforceGCContent(gc_min=0.3, gc_max=0.7, gc_window=50)
        ],
        objectives=[EnforceGCContent(gc_objective=0.4)]
    )

    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
