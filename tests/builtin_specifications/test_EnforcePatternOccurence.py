"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

from dnachisel import (EnforceTranslation, DnaOptimizationProblem,
                       random_dna_sequence, Location, EnforcePatternOccurence)

import numpy
#numpy.random.seed(123)

# Note: we are not providing a location for AvoidChanges: it applies globally

def test_enforce_pattern_basics():
    for seed in [2, 3, 123456]:
        # The seeds cover various cases:
        # 2: the problem has no occurences instead of 1 wanted
        # 3: the pattern has no occurences instead of 1 wanted
        # 123456: the pattern is over-represented (4 times instead of 1)
        sequence = random_dna_sequence(5000, seed=seed)

        constraints = [
            EnforceTranslation(Location(1000, 2500)),
            EnforceTranslation(Location(3000, 4500)),
            EnforcePatternOccurence("ANANANANTT",
                                    location=Location(1100, 2150))
        ]

        problem = DnaOptimizationProblem(
            sequence=sequence,
            constraints=constraints
        )
        assert not problem.all_constraints_pass()
        problem.resolve_constraints()
        assert problem.all_constraints_pass()
