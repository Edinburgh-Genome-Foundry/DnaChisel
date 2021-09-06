"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

from dnachisel import (
    EnforceTranslation,
    DnaOptimizationProblem,
    random_dna_sequence,
    Location,
    EnforcePatternOccurence,
)
import dnachisel as dc
import numpy

# numpy.random.seed(123)

# Note: we are not providing a location for AvoidChanges: it applies globally


def test_enforce_pattern_basics():
    numpy.random.seed(123)
    for seed in [2, 3, 123456]:
        # The seeds cover various cases:
        # 2: the problem has no occurences instead of 1 wanted
        # 3: the pattern has no occurences instead of 1 wanted
        # 123456: the pattern is over-represented (4 times instead of 1)
        sequence = random_dna_sequence(5000, seed=seed)

        constraints = [
            EnforceTranslation(location=Location(1000, 2500)),
            EnforceTranslation(location=Location(3000, 4500)),
            EnforcePatternOccurence("ANANANANTT", location=Location(1100, 2150)),
        ]

        problem = DnaOptimizationProblem(
            sequence=sequence, constraints=constraints, logger=None
        )
        assert not problem.all_constraints_pass()
        problem.resolve_constraints()
        assert problem.all_constraints_pass()


def test_insert_and_erase_pattern():
    numpy.random.seed(123)
    protein = dc.random_protein_sequence(100)
    pattern = "ATGC"

    # CREATE A SEQUENCE WITH 0 PATTERN OCCURENCES

    sequence = dc.random_compatible_dna_sequence(
        sequence_length=300,
        constraints=[
            dc.EnforceTranslation(translation=protein),
            dc.AvoidPattern(pattern),
        ],
        logger=None,
    )

    # NOW INCREASE PATTERN OCCURENCES FROM 0 TO 5

    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            dc.EnforcePatternOccurence(pattern, occurences=5),
            dc.EnforceTranslation(),
        ],
        logger=None,
    )
    assert problem.constraints[0].evaluate(problem).score == -5
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    sequence = problem.sequence

    # NOW DECREASE THE NUMBER OF OCCURENCES FROM 5 TO 2

    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            dc.EnforcePatternOccurence(pattern, occurences=2),
            dc.EnforceTranslation(),
        ],
        logger=None,
    )
    assert problem.constraints[0].evaluate(problem).score == -3
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
