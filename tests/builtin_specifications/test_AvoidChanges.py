"""Example of use of the AvoidChanges as an objective to minimize modifications
of a sequence."""

from dnachisel import (
    DnaOptimizationProblem,
    CircularDnaOptimizationProblem,
    random_dna_sequence,
    AvoidChanges,
    AvoidPattern,
    EnforcePatternOccurence,
    sequences_differences,
    EnforceGCContent,
    EnforceChanges,
)
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally


def test_avoid_change_as_objectives_basics():
    numpy.random.seed(123)
    results = []
    for boost in (0, 0.1, 0.2, 1):
        sequence = random_dna_sequence(1000, seed=123)
        problem = DnaOptimizationProblem(
            sequence=sequence,
            objectives=[
                EnforceGCContent(
                    mini=0.45, maxi=0.55, window=80
                ).copy_with_changes(locations_span=300),
                AvoidChanges(boost=boost).as_passive_objective(),
            ],
            logger=None,
        )

        problem.optimize()
        differences = sequences_differences(
            problem.sequence, problem.sequence_before
        )
        results.append(differences)
    assert results[0] > 40
    assert results[0] > results[1] > results[2] > results[3]
    assert results[-1] == 0


def test_avoid_change_circular():
    numpy.random.seed(123)
    results = []
    sequence = random_dna_sequence(500, seed=123)
    for boost in (0, 0.1, 0.2, 1):
        problem = CircularDnaOptimizationProblem(
            sequence=sequence,
            objectives=[
                EnforceGCContent(
                    mini=0.45, maxi=0.55, window=80
                ).copy_with_changes(locations_span=300),
                AvoidChanges(boost=boost).as_passive_objective(),
            ],
            logger=None,
        )

        problem.optimize()
        differences = sequences_differences(
            problem.sequence, problem.sequence_before
        )
        results.append(differences)
    assert results[0] > 40
    assert results[0] > results[1] > results[2] > results[3]
    assert results[-1] == 0


def test_avoid_changes_with_indices_as_constraint():
    numpy.random.seed(123)

    indices = [10, 20] + list(range(30, 40)) + [44, 45, 46]
    sequence = random_dna_sequence(50)

    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidChanges(indices=indices)],
        objectives=[EnforceChanges()],
        logger=None,
    )
    problem.optimize()
    assert problem.number_of_edits() == 50 - 15


def test_avoid_changes_with_indices_as_objectives():
    numpy.random.seed(123)

    indices = [10, 20] + list(range(30, 40)) + [44, 45, 46]
    sequence = random_dna_sequence(50)

    problem = DnaOptimizationProblem(
        sequence=sequence,
        objectives=[EnforceChanges(boost=0.5), AvoidChanges(indices=indices)],
        logger=None,
    )
    problem.optimize()
    assert problem.number_of_edits() == 50 - 15  # 15 == len(indices)

    problem = DnaOptimizationProblem(
        sequence=sequence,
        objectives=[EnforceChanges(boost=1.5), AvoidChanges(indices=indices)],
        logger=None,
    )
    problem.optimize()
    assert problem.number_of_edits() == 50


def test_AvoidChanges_with_max_edits():
    numpy.random.seed(1)
    problem = DnaOptimizationProblem(
        sequence="ATATATATATA",
        constraints=[
            AvoidChanges(max_edits=2),
            AvoidPattern("ATATA"),
            EnforcePatternOccurence("A", occurences=6, location=(0, 11, 1)),
            EnforcePatternOccurence("T", occurences=4, location=(0, 11, 1)),
        ],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
