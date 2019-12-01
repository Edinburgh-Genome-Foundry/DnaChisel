from dnachisel import (
    DnaOptimizationProblem,
    EnforceTranslation,
    translate,
    AvoidRareCodons,
    reverse_complement,
)
import numpy


def test_AvoidRareCodons_as_constraint():
    numpy.random.seed(123)

    sequence = "ATG" "TTT" "ATA" "CCA" "CTT" "TAG"
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation(), AvoidRareCodons(0.11, "e_coli")],
    )
    assert problem.all_constraints_pass()
    assert problem.sequence_edits_as_array().sum() == 4
    assert translate(problem.sequence) == translate(sequence)


def test_AvoidRareCodons_as_constraint_reversed():
    numpy.random.seed(123)

    sequence = "ATG" "TTT" "ATA" "CCA" "CTT" "TAG"
    rev_sequence = reverse_complement(sequence)
    location = (0, len(sequence), -1)
    problem = DnaOptimizationProblem(
        sequence=rev_sequence,
        constraints=[
            EnforceTranslation(location=location),
            AvoidRareCodons(0.11, "e_coli", location=location),
        ],
    )
    assert problem.all_constraints_pass()
    assert problem.sequence_edits_as_array().sum() == 4
    new_sequence = reverse_complement(problem.sequence)
    assert translate(new_sequence) == translate(sequence)


def test_AvoidRareCodons_as_objective():
    numpy.random.seed(123)

    sequence = "ATG" "TTT" "ATA" "CCA" "CTT" "TAG"
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[AvoidRareCodons(0.11, "e_coli")],
    )
    assert abs(problem.objective_scores_sum() + 0.09) < 0.001
    problem.optimize()
    assert problem.objective_scores_sum() == 0


def test_AvoidRareCodons_as_objective_reversed():
    numpy.random.seed(123)

    sequence = "ATG" "TTT" "ATA" "CCA" "CTT" "TAG"
    rev_sequence = reverse_complement(sequence)
    location = (0, len(sequence), -1)
    problem = DnaOptimizationProblem(
        sequence=rev_sequence,
        constraints=[EnforceTranslation(location=location)],
        objectives=[AvoidRareCodons(0.11, "e_coli", location=location)],
    )
    assert abs(problem.objective_scores_sum() + 0.09) < 0.001
    problem.optimize()
    assert problem.objective_scores_sum() == 0
