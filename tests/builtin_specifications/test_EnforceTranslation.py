from dnachisel import (
    DnaOptimizationProblem,
    EnforceTranslation,
    EnforceGCContent,
    random_protein_sequence,
    reverse_translate,
    reverse_complement,
    translate,
    AvoidPattern,
    EnforceChanges,
    Location,
)
import pytest
import numpy


# Note: we are not providing a location for AvoidChanges: it applies globally


def test_EnforceTranslation():
    numpy.random.seed(1234)
    sequence = reverse_translate(random_protein_sequence(50, seed=123))
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidPattern("AAA"), EnforceTranslation()],
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_EnforceTranslationReversed():
    numpy.random.seed(1234)
    sequence = reverse_translate(random_protein_sequence(50, seed=123))
    rev_sequence = reverse_complement(sequence)
    problem = DnaOptimizationProblem(
        sequence=rev_sequence,
        constraints=[
            AvoidPattern("AGC"),
            EnforceTranslation(location=(0, len(sequence), -1)),
        ],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_EnforceTranslation_error_location_not_3x():
    """Providing a location that is not multiple of 3 raises an error"""
    numpy.random.seed(1234)
    sequence = reverse_translate(random_protein_sequence(50, seed=123))
    with pytest.raises(ValueError) as err:
        _ = DnaOptimizationProblem(
            sequence=sequence,
            constraints=[EnforceTranslation(location=(0, 16))],
            logger=None,
        )
    assert "Location 0-16(+) has length 16" in str(err.value)


def test_EnforceTranslation_error_location_smaller_than_translation():
    """Providing a location that is not multiple of 3 raises an error"""
    numpy.random.seed(1234)
    sequence = reverse_translate(random_protein_sequence(15, seed=123))
    with pytest.raises(ValueError) as err:
        _ = DnaOptimizationProblem(
            sequence=sequence,
            constraints=[
                EnforceTranslation(
                    translation=random_protein_sequence(30, seed=111)
                )
            ],
            logger=None,
        )
    assert str(err.value).startswith("Window size")


def test_EnforceTranslation_bacterial_valine():
    table_name = "Bacterial"
    protein = "LLTMMVTTTTVMVL"
    protein_sequence = reverse_translate(protein, table=table_name)

    for first_codon_before, first_codon_after in [
        ("ATG", "ATG"),  # methionine stays the only methionine codon
        ("GTG", "GTG"),  # valine-start-codon stays the only valine-start-codon
    ]:
        sequence = first_codon_before + protein_sequence
        cds_constraint = EnforceTranslation(
            genetic_table="Bacterial", start_codon="keep"
        )
        problem = DnaOptimizationProblem(
            sequence=sequence,
            constraints=[cds_constraint],
            objectives=[EnforceChanges()],
            logger=None,
        )
        assert problem.constraints[0].translation == "MLLTMMVTTTTVMVL"
        problem.optimize()
        protein_after = translate(
            problem.sequence, table_name, assume_start_codon=True
        )
        assert protein_after == "M" + protein
        assert problem.sequence[:3] == first_codon_after


def test_EnforceTranslation_bacterial_valine_antisense():
    table_name = "Bacterial"
    protein = "LLTMMVTTTTVMVL"
    protein_sequence = reverse_translate(protein, table=table_name)
    for first_codon_before, first_codon_after in [
        ("ATG", "ATG"),  # methionine stays the only methionine codon
        ("GTG", "GTG"),  # valine-start-codon stays the only valine-start-codon
    ]:
        sequence = first_codon_before + protein_sequence
        cds_constraint = EnforceTranslation(
            genetic_table="Bacterial",
            start_codon="keep",
            location=Location(0, len(sequence), -1),
        )
        problem = DnaOptimizationProblem(
            sequence=reverse_complement(sequence),
            constraints=[cds_constraint],
            objectives=[EnforceChanges()],
            logger=None,
        )
        assert problem.constraints[0].translation == "MLLTMMVTTTTVMVL"
        problem.optimize()
        problem_sequence_rv = reverse_complement(problem.sequence)
        protein_after = translate(
            problem_sequence_rv, table_name, assume_start_codon=True
        )
        assert protein_after == "M" + protein
        assert problem_sequence_rv[:3] == first_codon_after
