from dnachisel import (DnaOptimizationProblem, EnforceTranslation,
                       EnforceGCContent, random_protein_sequence,
                       reverse_translate, reverse_complement, AvoidPattern)
import pytest
import numpy


# Note: we are not providing a location for AvoidChanges: it applies globally

def test_EnforceTranlation():
    numpy.random.seed(1234)
    sequence = reverse_translate(random_protein_sequence(50, seed=123))
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidPattern("AAA"), EnforceTranslation()],
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()

def test_EnforceTranlationReversed():
        numpy.random.seed(1234)
        sequence = reverse_translate(random_protein_sequence(50, seed=123))
        rev_sequence = reverse_complement(sequence)
        problem = DnaOptimizationProblem(
            sequence=rev_sequence,
            constraints=[AvoidPattern("AGC"),
                         EnforceTranslation(location=(0, len(sequence), -1))],
        )
        assert not problem.all_constraints_pass()
        problem.resolve_constraints()
        assert problem.all_constraints_pass()

def test_EnforceTranlationError():
    """Providing a location that is not multiple of 3 raises an error"""
    numpy.random.seed(1234)
    sequence = reverse_translate(random_protein_sequence(50, seed=123))
    with pytest.raises(ValueError) as err:
        problem = DnaOptimizationProblem(
            sequence=sequence,
            constraints=[EnforceTranslation(location=(0, 16))],
        )
    assert "Location 0-16(+) has length 16" in str(err.value)
