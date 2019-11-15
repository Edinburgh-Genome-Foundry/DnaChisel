from dnachisel import (
    DnaOptimizationProblem,
    load_record,
    EnforceGCContent,
    AvoidPattern,
)
import pytest
import os

example_sequence_path = os.path.join("tests", "data", "example_sequence.gbk")


def test_genbank_import_from_record():
    record = load_record(example_sequence_path)
    problem = DnaOptimizationProblem.from_record(record)
    assert len(problem.constraints) == 5
    assert len(problem.objectives) == 3


def test_genbank_import_from_record_unknown_specs():
    record = load_record(example_sequence_path)
    with pytest.raises(TypeError):
        _ = DnaOptimizationProblem.from_record(
            record, specifications_dict={}
        )


def test_genbank_import_from_filepath():
    problem = DnaOptimizationProblem.from_record(example_sequence_path)
    assert len(problem.constraints) == 5
    assert len(problem.objectives) == 3


def test_constraints_text_summary():
    problem = DnaOptimizationProblem(
        sequence="ATTGCCATATGCGC",
        constraints=[
            EnforceGCContent(mini=0.4, maxi=0.6),
            AvoidPattern("ATT"),
        ],
    )
    text = problem.constraints_text_summary()
    assert "FAILURE: 1 constraints evaluations failed" in text
