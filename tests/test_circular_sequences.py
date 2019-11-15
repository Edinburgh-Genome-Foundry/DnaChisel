import dnachisel as dc
import numpy as np
import os


def test_circular_sequence_basic():
    np.random.seed(123)
    # Until the feature gets more battle-test, we're making sure it works
    # across a range of sequences.
    for i in range(4):
        dna_sequence = (
            "CTC"
            + dc.random_dna_sequence(100)
            + "CGTCTC"
            + dc.random_dna_sequence(100)
            + "CGT"
        )
        problem = dc.CircularDnaOptimizationProblem(
            sequence=dna_sequence,
            constraints=[
                dc.AvoidPattern("BsmBI_site"),
                dc.EnforceGCContent(
                    mini=0.4, maxi=0.6, location=(150, 250), window=50
                ),
                dc.UniquifyAllKmers(k=9, location=(10, 100)),
            ],
            logger=None,
        )
        assert not problem.all_constraints_pass()
        problem.resolve_constraints()
        assert problem.all_constraints_pass()


def test_circular_sequence_optimize_with_report(tmpdir):
    """Test that the custom function of CircularDnaOptimizationProblems works.
    """
    np.random.seed(123)
    # Until the feature gets more battle-test, we're making sure it works
    # across a range of sequences.
    dna_sequence = (
        "CTC"
        + dc.random_dna_sequence(100)
        + "CGTCTC"
        + dc.random_dna_sequence(100)
        + "CGT"
    )
    problem = dc.CircularDnaOptimizationProblem(
        sequence=dna_sequence,
        constraints=[
            dc.AvoidPattern("BsmBI_site"),
            dc.EnforceGCContent(
                mini=0.4, maxi=0.6, location=(150, 250), window=50
            ),
            dc.UniquifyAllKmers(k=9, location=(10, 100)),
        ],
        logger=None,
    )

    target = os.path.join(str(tmpdir), "circular_with_solution")
    os.mkdir(target)
    assert os.listdir(target) == []
    assert not problem.all_constraints_pass()
    success, message, data = problem.optimize_with_report(target)
    assert problem.all_constraints_pass()
    record = problem.to_record()
    assert str(record.seq) != dna_sequence
