import dnachisel as dc
import numpy as np


def test_circular_sequence_basic():
    np.random.seed(123)
    # Until the feature gets more battle-test, we're making sure it works
    # across a range of sequences.
    for i in range(10):
        dna_sequence = (
            "CTC"
            + dc.random_dna_sequence(1000)
            + "CGTCTC"
            + dc.random_dna_sequence(1000)
            + "CGT"
        )
        problem = dc.CircularDnaOptimizationProblem(
            sequence=dna_sequence,
            constraints=[
                dc.AvoidPattern("BsmBI_site"),
                dc.EnforceGCContent(
                    mini=0.4, maxi=0.6, location=(1500, 2500), window=50
                ),
                dc.AvoidNonUniqueSegments(
                    min_length=9,
                    location=(10, 1000),
                    extended_location=(20, 500),
                ),
            ],
            logger=None,
        )
        assert not problem.all_constraints_pass()
        problem.resolve_constraints()
        assert problem.all_constraints_pass()
