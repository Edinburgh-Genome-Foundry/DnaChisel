import itertools
from dnachisel import (
    DnaOptimizationProblem,
    EnforceGCContent,
    EnforceRegionsCompatibility,
    sequences_differences,
    random_dna_sequence,
)
import numpy


# Note: we are not providing a location for AvoidChanges: it applies globally


def test_EnforceRegionsCompatibility():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    numpy.random.seed(123)

    def compatibility_condition(location1, location2, problem):
        seq1 = location1.extract_sequence(problem.sequence)
        seq2 = location2.extract_sequence(problem.sequence)
        return sequences_differences(seq1, seq2) >= 2

    locations = [(0, 4), (50, 54), (100, 104), (150, 154)]
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(200, seed=123),
        constraints=[
            EnforceRegionsCompatibility(
                locations=locations,
                compatibility_condition=compatibility_condition,
                condition_label="2bp difference",
            ),
            EnforceGCContent(mini=0.4, maxi=0.6, window=40),
        ],
        logger=None,
    )
    assert not any([e.passes for e in problem.constraints_evaluations()])
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    seq = problem.sequence
    assert [
        sequences_differences(seq[s1:e1], seq[s2:e2]) >= 2
        for (s1, e1), (s2, e2) in itertools.combinations(locations, 2)
    ]
