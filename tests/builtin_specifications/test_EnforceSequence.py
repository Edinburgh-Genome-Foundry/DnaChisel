from dnachisel import (
    DnaOptimizationProblem,
    AvoidPattern,
    EnforceSequence,
    EnforcePatternOccurence,
    Location,
)
import numpy


# Note: we are not providing a location for AvoidChanges: it applies globally


def test_EnforceSequence():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    numpy.random.seed(1234)
    for symbol, nucleotides in [("W", "AT"), ("S", "GC")]:
        n_nucleotides = 15
        start = 50
        location = (start, start + n_nucleotides)
        problem = DnaOptimizationProblem(
            sequence=25 * "ATGC",
            constraints=[
                AvoidPattern("ATGC"),
                AvoidPattern("AAA"),
                AvoidPattern("GGG"),
                EnforceSequence(n_nucleotides * symbol, location=location),
            ],
        )
        problem.max_random_iters = 10000
        problem.resolve_constraints()
        s, e = start, start + n_nucleotides
        assert all([n in nucleotides for n in problem.sequence[s:e]])

    # Test -1 strand:
    seq = "ATG" + "CAG" + "AGCAAGGTGCTGCT"
    problem = DnaOptimizationProblem(
        sequence=seq,
        constraints=[
            EnforcePatternOccurence(
                pattern="CTG",  # CAG on strand +1
                occurences=2,
                strand=-1,
                location=Location(start=0, end=50),
            )
        ],
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_EnforceSequence_as_objective():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    numpy.random.seed(1234)
    n_nucleotides = 15
    start = 50
    location = (start, start + n_nucleotides)
    problem = DnaOptimizationProblem(
        sequence=25 * "ATGC",
        constraints=[AvoidPattern("ATGC")],
        objectives=[EnforceSequence("W" * n_nucleotides, location=location)],
        logger=None,
    )
    assert problem.objective_scores_sum() < 0
    problem.resolve_constraints()
    problem.optimize()
    assert problem.objective_scores_sum() == 0
