from dnachisel import (
    DnaOptimizationProblem,
    random_compatible_dna_sequence,
    AvoidPattern,
    EnforceGCContent,
)


def test_random_compatible_dna_sequence():
    constraints = [
        EnforceGCContent(mini=0.4, maxi=0.6, window=50),
        AvoidPattern("ATC"),
    ]
    seq = random_compatible_dna_sequence(1000, constraints=constraints)
    problem = DnaOptimizationProblem(
        sequence=seq, constraints=constraints, logger=None
    )
    assert "ATC" not in seq
    assert problem.all_constraints_pass()
