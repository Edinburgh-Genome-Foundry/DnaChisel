import dnachisel as dc


def test_SequenceLengthBounds():
    for length, expected in [(750, True), (400, False), (1200, False)]:
        problem = dc.DnaOptimizationProblem(
            sequence=dc.random_dna_sequence(length),
            constraints=[dc.SequenceLengthBounds(500, 800)],
            logger=None,
        )
        assert problem.all_constraints_pass() == expected
