from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceSequence
import numpy


# Note: we are not providing a location for AvoidChanges: it applies globally

def test_EnforceSequence():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    numpy.random.seed(1234)
    for symbol, nucleotides in [('W', 'AT'), ('S', 'GC')]:
        n_nucleotides = 15
        start = 50
        location = (start, start + n_nucleotides)
        problem = DnaOptimizationProblem(
            sequence=25*"ATGC",
            constraints=[AvoidPattern("ATGC"), AvoidPattern("AA"),
                         AvoidPattern("GG"),
                         EnforceSequence(n_nucleotides*symbol,
                                         location=location)]
        )
        problem.max_random_iters = 10000
        problem.resolve_constraints()
        s, e = start, start + n_nucleotides
        assert all([n in nucleotides for n in problem.sequence[s:e]])
