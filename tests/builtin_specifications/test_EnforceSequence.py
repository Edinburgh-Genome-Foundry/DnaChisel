from dnachisel import DnaOptimizationProblem, AvoidPattern, EnforceSequence
import numpy
numpy.random.seed(123)

# Note: we are not providing a location for AvoidChanges: it applies globally

def test_EnforceSequence():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    for symbol, nucleotides in [('W', 'AT'), ('S', 'GC')]:
        problem = DnaOptimizationProblem(
            sequence=50*"ATGC",
            constraints=[AvoidPattern("ATGC"), AvoidPattern("AA"),
                         AvoidPattern("GG"), 
                         EnforceSequence(30*symbol, location=(50, 80))]
        )
        problem.resolve_constraints()
        assert all([n in nucleotides for n in problem.sequence[50:80]])
