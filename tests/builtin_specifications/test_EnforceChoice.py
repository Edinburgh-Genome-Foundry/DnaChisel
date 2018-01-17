from dnachisel import DnaOptimizationProblem, EnforceChoice, EnforceGCContent
import numpy
numpy.random.seed(123)

# Note: we are not providing a location for AvoidChanges: it applies globally

def test_EnforceChoice():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    spec = EnforceChoice(enzymes=['BsmBI', 'EcoRI'], location=(2, 8))

    problem = DnaOptimizationProblem(
        sequence="AGCCCCCCGT",  constraints=[spec, EnforceGCContent(maxi=0.3)])
    problem.resolve_constraints()
    assert 'GAATTC' in problem.sequence

    problem = DnaOptimizationProblem(
        sequence="AGCCCCCCGT", constraints=[spec, EnforceGCContent(mini=0.7)])
    problem.resolve_constraints()
    assert 'CGTCTC' in problem.sequence
