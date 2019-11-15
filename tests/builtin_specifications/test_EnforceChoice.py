from dnachisel import DnaOptimizationProblem, EnforceChoice, EnforceGCContent
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally


def test_EnforceChoice():
    # Two enzymes, BsmBI(CGTCTC) is GC-rich, EcoRI(GAATTC) is GC-poor, which
    # enzyme will be chosen and inserted in the sequence depends on the other
    # constraint on GC content
    numpy.random.seed(123)
    spec = EnforceChoice(choices=["BsmBI_site", "EcoRI_site"], location=(2, 8))

    problem = DnaOptimizationProblem(
        sequence="AGCCCCCCGT",
        constraints=[spec, EnforceGCContent(maxi=0.3)],
        logger=None,
    )
    problem.resolve_constraints()
    assert "GAATTC" in problem.sequence

    problem = DnaOptimizationProblem(
        sequence="AGCCCCCCGT",
        constraints=[spec, EnforceGCContent(mini=0.7)],
        logger=None,
    )
    problem.resolve_constraints()
    assert "CGTCTC" in problem.sequence
