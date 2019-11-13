from dnachisel import (
    DnaOptimizationProblem,
    random_dna_sequence,
    AvoidNonUniqueSegments,
)
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally
def test_AvoidNonUniqueSegments_as_constraint():
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000, seed=123)
    problem = DnaOptimizationProblem(
        sequence=sequence, constraints=[AvoidNonUniqueSegments(8)]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_AvoidNonUniqueSegments_from_polyAs():
    problem = DnaOptimizationProblem(
        sequence=40 * "A",
        constraints=[AvoidNonUniqueSegments(3, location=(10, 30))],
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_AvoidNonUniqueSegments_as_objective():
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000, seed=123)
    specification = AvoidNonUniqueSegments(8)
    problem = DnaOptimizationProblem(
        sequence=sequence, objectives=[specification]
    )
    problem.optimize()
    assert problem.objectives[0].evaluate(problem).passes


def test_AvoidNonUniqueSegments_from_polyAs_uncached():
    """Uncaching actually calls another function get_kmer_extractor."""
    constraint = AvoidNonUniqueSegments(3, location=(10, 30))
    constraint.use_cache = False
    problem = DnaOptimizationProblem(
        sequence=40 * "A", constraints=[constraint]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
