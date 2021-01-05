from dnachisel import (
    DnaOptimizationProblem,
    random_dna_sequence,
    UniquifyAllKmers,
)
import numpy

# Note: we are not providing a location for AvoidChanges: it applies globally
def test_UniquifyAllKmers_as_constraint():
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000, seed=123)
    problem = DnaOptimizationProblem(
        sequence=sequence, constraints=[UniquifyAllKmers(8)], logger=None
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_UniquifyAllKmers_from_polyAs():
    problem = DnaOptimizationProblem(
        sequence=40 * "A",
        constraints=[UniquifyAllKmers(3, location=(10, 30))],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_UniquifyAllKmers_as_objective():
    numpy.random.seed(123)
    sequence = random_dna_sequence(1000, seed=123)
    specification = UniquifyAllKmers(8)
    problem = DnaOptimizationProblem(
        sequence=sequence, objectives=[specification], logger=None
    )
    problem.optimize()
    assert problem.objectives[0].evaluate(problem).passes


def test_UniquifyAllKmers_from_polyAs_uncached():
    """Uncaching actually calls another function get_kmer_extractor."""
    constraint = UniquifyAllKmers(3, location=(10, 30))
    constraint.use_cache = False
    problem = DnaOptimizationProblem(
        sequence=40 * "A", constraints=[constraint], logger=None
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_UniquifyAllKmers_include_reverse_complement_false():
    constraint = UniquifyAllKmers(10, include_reverse_complement=False)
    problem = DnaOptimizationProblem(sequence=40 * "A", constraints=[constraint])
    problem.constraints_text_summary()
