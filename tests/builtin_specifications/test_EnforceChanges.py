import dnachisel as dc
import numpy as np


def test_total_sequence_change():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(100), objectives=[dc.EnforceChanges()]
    )
    problem.optimize()
    assert problem.number_of_edits() == 100


def test_maximal_protein_sequence_change():
    np.random.seed(123)
    protein = dc.random_protein_sequence(200)
    sequence = dc.reverse_translate(protein)
    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[dc.EnforceTranslation()],
        objectives=[dc.EnforceChanges()],
    )
    problem.resolve_constraints()
    problem.optimize()
    assert problem.number_of_edits() == 238
    assert dc.translate(problem.sequence) == protein


def test_enforce_changes_as_constraint():
    # This test checks that using EnforceChanges as constraint will
    # automatically change the whole sequence at initialization time.
    np.random.seed(123)

    sequence = dc.random_dna_sequence(500)
    problem = dc.DnaOptimizationProblem(
        sequence=sequence, constraints=[dc.EnforceChanges()]
    )
    assert problem.all_constraints_pass()
    assert problem.number_of_edits() == 500
    # check the randomness / non-bias of the new selected sequences
    assert 0.45 < dc.biotools.gc_content(problem.sequence) < 0.55


def test_enforce_changes_with_indices_as_constraint():
    np.random.seed(123)
    indices = [10, 20] + list(range(30, 40)) + [44, 45, 46]
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50),
        constraints=[dc.EnforceChanges(indices=indices)],
    )
    assert problem.number_of_edits() == 15


def test_enforce_changes_with_indices_vs_avoid_changes():
    np.random.seed(123)
    indices = [10, 20] + list(range(30, 40)) + [44, 45, 46]
    sequence = dc.random_dna_sequence(50)
    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        objectives=[
            dc.EnforceChanges(indices=indices),
            dc.AvoidChanges(boost=0.5),
        ],
    )
    problem.optimize()
    assert problem.number_of_edits() == 15

    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        objectives=[
            dc.EnforceChanges(indices=indices),
            dc.AvoidChanges(boost=1.5),
        ],
    )
    problem.optimize()
    assert problem.number_of_edits() == 0
