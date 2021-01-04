import dnachisel as dc
import numpy as np

# TEST PARAMETERS AUTO-SETTINGS

def test_parameterization():
    def all_none(variables):
        return all([c is None for c in variables])

    problem1 = dc.DnaOptimizationProblem(
        sequence=200 * "A",
        constraints=[
            dc.EnforceChanges(),
            dc.EnforceChanges(minimum=20),
            dc.EnforceChanges(minimum_percent=5),
        ],
        objectives=[
            dc.EnforceChanges(),
            dc.EnforceChanges(amount=20),
            dc.EnforceChanges(amount_percent=5),
        ],
    )

    record = dc.sequence_to_biopython_record(200 * "A")
    dc.annotate_record(record, label="@change")
    dc.annotate_record(record, label="@change(minimum=20)")
    dc.annotate_record(record, label="@change(minimum=5%)")
    dc.annotate_record(record, label="~change")
    dc.annotate_record(record, label="~change(amount=20)")
    dc.annotate_record(record, label="~change(5%)")
    problem2 = dc.DnaOptimizationProblem.from_record(record)

    for problem in [problem1, problem2]:

        # CHECK CONSTRAINTS

        c100 = problem.constraints[0]
        assert c100.minimum == 200
        assert c100.minimum_percent == 100
        assert all_none([c100.amount, c100.amount_percent])

        c20 = problem.constraints[1]
        assert c20.minimum == 20
        assert all_none([c20.minimum_percent, c20.amount, c20.amount_percent])

        c5 = problem.constraints[2]
        assert c5.minimum == 10
        assert c5.minimum_percent == 5
        assert all_none([c5.amount, c5.amount_percent])

        # CHECK OBJECTIVES

        o100 = problem.objectives[0]
        assert o100.amount == 200
        assert o100.amount_percent == 100
        assert all_none([o100.minimum, o100.minimum_percent])

        o20 = problem.objectives[1]
        assert o20.amount == 20
        assert all_none([o20.minimum_percent, o20.minimum, o20.amount_percent])

        o5 = problem.objectives[2]
        assert o5.amount == 10
        assert o5.amount_percent == 5
        assert all_none([o5.minimum, o5.minimum_percent])

# TEST OBJECTIVES

def test_whole_sequence_change_objective_100():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50), objectives=[dc.EnforceChanges()]
    )
    problem.optimize()
    assert problem.number_of_edits() == 50


def test_whole_sequence_change_objective_40():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50),
        objectives=[dc.EnforceChanges(amount_percent=40)],
    )
    problem.optimize()
    assert problem.number_of_edits() == 20

def test_whole_sequence_change_objective_4():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50),
        objectives=[dc.EnforceChanges(amount=4)],
    )
    problem.optimize()
    assert problem.number_of_edits() == 4

def test_whole_sequence_change_objective_20_going_down():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=20*"AT",
        constraints=[dc.AvoidPattern("ATA")],
        objectives=[dc.EnforceChanges(amount=20)],
    )
    problem.mutations_per_iteration = 2
    problem.resolve_constraints()
    assert problem.number_of_edits() >= 24
    problem.optimize()
    assert problem.number_of_edits() == 20

# TEST CONSTRAINTS

def test_whole_sequence_change_constraint_100():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50), constraints=[dc.EnforceChanges()]
    )
    assert problem.all_constraints_pass()  # due to initial seq. constraining
    assert problem.number_of_edits() == 50


def test_whole_sequence_change_constraint_40():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50),
        constraints=[dc.EnforceChanges(minimum_percent=40)],
    )
    print(problem.number_of_edits())
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert 25 >= problem.number_of_edits() >= 20

def test_whole_sequence_change_constraint_4():
    np.random.seed(123)
    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(50),
        constraints=[dc.EnforceChanges(minimum=4)],
    )
    print(problem.number_of_edits())
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert 6 >= problem.number_of_edits() >= 4

# TEST SCENARIOS

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

# TEST WITH INDICES

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
