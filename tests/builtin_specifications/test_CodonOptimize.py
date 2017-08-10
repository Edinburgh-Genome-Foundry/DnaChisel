"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_protein_sequence,
                       random_dna_sequence, Location, CodonOptimize,
                       reverse_translate, EnforceTranslation)


def test_codon_optimize_basics():
    protein = random_protein_sequence(3000, seed=123)
    sequence = reverse_translate(protein)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species='e_coli')]
    )
    assert problem.objective_scores_sum() < 0
    problem.optimize()
    assert problem.objective_scores_sum() == 0


def test_codon_optimize_as_hard_constraint():
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(2000, seed=123),
        constraints=[
            EnforceTranslation(location=Location(1000, 1300)),
            CodonOptimize(location=Location(1000, 1300), species='e_coli')
        ]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
