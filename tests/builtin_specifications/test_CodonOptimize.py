"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_protein_sequence,
                       CodonOptimize, reverse_translate, EnforceTranslation)

def test_basics():
    protein = random_protein_sequence(3000, seed=123)
    sequence = reverse_translate(protein)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species='e_coli')]
    )
    assert problem.objective_scores_sum() < 0
    problem.optimize(progress_bars=2)
    assert problem.objective_scores_sum() == 0
