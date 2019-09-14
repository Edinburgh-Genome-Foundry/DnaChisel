"""Example of use of the AvoidPAttern specification"""

from dnachisel import (
    DnaOptimizationProblem,
    random_protein_sequence,
    reverse_translate,
    CodonOptimize,
    EnforceTranslation,
    AvoidPattern,
    EnforceGCContent,
)

protein = random_protein_sequence(1000, seed=123)
sequence = reverse_translate(protein)
problem = DnaOptimizationProblem(
    sequence=sequence,
    constraints=[
        EnforceTranslation(),
        AvoidPattern("BsmBI_site"),
        EnforceGCContent(mini=0.4, maxi=0.6, window=60),
    ],
    objectives=[CodonOptimize(species="s_cerevisiae")],
)

print("\nBefore optimization:\n")
print(problem.constraints_text_summary())
print(problem.objectives_text_summary())

problem.resolve_constraints(final_check=True)
problem.optimize()

print("\nAfter optimization:\n")
print(problem.constraints_text_summary())
print(problem.objectives_text_summary())
