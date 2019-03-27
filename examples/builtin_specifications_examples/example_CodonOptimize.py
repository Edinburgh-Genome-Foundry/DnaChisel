"""Example of use of the CodonOptimize specification."""

from dnachisel import (DnaOptimizationProblem, random_protein_sequence,
                       CodonOptimize, reverse_translate, EnforceTranslation)

protein = random_protein_sequence(3000, seed=123)
sequence = reverse_translate(protein)
problem = DnaOptimizationProblem(
    sequence=sequence,
    constraints=[EnforceTranslation()],
    objectives=[CodonOptimize('e_coli')])

print ("\nBefore optimization:\n")
print (problem.objectives_text_summary())

problem.optimize()

print ("\nAfter optimization:\n")
print (problem.objectives_text_summary())
