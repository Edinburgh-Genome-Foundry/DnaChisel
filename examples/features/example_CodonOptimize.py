"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_protein_sequence,
                       CodonOptimize, reverse_translate, EnforceTranslation)

protein = random_protein_sequence(3000, seed=123)
sequence = reverse_translate(protein)
problem = DnaOptimizationProblem(sequence=sequence,
                                 constraints=[EnforceTranslation()],
                                 objectives=[CodonOptimize(species='e_coli')])

print ("\nBefore optimization:\n")
print (problem.objectives_text_summary())

import cProfile
problem.optimize(progress_bars=2)

print ("\nAfter optimization:\n")
print (problem.objectives_text_summary())
