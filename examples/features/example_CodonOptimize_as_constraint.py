"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       CodonOptimize, Location, EnforceTranslation)

problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(2000, seed=123),
    constraints=[
        EnforceTranslation(location=Location(1000, 1300)),
        CodonOptimize(location=Location(1000, 1300), species='e_coli')
    ]
)

print ("\nBefore resolution:\n")
print (problem.constraints_text_summary())

problem.resolve_constraints()

print ("\nAfter resolution:\n")
print (problem.constraints_text_summary())
