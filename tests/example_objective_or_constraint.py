from dnachisel import *
import numpy

# We setup the randomizer to always get the same sequence
from dnachisel import (random_dna_sequence, DnaOptimizationProblem,
                       AvoidPattern)
numpy.random.seed(123)
sequence = random_dna_sequence(10000)
spec = AvoidPattern("CG")
problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(10000),
    constraints=[
        AvoidPattern(enzyme_pattern("BsaI")),
        EnforceGCContent(mini=0.3, maxi=0.7, window=50)
    ],
    objectives=[EnforceGCContent(gc_objective=0.4)]
)

print ("\n\n=== Status before optimization ===")
print(problem.constraints_text_summary())
print(problem.objectives_text_summary())

print ("Now solving constraints...")
problem.solve_constraints(progress_bars=1)
print ("Done. Now optimizing objectives...")
problem.maximize_objectives(max_random_iters=10000)

print ("\n\n=== Status after optimization ===\n")
print (problem.constraints_text_summary())
print (problem.objectives_text_summary())
