from dnachisel import *
import numpy

# We setup the randomizer to always get the same sequence
numpy.random.seed(123)

problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(10000),
    constraints=[
        AvoidPattern("BsaI_site"),
        EnforceGCContent(mini=0.3, maxi=0.7, window=50)
    ],
    objectives=[EnforceGCContent(target=0.4)]
)

print ("\n\n=== Status before optimization ===")
print(problem.constraints_text_summary())
print(problem.objectives_text_summary())

print ("Now solving constraints...")
problem.resolve_constraints()
print ("Done. Now optimizing objectives...")
problem.optimize()

print ("\n\n=== Status after optimization ===\n")
print (problem.constraints_text_summary())
print (problem.objectives_text_summary())
