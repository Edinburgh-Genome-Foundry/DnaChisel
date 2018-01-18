from dnachisel import (DnaOptimizationProblem, EnforceRegionsCompatibility,
                       EnforceGCContent, random_dna_sequence,
                       sequences_differences)
import numpy

# We setup the randomizer to always get the same sequence
numpy.random.seed(123)


def compatibility_condition(location1, location2, problem):

    seq1 = location1.extract_sequence(problem.sequence)
    seq2 = location2.extract_sequence(problem.sequence)
    return sequences_differences(seq1, seq2) >= 2


problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(200, seed=123),
    constraints=[
        EnforceRegionsCompatibility(
            locations=[(0, 4), (50, 54), (100, 104), (150, 154)],
            compatibility_condition=compatibility_condition,
            condition_label='2bp difference'
        ),
        EnforceGCContent(mini=0.4, maxi=0.6, window=40)
    ]
)

print ("\n\n=== Status before optimization ===")
print(problem.constraints_text_summary())

print ("Now solving constraints...")
problem.resolve_constraints()
print ("Done.")

print ("\n\n=== Status after optimization ===\n")
print (problem.constraints_text_summary())
