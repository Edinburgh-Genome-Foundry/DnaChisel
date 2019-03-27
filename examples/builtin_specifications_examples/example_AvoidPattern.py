"""Example of use of the AvoidPAttern specification"""

from dnachisel import DnaOptimizationProblem, random_dna_sequence, AvoidPattern

problem = DnaOptimizationProblem(sequence=random_dna_sequence(10000, seed=123),
                                 constraints=[AvoidPattern("BsaI_site")])

print ("\nBefore optimization, the sequence has BsaI sites (fail):\n")
print (problem.constraints_text_summary())

problem.resolve_constraints(final_check=False)

print ("\nAfter optimization, all sites are removed:\n")
print (problem.constraints_text_summary())
