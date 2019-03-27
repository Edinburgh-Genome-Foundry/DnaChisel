"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       reverse_complement, AvoidHairpins)

random_sequences = [random_dna_sequence(30) for i in range(10)]

full_sequence = "".join([
    seq
    for sequence in random_sequences
    for seq in (random_dna_sequence(50),
                sequence,
                random_dna_sequence(50),
                reverse_complement(sequence),
                random_dna_sequence(50))
])

problem = DnaOptimizationProblem(full_sequence, constraints=[AvoidHairpins()])

print ("\nBefore optimization:\n")
print (problem.constraints_text_summary())

problem.resolve_constraints()

print ("\nAfter optimization:\n")

print (problem.constraints_text_summary())
