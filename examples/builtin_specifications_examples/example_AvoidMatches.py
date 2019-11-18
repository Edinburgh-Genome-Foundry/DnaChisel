"""Example of use for AvoidBlastMatches.

In this example we create a 1000bp random sequence, then edit out every match
with E. coli that is 14bp or longer.

"""
from dnachisel import DnaOptimizationProblem, random_dna_sequence, AvoidMatches
from genome_collector import GenomeCollection

# THIS CREATES THE ECOLI BLAST DATABASE ON YOUR MACHINE IF NOT ALREADY HERE

collection = GenomeCollection()
ecoli_index = collection.get_taxid_bowtie_index_path(511145, version="1")

# DEFINE AND SOLVE THE PROBLEM

problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(500, seed=123),
    constraints=[
        AvoidMatches(bowtie_index=ecoli_index, match_length=15, mismatches=1)
    ],
)

print(
    "Constraints validity before optimization\n",
    problem.constraints_text_summary(),
)

print("\nNow resolving the problems\n")
problem.resolve_constraints(final_check=True)

print(
    "Constraints validity after optimization\n",
    problem.constraints_text_summary(),
)
