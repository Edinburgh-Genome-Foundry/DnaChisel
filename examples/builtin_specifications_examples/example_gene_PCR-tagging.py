"""Example of use for AvoidBlastMatches.

In this example we create a 1000bp random sequence, then edit out every match
with E. coli that is 14bp or longer.

"""
import os
from _utils import download_ecoli_genome_if_not_already_there

from dnachisel import (
    DnaOptimizationProblem,
    random_dna_sequence,
    AvoidBlastMatches,
)

# THIS CREATES THE ECOLI BLAST DATABASE ON YOUR MACHINE IF NOT ALREADY HERE
download_ecoli_genome_if_not_already_there()

# DEFINE AND SOLVE THE PROBLEM

problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(1000, seed=123),
    constraints=[
        AvoidBlastMatches(
            blast_db=os.path.join("downloaded_data", "ecoli_genome"),
            min_align_length=13,
            perc_identity=95,
        )
    ],
)
print("\n\nWarning! BLASTing can take a long time! Be patient.", end="\n\n")
print(
    "Constraints validity before optimization\n",
    problem.constraints_text_summary(),
)
print("\nLet's resolve these problems:\n")
problem.resolve_constraints()
print("\nAfter optimization\n")
print(problem.constraints_text_summary())
