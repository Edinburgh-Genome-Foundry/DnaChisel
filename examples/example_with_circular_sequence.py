"""This example demonstrates the use of DNA Chisel to solve a circular DNA
optimization problem:

- The sequence is designed to have a cross-origin BsmBI site that will need
  to be removed, because the location-less specification ``AvoidPattern``
  is interpreted as applying to the full circle.

- The specification ``EnforceGCContent`` is cross-origin since its location is
  1500-2500, and the sequence is ~2000bp long.

"""

import dnachisel as dc

dna_sequence = (
    "CTC"
    + dc.random_dna_sequence(1000)
    + "CGTCTC"
    + dc.random_dna_sequence(1000)
    + "CGT"
)

problem = dc.CircularDnaOptimizationProblem(
    sequence=dna_sequence,
    constraints=[
        dc.AvoidPattern("BsmBI_site"),
        dc.EnforceGCContent(
            mini=0.4, maxi=0.6, location=(1500, 2500), window=50
        ),
        dc.AvoidNonUniqueSegments(
            min_length=9, location=(10, 1000), extended_location=(20, 500)
        )
    ],
    logger=None,
)

print("BEFORE OPTIMIZATION:\n\n", problem.constraints_text_summary())
problem.resolve_constraints()
print("AFTER OPTIMIZATION:\n\n", problem.constraints_text_summary())

