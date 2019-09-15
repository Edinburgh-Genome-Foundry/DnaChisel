"""We will build a collection of 20 non-interannealing primers.

DNA Chisel is not originally meant for creating collections of sequences
(frameworks such as D-tailor were written with this purpose in mind), but
it is still possible to create collections of inter-compatible sequences.

Here we create 20 short sequences one after the other, while making sure
that each new primer is "compatible" with the other ones made so far.
Compatibility in this example means low heterodimerization (=the primers
could be used together without undesired annealing between them).

"""

from dnachisel import DnaOptimizationProblem, random_dna_sequence
from dnachisel.builtin_specifications import (
    AvoidHeterodimerization,
    EnforceGCContent,
    AvoidPattern,
)


def create_new_primer(existing_primers):
    """Create a new primer based on the primers created so far"""
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(length=20),
        constraints=[
            AvoidHeterodimerization(existing_primers, tmax=3),
            AvoidPattern("3x3mer"),
            AvoidPattern("4xG"),
        ],
        objectives=[EnforceGCContent(target=0.6)],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence


# MAIN LOOP, WHERE PRIMERS ARE CREATED ONE BY ONE

existing_primers = []
for i in range(20):
    new_primer = create_new_primer(existing_primers)
    existing_primers.append(new_primer)

print("PRIMERS GENERATED: \n\n%s\n" % "\n".join(existing_primers))

# (OPTIONAL VERIFICATION OF THE COLLECTION)

import itertools
import primer3
from dnachisel.biotools import gc_content

max_tm = max(
    primer3.calcHeterodimer(seq1, seq2).tm
    for seq1, seq2 in itertools.combinations(existing_primers, 2)
)
print("Max Tm heterodimerization between any 2 primers: %.2f" % max_tm)

gc_contents = [gc_content(p) for p in existing_primers]
print("GC content range %.2f-%.2f" % (min(gc_contents), max(gc_contents)))

