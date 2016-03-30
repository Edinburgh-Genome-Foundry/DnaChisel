"""Minimization of non-unique 9-mers in a coding sequence.

This example shows how to create you own ``Objectives`` in DNAChisel.

The DNA provider Gen9 indicates that repeated 9-mers (i.e. 9-base-pair segments
of a sequence which have at least another identical twin somewhere else in the
sequence) create problems for synthesis (maybe due to misannealing ?) and should
be avoided.
So too much 9-mers, and synthesis becomes a 9-mer (nightmare), ah ah ah !
Anyways, Gen9 wants us to get read of these repeated 9-mers and has created
the following score that we should minimize:

    score = 9.0 * number_of_non_unique_9mers / sequence_length

In this script we show how to create a new objective ``MinimizeNinemersScore``
to minimize this score, and we use it to optimize a coding sequence.
"""

from collections import Counter
from dnachisel.objectives import ObjectiveEvaluation, Objective
from dnachisel import *


class MinimizeNinemersScore(Objective):
    """Minimize Gen9's "no-9mers" score."""

    def evaluate(self, canvas):
        """Return Gen9's ninemer score for the canvas' sequence"""
        seq = canvas.sequence
        all_9mers = [seq[i:i + 9] for i in range(len(seq) - 9)]
        number_of_non_unique_9mers = sum([
            count
            for ninemer, count in Counter(all_9mers).items()
            if count > 1
        ])
        score = - 9.0 * number_of_non_unique_9mers / len(seq)
        return ObjectiveEvaluation(
            self, canvas, score=score, windows=[[0, len(seq)]],
            message="Score: %.02f (%d non-unique ninemers)" % (
                score, number_non_unique_9mers
            )
        )

    def __str__(self):
        """String representation."""
        return "MinimizeNinemersScore"


sequence = reverse_translate(random_protein_sequence(300))
canvas = DNACanvas(
    sequence=sequence,
    constraints=[EnforceTranslationConstraint([0, len(sequence)],
                                              sequence=sequence)],
    objectives=[MinimizeNinemersScore()]
)

print ("\n=== Status before optimization ===")
canvas.print_objectives_summary()

canvas.maximize_all_objectives_one_by_one()

print ("\n=== Status after optimization ===")
canvas.print_objectives_summary()
canvas.print_constraints_summary(failed_only=True)
