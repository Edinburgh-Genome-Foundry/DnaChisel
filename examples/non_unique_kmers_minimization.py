"""Minimization of non-unique 9-mers in a coding sequence.

This example shows how to create you own ``Objectives`` in DnaChisel.

The DNA provider Gen9 indicates that repeated 9-mers (i.e. 9-base-pair segments
of a sequence which have at least another identical twin somewhere else in the
sequence) create problems for synthesis (maybe due to misannealing ?) and
should be avoided.
So too much 9-mers, and synthesis becomes a 9-mer (nightmare), ah ah ah !
Anyways, Gen9 wants us to get read of these repeated 9-mers and has created
the following score that we should minimize:

    score = 9.0 * number_of_non_unique_9mers / sequence_length

In this script we show how to create a new objective ``MinimizeNinemersScore``
to minimize this score, and we use it to optimize a coding sequence.
"""

from collections import Counter
from dnachisel import (EnforceTranslation, ObjectiveEvaluation, Objective,
                       reverse_translate, random_protein_sequence,
                       DnaOptimizationProblem)


class MinimizeNinemersScore(Objective):
    """Minimize Gen9's "no-9mers" score."""

    def evaluate(self, canvas):
        """Return Gen9's ninemer score for the canvas' sequence"""
        sequence = canvas.sequence
        all_9mers = [sequence[i:i + 9] for i in range(len(sequence) - 9)]
        number_of_non_unique_9mers = sum([
            count
            for ninemer, count in Counter(all_9mers).items()
            if count > 1
        ])
        score = - 9.0 * number_of_non_unique_9mers / len(sequence)
        return ObjectiveEvaluation(
            self, canvas,
            score=score,
            windows=[[0, len(sequence)]],
            message="Score: %.02f (%d non-unique ninemers)" % (
                score, number_of_non_unique_9mers
            )
        )

    def __str__(self):
        """String representation."""
        return "MinimizeNinemersScore"


sequence = reverse_translate(random_protein_sequence(300))
canvas = DnaOptimizationProblem(
    sequence=sequence,
    constraints=[EnforceTranslation([0, len(sequence)],
                                    sequence=sequence)],
    objectives=[MinimizeNinemersScore()]
)

print ("\n=== Status before optimization ===")
print (canvas.objectives_summary())

canvas.maximize_all_objectives_one_by_one()

print ("\n=== Status after optimization ===")
print (canvas.objectives_summary())
print (canvas.constraints_summary(failed_only=True))
