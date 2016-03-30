"""Minimization of non-unique ninemers in a coding sequence.

This example shows how to create 
"""

from collections import Counter
from dnachisel.objectives import ObjectiveEvaluation, Objective
from dnachisel import *


class MinimizeNinemersScore(Objective):
    """
    """
    def evaluate(self, canvas):
        seq = canvas.sequence
        all_9mers = [seq[i:i + 9] for i in range(len(seq) - 9)]
        number_non_unique_9mers = sum([
            count
            for ninemer, count in Counter(all_9mers).items()
            if count>1
        ])
        score = - 9.0 * number_non_unique_9mers / len(seq)
        return ObjectiveEvaluation(
            self, canvas, score=score, windows=[[0, len(seq)]],
            message="Score: %.02f (%d non-unique ninemers)" % (
                score, number_non_unique_9mers
            )
        )
    def __str__(self):
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
