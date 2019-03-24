"""Example of DNA Chisel user-defined specification: Gen9's 9mer score.

Once upon the time there was a DNA synthesis company called Gen9, who accepted
or refused orders based on sequence features.

One of them was the "9mer score", defined as 9 * N / len(sequence) where N is
the number of non-unique kmers in the sequence.

The score can be interpreted as the density of nucleotides which are part of
a 9mer also present somewhere else in the sequence. These should be avoided
as they impede the synthesis (certainly at the oligo assembly stage). Too much
9mers are a... nightmare (ah ah ah pun).

In this script we show how to create a new objective ``MinimizeNinemersScore``
to minimize this score, and we use it to optimize a coding sequence.
"""

from collections import Counter
from dnachisel import (EnforceTranslation, Specification, SpecEvaluation,
                       reverse_translate, random_protein_sequence, Location,
                       DnaOptimizationProblem)


class MinimizeNinemersScore(Specification):
    """Minimize Gen9's "no-9mers" score."""

    def evaluate(self, problem):
        """Return Gen9's ninemer score for the problem' sequence"""
        sequence = problem.sequence
        all_9mers = [sequence[i:i + 9] for i in range(len(sequence) - 9)]
        number_of_non_unique_9mers = sum([
            count
            for ninemer, count in Counter(all_9mers).items()
            if count > 1
        ])
        score = - (9.0 * number_of_non_unique_9mers) / len(sequence)
        return SpecEvaluation(
            self, problem,
            score=score,
            locations=[Location(0, len(sequence))],
            message="Score: %.02f (%d non-unique ninemers)" % (
                score, number_of_non_unique_9mers
            )
        )

    def __str__(self):
        """String representation."""
        return "MinimizeNinemersScore"


sequence = reverse_translate(random_protein_sequence(300))
problem = DnaOptimizationProblem(
    sequence=sequence,
    constraints=[EnforceTranslation()],
    objectives=[MinimizeNinemersScore()]
)

print ("\n=== Status before optimization ===")
print (problem.objectives_text_summary())

problem.optimize()

print ("\n=== Status after optimization ===")
print (problem.objectives_text_summary())
print (problem.constraints_text_summary(failed_only=True))
