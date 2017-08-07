"""Implement SequenceLengthBounds."""

from ..Specification import Specification
from ..SpecEvaluation import SpecEvaluation


class SequenceLengthBounds(Specification):
    """Checks that the sequence length is between bounds.

    Quite an uncommon specification as it can't really be solved or optimized.
    But practical as part of a list of constraints to verify.

    Parameters
    ----------

    min_length
      Minimal allowed sequence length in nucleotides

    max_length
      Maximal allowed sequence length in nucleotides. None means no bound.
    """
    best_possible_score = 0

    def __init__(self, min_length=0, max_length=None):
        self.min_length = min_length
        self.max_length = max_length

    def evaluate(self, problem):
        """Return 0 if the sequence length is between the bounds, else -1"""
        L, mini, maxi = len(problem.sequence), self.min_length, self.max_length
        if maxi is None:
            score = (L >= mini)
        else:
            score = (mini <= L <= maxi)
        return SpecEvaluation(self, problem, score - 1)

    def __repr__(self):
        return "Length(%d < L < %d)" % (self.min_length, self.max_length)
