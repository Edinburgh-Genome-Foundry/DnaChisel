"""Implement AvoidPattern"""

from ..Specification import PatternSpecification
from ..SpecEvaluation import SpecEvaluation

from dnachisel.Location import Location


class AvoidPattern(PatternSpecification):
    """Enforce that the given pattern is absent in the sequence.

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    enzyme
      Enzyme name, can be provided instead of pattern or dna_pattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``

    """
    best_possible_score = 0

    def evaluate(self, problem):
        """Return score=-number_of_occurences. And patterns locations."""
        locations = self.pattern.find_matches(problem.sequence, self.location)
        score = -len(locations)
        if score == 0:
            message = "Passed. Pattern not found !"
        else:
            message = "Failed. Pattern found at positions %s" % locations
        return SpecEvaluation(
            self, problem, score, locations=locations, message=message
        )

    def __repr__(self):
        """Represent."""
        return "AvoidPattern(%s, %s)" % (self.pattern, self.location)

    def __str__(self):
        """Represent."""
        return "AvoidPattern(%s, %s)" % (self.pattern, self.location)
