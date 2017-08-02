"""Implement AvoidPattern"""

from ..Specification import PatternSpecification, VoidSpecification
from ..SpecEvaluation import SpecEvaluation

from dnachisel.Location import Location


class EnforcePattern(PatternSpecification):
    """Enforce a number of occurences of the given pattern in the sequence.

    Parameters
    ----------
    pattern
      A SequencePattern or DnaNotationPattern

    enzyme
      Enzyme name, can be provided instead of a pattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``

    """

    best_possible_score = 0
    shrink_when_localized = False

    def __init__(self, pattern=None, enzyme=None,
                 location=None, occurences=1, boost=1.0):
        """Initialize."""
        PatternSpecification.__init__(self, pattern=pattern, location=location,
                                      enzyme=enzyme)
        self.occurences = occurences
        self.boost = boost

    def evaluate(self, problem):
        """Score the difference between expected and observed n_occurences."""
        location = (self.location if (self.location is not None) else
                    Location(0, len(problem.sequence)))
        locations = self.pattern.find_matches(problem.sequence, location)
        score = -abs(len(locations) - self.occurences)

        if score == 0:
            message = "Passed. Pattern found at positions %s" % locations
        else:
            if self.occurences == 0:
                message = "Failed. Pattern not found."
            else:
                message = ("Failed. Pattern found %d times instead of %d"
                           " wanted at positions %s") % (len(locations),
                                                         self.occurences,
                                                         location)
        return SpecEvaluation(
            self, problem, score, message=message,
            locations=None if location is None else [location],
        )

    def __repr__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.location)

    def __str__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.location)
