"""Implement AvoidPattern"""

from ..SequencePattern import SequencePattern
from ..Location import Location
from ..Specification.Specification import Specification
from ..Specification.SpecEvaluation import SpecEvaluation



class AvoidPattern(Specification):
    """Enforce that the given pattern is absent in the sequence.

    Shorthand for annotations: "no".

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``

    """

    best_possible_score = 0
    shrink_when_localized = True
    priority = 1
    shorthand_name = "no"

    def __init__(self, pattern=None, location=None, boost=1.0):
        """Initialize."""
        if isinstance(pattern, str):
            pattern = SequencePattern.from_string(pattern)
        self.pattern = pattern
        self.location = Location.from_data(location)
        self.boost = boost

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

    def short_label(self):
        if self.pattern.name is not None:
            return "No %s" % self.pattern.name
        else:
            return "No %s" % self.pattern

    def initialized_on_problem(self, problem, role="constraint"):
        return self._copy_with_full_span_if_no_location(problem)

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the pattern to the given location. Taking into account the
        specification's own location, and the size of the pattern."""
        pattern_size = self.pattern.size
        if self.location.overlap_region(location) is None:
            return None
        else:
            if not self.shrink_when_localized:
                return self
            extended_location = location.extended(
                pattern_size - 1, right=with_righthand
            )
            new_location = self.location.overlap_region(extended_location)

        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        return [("pattern", str(self.pattern))]
