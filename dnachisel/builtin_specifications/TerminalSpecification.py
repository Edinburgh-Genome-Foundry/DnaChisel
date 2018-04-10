"""Implements TerminalSpecification."""

from ..Specification import Specification
from ..SpecEvaluation import SpecEvaluation
from ..Location import Location

class TerminalSpecification(Specification):
    """Specifications that apply in the same way to both ends of the sequence.

    These are particularly useful for modeling constraints from providers
    who have terminal-ends constraints.

    Subclasses of these specifications should have a `location_size` and a
    `evaluate_end` method"""

    def evaluate(self, problem):
        """Apply method ``evaluate_end`` to both sides and compile results."""
        sequence = problem.sequence
        ends_evaluations = [
            self.evaluate_end(location.extract_sequence(sequence))
            for location in self.ends_locations
        ]

        locations = [
            location
            for (score, location) in zip(ends_evaluations, self.ends_locations)
            if score < 0
        ]

        score = sum(ends_evaluations)
        if locations == []:
            message = "Passed (no breach at the ends)"
        else:
            message = "Failed: breaches at ends %s" % str(locations)

        return SpecEvaluation(self, problem, score=score, locations=locations,
                              message=message)

    def initialize_on_problem(self, problem, role):
        """Find out what sequence it is that we are supposed to conserve."""
        if not hasattr(self, 'ends_locations') or self.ends_locations is None:
            L = len(problem.sequence)
            wsize = self.window_size
            ends_locations = [Location(0, wsize), Location(L - wsize, L)]
            return self.copy_with_changes(ends_locations=ends_locations)
        else:
            return self

    def localized(self, location, problem=None):
        new_end_locations = []
        for end_location in self.ends_locations:
            overlap = location.overlap_region(end_location)
            if overlap is not None:
                new_end_locations.append(end_location)

        if len(new_end_locations) == 2:
            return self
        else:
            return self.copy_with_changes(ends_locations=new_end_locations)
