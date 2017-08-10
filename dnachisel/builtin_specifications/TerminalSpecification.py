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
        L = len(sequence)
        wsize = self.window_size
        locations = [
            location
            for location in [Location(0, wsize), Location(L - wsize, L)]
            if not self.evaluate_end(location.extract_sequence(sequence))
        ]

        if locations == []:
            message = "Passed (no breach at the ends)"
        else:
            message = "Failed: breaches at ends %s" % str(locations)

        return SpecEvaluation(self, problem, score=-len(locations),
                              locations=locations, message=message)
