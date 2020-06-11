"""Implementation of AvoidHairpins."""

from ..Specification import Specification, SpecEvaluation
from ..biotools import reverse_complement, group_nearby_segments
from ..Location import Location


class AvoidHairpins(Specification):
    """Avoid Hairpin patterns as defined by the IDT guidelines.

    A hairpin is defined by a sequence segment which has a reverse complement
    "nearby" in a given window.

    Parameters
    ----------
    stem_size
      Size of the stem of a hairpin, i.e. the length of the sequence which
      should have a reverse complement nearby to be considered a hairpin.

    hairpin_window
      The window in which the stem's reverse complement should be searched for.

    boost
      Multiplicative factor, importance of this objective in a multi-objective
      optimization.
    """

    best_possible_score = 0

    def __init__(
        self, stem_size=20, hairpin_window=200, location=None, boost=1.0
    ):
        """Initialize."""
        self.stem_size = stem_size
        self.hairpin_window = hairpin_window
        self.location = Location.from_data(location)
        self.boost = boost

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Return the score (-number_of_hairpins) and hairpins locations."""
        sequence = self.location.extract_sequence(problem.sequence)
        reverse = reverse_complement(sequence)
        locations = []
        for i in range(len(sequence) - self.stem_size):
            word = sequence[i : i + self.stem_size]
            rest = reverse[-(i + self.hairpin_window) : -(i + self.stem_size)]
            if word in rest:
                index = rest.index(word)
                locations.append((i, i + self.hairpin_window - index - 1))
        score = -len(locations)
        locations = group_nearby_segments(locations, max_start_spread=10)
        locations = sorted([Location(l[0][0], l[-1][1]) for l in locations])

        return SpecEvaluation(self, problem, score, locations=locations)

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the spec, make sure no neighbouring hairpin is created."""
        new_location = self.location.overlap_region(location)
        if new_location is None:
            return None
        # VoidSpecification(parent_specification=self)
        else:
            new_location.start = max(
                self.location.start, new_location.start - self.hairpin_window
            )
            if with_righthand:
                new_location.end = min(
                    self.location.end, new_location.end + self.hairpin_window
                )
            return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        return [
            ("stem_size", str(self.stem_size)),
            ("hairpin_window", str(self.hairpin_window)),
        ]

    def short_label(self):
        stem = self.stem_size
        inside = self.hairpin_window - 2 * self.stem_size
        return "No %d-%d-%dbp hairpin" % (stem, inside, stem)
    
    def breach_label(self):
        stem = self.stem_size
        inside = self.hairpin_window - 2 * self.stem_size
        return "%d-%d-%dbp hairpin" % (stem, inside, stem)
