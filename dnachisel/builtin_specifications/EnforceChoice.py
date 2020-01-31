"""Implement EnforceSequence (DO NOT USE YET: Work in progress, stabilizing)"""

# TODO: factorize with self.sequence ?

from ..SequencePattern import SequencePattern
from ..Location import Location
from ..biotools import reverse_complement
from ..Specification import Specification, SpecEvaluation



class EnforceChoice(Specification):
    """The sequence at the given location must be one of several alternatives.

    Shorthand for annotations: "choice".

    Parameters
    ----------
    choices
      List of same-length ATGC sequences, e.g. ["ATT", "ATC", "ATA"].

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)`` or simply ``(10, 45, 1)``. the location length
      must match that of the sequence in the list of ``choices``

    """

    localization_interval_length = 6  # used when optimizing
    best_possible_score = 0
    enforced_by_nucleotide_restrictions = True
    shorthand_name = "choice"

    def __init__(self, choices=None, location=None, boost=1.0):
        """Initialize."""
        if isinstance(choices, str) and "|" in choices:
            choices = choices.split("|")
        choices = [
            SequencePattern.from_string(c) if isinstance(c, str) else c
            for c in choices
        ]
        # PRECOMPUTE ALL VARIANTS
        choices = [
            variant for choice in choices for variant in choice.all_variants()
        ]
        self.choices = choices
        self.location = Location.from_data(location)
        self.boost = boost

    def initialized_on_problem(self, problem, role="constraint"):
        """Find out what sequence it is that we are supposed to conserve."""
        result = self._copy_with_full_span_if_no_location(problem)
        if not all([len(c) == len(result.location) for c in result.choices]):
            raise ValueError(
                "All sequence choices should have the same length as the "
                "region on which the spec is applied."
            )
        return result

    def evaluate(self, problem):
        """Return a score equal to -number_of modifications.

        Locations are "binned" modifications regions. Each bin has a length
        in nucleotides equal to ``localization_interval_length`.`
        """
        sequence = self.location.extract_sequence(problem.sequence)
        score = 0 if (sequence in self.choices) else -1
        locations = [] if (score == 0) else [self.location]
        return SpecEvaluation(self, problem, score=score, locations=locations)

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the spec to the overlap of its location and the new."""
        return self

    def restrict_nucleotides(self, sequence, location=None):
        """As a constraint, put the choices in the mutation space."""

        if self.location.strand != -1:
            choices = set(self.choices)
        else:
            choices = set([reverse_complement(c) for c in self.choices])
        return [((self.location.start, self.location.end), choices)]

    def __repr__(self):
        """Represent."""
        return "EnforceChoice(%s)" % str(self.location)

    def __str__(self):
        """Represent."""
        return "EnforceChoice(%s)" % str(self.location)

