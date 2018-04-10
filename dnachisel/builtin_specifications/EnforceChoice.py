"""Implement EnforceSequence (DO NOT USE YET: Work in progress, stabilizing)"""

# TODO: factorize with self.sequence ?

import numpy as np

from ..Specification import Specification
from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from ..SequencePattern import DnaNotationPattern, enzyme_pattern
from dnachisel.Location import Location
from dnachisel.biotools import (group_nearby_indices,
                                reverse_complement,
                                IUPAC_NOTATION)


class EnforceChoice(Specification):
    """Constr

    Parameters
    ----------

    sequence
      An ATGC string representing the wanted sequence, possibly degenerated,
      for instance ATTCGCGTYTTKWNAA

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)`` or simply ``(10, 45, 1)``

    """
    localization_interval_length = 6  # used when optimizing
    best_possible_score = 0
    enforced_by_nucleotide_restrictions = True

    def __init__(self, choices=None, enzymes=None, location=None, boost=1.0):
        """Initialize."""
        if enzymes is not None:
            choices = [enzyme_pattern(e) for e in enzymes]
        if choices is None:
            raise ValueError('`choices` or `enzymes` should not be None.')
        choices = [
             choice if isinstance(choice, DnaNotationPattern)
             else DnaNotationPattern(choice)
             for choice in choices
        ]
        choices = [
            variant
            for choice in choices
            for variant in choice.all_variants()
        ]
        self.choices = choices
        if isinstance(location, tuple):
            location = Location.from_tuple(location, default_strand=+1)
        self.location = location
        self.boost = boost

    def initialize_on_problem(self, problem, role='constraint'):
        """Find out what sequence it is that we are supposed to conserve."""
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            result = self.copy_with_changes(location=location)
        else:
            result = self
        if not all([len(c) == len(result.location) for c in result.choices]):
            raise ValueError("All sequence choices should have the same "
                             "length as the region on which the spec is "
                             "applied.")
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

    def localized(self, location, problem=None):
        """Localize the spec to the overlap of its location and the new."""
        return self

    def restrict_nucleotides(self, sequence, location=None):
        """When localizing, forbid any nucleotide but the one already there."""
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
