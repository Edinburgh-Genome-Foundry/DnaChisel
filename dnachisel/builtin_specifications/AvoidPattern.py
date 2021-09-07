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
      ``Location(10, 45, 1)``. For patterns which are not palindromic,
      the strand matters! Use +1 for eliminating the pattern on the +1 strand
      only, -1 for eliminating the pattern on the -1 strand, and 0 for
      eliminating the pattern on both strands. Default ``None`` enforces on
      the whole sequence.

    strand
      Alternative way to set the strand, meant to be used in two cases only:
      (1) in a Genbank annotation by setting ``strand=both`` to indicate that
      the pattern should be avoided on both strands (otherwise, only the
      feature's strand will be considered).
      (2) if you want to create a specification without preset location, but
      with a set strand: ``AvoidPattern('BsmBI_site', strand=1)``.
      The default 'from_location' uses the strand specified in ``location``,
      or if that is ``None``, it sets both strands.
    """

    best_possible_score = 0
    priority = 1
    shorthand_name = "no"  # will appear as, for instance, @no(BsmBI_site)

    def __init__(self, pattern=None, location=None, strand="from_location", boost=1.0):
        """Initialize."""
        if isinstance(pattern, str):
            pattern = SequencePattern.from_string(pattern)
        self.pattern = pattern
        self.location = Location.from_data(location)

        if strand == "from_location":
            if self.location is None:
                self.strand = 0
            else:
                self.strand = self.location.strand
        elif strand == "both":
            self.strand = 0
        elif strand in [-1, 0, 1]:
            self.strand = strand
        else:
            raise ValueError("unknown strand: %s" % strand)

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

    def breach_label(self):
        if self.pattern.name is not None:
            return str(self.pattern.name)
        else:
            return str(self.pattern)

    def initialized_on_problem(self, problem, role="constraint"):
        copy_of_constraint = self._copy_with_full_span_if_no_location(problem)
        copy_of_constraint.location.strand = self.strand
        return copy_of_constraint

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the pattern to the given location. Taking into account the
        specification's own location, and the size of the pattern."""
        if self.location.overlap_region(location) is None:
            return None
        if self.pattern.size is None:
            return self
        extended_location = location.extended(
            self.pattern.size - 1, right=with_righthand
        )
        new_location = self.location.overlap_region(extended_location)
        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        return [("pattern", str(self.pattern))]
