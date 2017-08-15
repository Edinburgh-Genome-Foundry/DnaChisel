"""Implement EnforceGCContent."""

import numpy as np

from ..Specification import Specification
from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from ..biotools import group_nearby_segments
from dnachisel.biotools import gc_content
from dnachisel.Location import Location


class EnforceGCContent(Specification):
    """Specification on the local or global proportion of G/C nucleotides.

    Examples
    --------
    >>> # Enforce global GC content between 40 and 70 percent.
    >>> Specification = GCContentSpecification(0.4, 0.7)
    >>> # Enforce 30-80 percent local GC content over 50-nucleotides windows
    >>> Specification = GCContentSpecification(0.3, 0.8, window=50)


    Parameters
    ----------
    mini
      Minimal proportion of G-C (e.g. ``0.35``)

    maxi
      Maximal proportion of G-C (e.g. ``0.75``)

    window
      Length of the sliding window, in nucleotides, for local GC content.
      If not provided, the global GC content of the whole sequence is
      considered

    location
      Location objet indicating that the Specification only applies to a
      subsegment of the sequence. Make sure it is bigger than ``window``
      if both parameters are provided

    """

    best_possible_score = 0
    locations_span = 50  # The resolution will use locations of this size

    def __init__(self, mini=0, maxi=1.0, target=None,
                 window=None, location=None, boost=1.0):
        """Initialize."""
        if target is not None:
            mini = maxi = target
        self.target = target
        self.mini = mini
        self.maxi = maxi
        self.window = window
        self.location = location
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        if self.location is None:
            location = Location(0, len(problem.sequence))
            return self.copy_with_changes(location=location)
        else:
            return self

    def evaluate(self, problem):
        """Return the sum of breaches extent for all windowed breaches."""
        location = (self.location if self.location is not None
                    else Location(0, len(problem.sequence)))
        wstart, wend = location.start, location.end
        sequence = location.extract_sequence(problem.sequence)
        gc = gc_content(sequence, self.window)
        breaches = (np.maximum(0, self.mini - gc) +
                    np.maximum(0, gc - self.maxi))
        score = - (breaches.sum())
        breaches_starts = (breaches > 0).nonzero()[0]

        if len(breaches_starts) == 0:
            breaches_locations = []
        elif len(breaches_starts) == 1:
            if self.window is not None:
                start = breaches_starts[0]
                breaches_locations = [[start, start + self.window]]
            else:
                breaches_locations = [[wstart, wend]]
        else:
            segments = [(bs, bs + self.window) for bs in breaches_starts]
            groups = group_nearby_segments(
                segments,
                max_start_spread=max(1,  self.locations_span - self.window))
            breaches_locations = [
                (group[0][0], group[-1][-1])
                for group in groups
            ]

        if breaches_locations == []:
            message = "Passed !"
        else:
            breaches_locations = [Location(*loc) for loc in breaches_locations]
            message = ("Out of bound on segments " +
                       ", ".join([str(l) for l in breaches_locations]))
        return SpecEvaluation(self, problem, score,
                              locations=breaches_locations,
                              message=message)

    def localized(self, location, with_righthand=True):
        """Localize the GC content evaluation.

        For a location, the GC content evaluation will be restricted
        to [start - window, end + window]
        """
        # if self.location is not None:
        if self.window is None:
            # NOTE: this makes sense, but could be refined by computing
            # how much the local bounds should be in order not to outbound
            return self
        new_location = self.location.overlap_region(location)
        if new_location is None:
            return VoidSpecification(parent_specification=self)
        else:
            extension = 0 if self.window is None else self.window - 1
            extended_location = location.extended(
                extension, right=with_righthand)

            new_location = self.location.overlap_region(extended_location)
        # else:
        #     if self.window is not None:
        #         new_location = location.extended(self.window + 1)
        #     else:
        #         new_location = None
        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        show_mini = self.mini is not None
        show_maxi = self.maxi is not None
        show_target = not (show_mini or show_maxi)
        show_window = self.window is not None

        return (
            show_mini * [('mini', "%.2f" % self.mini)] +
            show_maxi * [('maxi', "%.2f" % self.maxi)] +
            show_target * [('target',
                           ("%.2f" % self.target) if self.target else '')] +
            show_window * [('window', "%d" % self.window)]
        )
