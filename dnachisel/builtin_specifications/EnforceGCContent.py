"""Implement EnforceGCContent."""

import numpy as np
import re

from ..biotools import gc_content, group_nearby_segments
from ..Location import Location
from ..Specification import Specification, SpecEvaluation



class EnforceGCContent(Specification):
    """Specification on the local or global proportion of G/C nucleotides.

    Shorthand for annotations: "gc".

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

    target
      Target proportion of GC (e.g. ``0.4``), which can be used instead of
      ``mini`` and ``maxi`` when using the specification as an optimization
      objective.

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
    locations_span = 50  # The resolution will use locations size
    shorthand_name = "gc"

    def __init__(
        self,
        mini=0,
        maxi=1.0,
        target=None,
        window=None,
        location=None,
        boost=1.0,
    ):
        """Initialize."""
        if isinstance(mini, str):
            mini, maxi, target, _window = self.string_to_parameters(mini)
            if _window is not None:
                window = _window
        if target is not None:
            mini = maxi = target
        self.target = target
        self.mini = mini
        self.maxi = maxi
        self.window = window
        location = Location.from_data(location)
        if location is not None and (location.strand == -1):
            location = Location(location.start, location.end, 1)
        self.location = location
        self.boost = boost

    @staticmethod
    def string_to_parameters(pattern):
        target = None
        gc_pattern = r"(\d+)(|-(\d+))%($|/(\d+)bp)"
        match = re.match(gc_pattern, pattern)
        mini, _, maxi, _, window = match.groups()
        if maxi is None:
            target = mini
            mini = None
        mini, maxi, target = [
            None if e is None else int(e) / 100.0 for e in (mini, maxi, target)
        ]
        window = None if window is None else int(window)
        return mini, maxi, target, window

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Return the sum of breaches extent for all windowed breaches."""
        wstart, wend = self.location.start, self.location.end
        sequence = self.location.extract_sequence(problem.sequence)
        gc = gc_content(sequence, window_size=self.window)
        breaches = np.maximum(0, self.mini - gc) + np.maximum(
            0, gc - self.maxi
        )
        score = -breaches.sum()
        breaches_starts = wstart + (breaches > 0).nonzero()[0]

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
                segments, max_start_spread=max(1, self.locations_span)
            )
            breaches_locations = [
                (group[0][0], group[-1][-1]) for group in groups
            ]

        if breaches_locations == []:
            message = "Passed !"
        else:
            breaches_locations = [Location(*loc) for loc in breaches_locations]
            message = "Out of bound on segments " + ", ".join(
                [str(l) for l in breaches_locations]
            )
        return SpecEvaluation(
            self, problem, score, locations=breaches_locations, message=message
        )

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the GC content evaluation.

        For a location, the GC content evaluation will be restricted
        to [start - window, end + window]
        """
        # if self.location is not None:
        if self.window is None:
            # No specific location policy at the moment for non-windowed GC
            # content. The same specification is returned.
            # NOTE: this could be refined by computing
            # how much the local bounds should be in order for the global
            # sequence not get out of bound. This would be faster but much
            # more complicated.
            return self
        new_location = self.location.overlap_region(location)
        if new_location is None:
            return None
        else:
            extension = 0 if self.window is None else self.window - 1
            extended_location = location.extended(
                extension, right=with_righthand
            )
            new_location = self.location.overlap_region(extended_location)
        return self.copy_with_changes(location=new_location)

    def label_parameters(self):
        show_mini = self.mini is not None and (self.target is None)
        show_maxi = self.maxi is not None and (self.target is None)
        show_target = not (show_mini or show_maxi)
        show_window = self.window is not None

        return (
            show_mini * [("mini", "%.2f" % self.mini)]
            + show_maxi * [("maxi", "%.2f" % self.maxi)]
            + show_target
            * [("target", ("%.2f" % self.target) if self.target else "")]
            + show_window
            * [("window", ("%d" % self.window) if self.window else "")]
        )

    def short_label(self):
        if self.target is not None:
            result = "%d%% GC" % (np.round(100 * self.target))
        else:
            result = "%d-%d%% GC" % (
                np.round(100 * self.mini),
                np.round(100 * self.maxi),
            )
        if self.window is not None:
            result += "/%dbp" % self.window
        return result
    
    def breach_label(self):
        if self.target is not None:
            result = "GC != %d%% " % (np.round(100 * self.target))
        else:
            result = "GC outside %d-%d%%" % (
                np.round(100 * self.mini),
                np.round(100 * self.maxi),
            )
        if self.window is not None:
            result += "/%dbp" % self.window
        return result
