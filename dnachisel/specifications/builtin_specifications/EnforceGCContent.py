"""Implement EnforceGCContent."""

import numpy as np

from ..Specification import Specification, VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from dnachisel.biotools import gc_content
from dnachisel.Location import Location


class EnforceGCContent(Specification):
    """Specification on the local or global proportion of G/C nucleotides.

    Examples
    --------
    >>> # Enforce global GC content between 40 and 70 percent.
    >>> Specification = GCContentSpecification(0.4, 0.7)
    >>> # Enforce 30-80 percent local GC content over 50-nucleotides windows
    >>> Specification = GCContentSpecification(0.3, 0.8, gc_window=50)


    Parameters
    ----------
    gc_min
      Minimal proportion of G-C (e.g. ``0.35``)

    gc_max
      Maximal proportion of G-C (e.g. ``0.75``)

    gc_window
      Length of the sliding window, in nucleotides, for local GC content.
      If not provided, the global GC content of the whole sequence is
      considered

    location
      Location objet indicating that the Specification only applies to a
      subsegment of the sequence. Make sure it is bigger than ``gc_window``
      if both parameters are provided

    """

    best_possible_score = 0.0

    def __init__(self, gc_min=0, gc_max=1.0, gc_objective=None,
                 gc_window=None, location=None, boost=1.0):
        """Initialize."""
        if gc_objective is not None:
            gc_min = gc_max = gc_objective
        self.gc_objective = gc_objective
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.gc_window = gc_window
        self.location = location
        self.boost = boost

    def evaluate(self, problem):
        """Return the sum of breaches extent for all windowed breaches."""
        location = (self.location if self.location is not None
                    else Location(0, len(problem.sequence)))
        wstart, wend = location.start, location.end
        sequence = location.extract_sequence(problem.sequence)
        gc = gc_content(sequence, self.gc_window)
        breaches = (np.maximum(0, self.gc_min - gc) +
                    np.maximum(0, gc - self.gc_max))
        score = - (breaches.sum())
        breaches_starts = (breaches > 0).nonzero()[0]

        if len(breaches_starts) == 0:
            breaches_locations = []
        elif len(breaches_starts) == 1:
            if self.gc_window is not None:
                start = breaches_starts[0]
                breaches_locations = [[start, start + self.gc_window]]
            else:
                breaches_locations = [[wstart, wend]]
        else:
            breaches_locations = []
            current_start = breaches_starts[0]
            last_end = current_start + self.gc_window
            for i in breaches_starts[1:]:
                if (i > last_end + self.gc_window):
                    breaches_locations.append([
                        wstart + current_start, wstart + last_end]
                    )
                    current_start = i
                    last_end = i + self.gc_window

                else:
                    last_end = i + self.gc_window
            breaches_locations.append(
                [wstart + current_start, wstart + last_end])

        if breaches_locations == []:
            message = "Passed !"
        else:
            breaches_locations = [Location(*loc) for loc in breaches_locations]
            message = ("Failed: GC content out of bound on segments " +
                       ", ".join([str(l) for l in breaches_locations]))
        return SpecEvaluation(self, problem, score,
                              locations=breaches_locations,
                              message=message)

    def localized(self, location):
        """Localize the GC content evaluation.

        For a location, the GC content evaluation will be restricted
        to [start - gc_window, end + gc_window]
        """
        if self.location is not None:
            if self.gc_window is None:
                return self
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return VoidSpecification(parent_specification=self)
            else:
                extension = 0 if self.gc_window is None else self.gc_window - 1
                extended_location = location.extended(extension)

                new_location = self.location.overlap_region(extended_location)
        else:
            if self.gc_window is not None:
                new_location = location.extended(self.gc_window + 1)
            else:
                new_location = None
        return self.copy_with_changes(location=new_location)

    def __repr__(self):
        """Represent."""
        return (
            "EnforceGCContent(min %.02f, max %.02f, gc_win %s, location %s)" %
            (self.gc_min, self.gc_max,
             "global" if (self.gc_window is None) else self.gc_window,
             self.location))

    def __str__(self):
        """Represent."""
        return (
            "EnforceGCContent(min %.02f, max %.02f, gc_win %s, location %s)" %
            (self.gc_min, self.gc_max,
             "global" if (self.gc_window is None) else self.gc_window,
             self.location))
