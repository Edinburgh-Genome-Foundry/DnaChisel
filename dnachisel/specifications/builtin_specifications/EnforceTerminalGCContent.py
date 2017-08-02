"""Implement EnforceTerminalGCContent"""

from ..Specification import TerminalSpecification
from dnachisel.biotools import gc_content

class EnforceTerminalGCContent(TerminalSpecification):
    """Enforce bounds for the GC content at the sequence's terminal ends.

    Parameters
    ----------
    window_size
      Size in basepair of the two terminal ends to consider

    gc_min
      A float between 0 and 1, minimal proportion of GC that the ends should
      contain

    gc_max
      Float between 0 and 1, maximal proportion of GC that the ends should
      contain

    boost
      Multiplicatory factor applied to this specification.
    """

    def __init__(self, window_size, gc_min=0, gc_max=1, boost=1.0):
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.window_size = window_size
        self.boost = boost

    def evaluate_end(self, sequence):
        return (self.gc_min < gc_content(sequence) < self.gc_max)

    def __repr__(self):
        return "Terminal(%.02f < gc < %.02f, window: %d)" % \
            (self.gc_min, self.gc_max, self.window_size)
