"""Implement EnforceTerminalGCContent"""

from .TerminalSpecification import TerminalSpecification
from dnachisel.biotools import gc_content

class EnforceTerminalGCContent(TerminalSpecification):
    """Enforce bounds for the GC content at the sequence's terminal ends.

    Parameters
    ----------
    window_size
      Size in basepair of the two terminal ends to consider

    mini
      A float between 0 and 1, minimal proportion of GC that the ends should
      contain

    maxi
      Float between 0 and 1, maximal proportion of GC that the ends should
      contain

    boost
      Multiplicatory factor applied to this specification.
    """

    def __init__(self, window_size, mini=0, maxi=1, boost=1.0):
        self.mini = mini
        self.maxi = maxi
        self.window_size = window_size
        self.boost = boost

    def evaluate_end(self, sequence):
        gc = gc_content(sequence)
        return -(max(0, self.mini - gc) + max(0, gc - self.maxi))

    def __repr__(self):
        return "Terminal(%.02f < gc < %.02f, window: %d)" % \
            (self.mini, self.maxi, self.window_size)
