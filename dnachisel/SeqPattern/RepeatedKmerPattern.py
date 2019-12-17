import re
from .SequencePattern import SequencePattern

class RepeatedKmerPattern(SequencePattern):
    """Direct repeats like ATT-ATT, ATGC-ATGC-ATGC, etc.

    Shorthand string version: "3x4mer", "5x2mer", etc.

    Examples
    --------

    >>> RepeatedKmerPattern(3, 2) # dimers repeated 3 times

    """

    def __init__(self, n_repeats, kmer_size):
        self.n_repeats = n_repeats
        self.kmer_size = kmer_size
        SequencePattern.__init__(
            self,
            size=kmer_size * n_repeats,
            expression=r"([ATGC]{%d})\1{%d}" % (kmer_size, n_repeats - 1),
            name="%d-repeats %d-mers" % (n_repeats, kmer_size),
            is_palyndromic=True,  # a repeat on a strand is also on the other
            lookahead="loop",
        )

    @staticmethod
    def from_string(string):
        match = re.match(r"(\d+)x(\d+)mer$", string)
        if match is not None:
            n_repeats, kmer_size = match.groups()
            return RepeatedKmerPattern(int(n_repeats), int(kmer_size))

    def __str__(self):
        return "%sx%smer" % (self.n_repeats, self.kmer_size)