import re
from .DnaNotationPattern import DnaNotationPattern


class HomopolymerPattern(DnaNotationPattern):
    """Homopolymer of the form AAAAAAA, TTTTT, etc.

    Shorthand string version: "7xA", "9xC", etc.

    Examples
    --------

    >>> pattern = HomopolymerPattern("A", 6)
    >>> constraint = AvoidPattern(pattern)

    """

    def __init__(self, nucleotide, number):
        self.nucleotide = nucleotide
        self.number = number
        DnaNotationPattern.__init__(self, number * nucleotide)

    @staticmethod
    def from_string(string):
        match = re.match(r"(\d+)x(\S)$", string)
        if match is not None:
            number, nucleotide = match.groups()
            return HomopolymerPattern(nucleotide, int(number))

    def __str__(self):
        return "%sx%s" % (self.number, self.nucleotide)
