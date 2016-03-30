import re
from biotools import dna_pattern_to_regexpr
from Bio.Restriction.Restriction_Dictionary import rest_dict

class DNAPattern:
    """Pattern that will be looked for in a DNA sequence.

    Example
    -------
    >>> pattern =
    >>> constraint = NoPatternConstraint(DNAPattern()

    Parameters
    ----------

    pattern
      A string representing a DNA sequence (like "ATTGC"), an DNA sequence in

    size
      Size of the pattern, in number of characters (if none provided, the size
      of the ``pattern`` string is used).
      The ``size`` is used to determine the size of windows when performing
      local optimization and constraint solving.
      It can be important to provide the size when the
      ``pattern`` string provided represents a complex regular expression whose
      maximal matching size cannot be easily evaluated.


    name
      Name of the pattern (will be displayed e.g. when the pattern is printed)

    """

    def __init__(self, pattern, size=None, name=None):
        if size is None:
            size = len(pattern)
        self.pattern = pattern
        self.regexpr = dna_pattern_to_regexpr(pattern)
        self.size = size
        self.name = name

    def find_matches(self, sequence, window=None):
        """Return the locations where the sequence matches the pattern.

        Parameters
        ----------

        sequence
          A string of "ATGC..."

        window
          Pair ``(start, end)`` indicating a segment to which to restrict
          the search. Only patterns entirely included in the segment will be
          returned

        Returns
        -------

        matches
          A list of the locations of matches, of the form
          ``[(start1, end1), (start2, end2),...]``.

        """


        if window is not None:
            wstart, wend = window
            subsequence = sequence[wstart:wend]
            return [
                (wstart + start, wstart + end)
                for (start, end) in self.find_matches(subsequence)
            ]
        return [
            (match.start(), match.end())
            for match in re.finditer(self.regexpr, sequence)
        ]

    def __str__(self):
        return self.pattern + ("" if self.name is None else
                               " (%s)" % self.name)


def homopolymer_pattern(letter, number):
    """Return a DNAPattern with the sequence of a homopolymer.

    Examples
    --------

    >>> homopolymer("A", 6) # returns DNAPattern("AAAAAA")

    """
    return DNAPattern(number * letter)


def enzyme_pattern(enzyme_name):
    """Return a DNAPattern with the sequence of a homopolymer.

    Examples
    --------

    >>> enzyme_pattern("BsaI") # returns DNAPattern("GGTCTC")

    """
    enzyme_site = rest_dict[enzyme_name]["site"]
    return DNAPattern(enzyme_site, name=enzyme_name)
