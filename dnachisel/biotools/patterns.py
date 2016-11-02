import re
from .biotools import (dna_pattern_to_regexpr, is_palyndromic,
                       reverse_complement)
from Bio.Restriction.Restriction_Dictionary import rest_dict

class SequencePattern:
    """Pattern/ that will be looked for in a DNA sequence.

    Use this class for matching regular expression patterns, and
    DnaNotationPattern for matching explicit sequences or sequences using Ns
    etc.

    Example
    -------
    >>> expression = "A[ATGC]{3,}"
    >>> pattern = SequencePattern(expression)
    >>> constraint = AvoidPattern(pattern)

    Parameters
    ----------

    expression
      Any string or regular expression for matching ATGC nucleotides.
      Note that multi-nucleotides symbols such as "N" (for A-T-G-C), or "K"
      are not supported by this class, see DnaNotationPattern instead.

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

    in_both_strands
      Set to True (default) if the pattern should also be looked for on the
      reverse-complement of sequences.

    """

    def __init__(self, expression, size=None, name=None, in_both_strands=True):
        if size is None:
            size = len(expression)
        self.expression = expression
        self.size = size
        self.name = name
        self.in_both_strands = in_both_strands


    def find_matches(self, sequence, window=None):
        """Return the locations where the sequence matches the expression.

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
        matches = [
            (match.start(), match.end())
            for match in re.finditer(self.expression, sequence)
        ]

        if self.in_both_strands:
            reverse = reverse_complement(sequence)
            L = len(sequence)
            matches += [
                (L - match.end(), L - match.start())
                for match in re.finditer(self.expression, reverse)
            ]

        return matches

    def __str__(self):
        return self.expression + ("" if self.name is None else
                                  " (%s)" % self.name)


class DnaNotationPattern(SequencePattern):
    """Class for patterns in plain DNA notation: ATTGCCA, GCNNKTA, etc.

    If the sequence is not palyndromic, the pattern will be looked for in
    both strands of sequences.
    """

    def __init__(self, pattern, name=None):
        SequencePattern.__init__(
            self,
            size=len(pattern),
            expression=dna_pattern_to_regexpr(pattern),
            name=name,
            in_both_strands=not is_palyndromic(pattern)
        )
        self.pattern = pattern

    def __repr__(self):
        return self.pattern + ("" if self.name is None else
                               " (%s)" % self.name)

# DEFINITION OF COMMON PATTERNS


def homopolymer_pattern(letter, number):
    """Return a DnaNotationPattern with the sequence of a homopolymer.

    Examples
    --------

    >>> homopolymer("A", 6) # returns DnaNotationPattern("AAAAAA")

    """
    return DnaNotationPattern(number * letter)


def enzyme_pattern(enzyme_name):
    """Return a DnaNotationPattern with the sequence of a homopolymer.

    Examples
    --------

    >>> pattern = enzyme_pattern("BsaI") # returns DnaNotationPattern("GGTCTC")
    >>> constraint = AvoidPattern(pattern)

    """
    enzyme_site = rest_dict[enzyme_name]["site"]
    return DnaNotationPattern(enzyme_site, name=enzyme_name)


def repeated_kmers(kmer_size, n_repeats):
    """Return a SequencePattern matching all k-mers repeated n times.

    Examples
    --------

    >>> pattern = repeated_kmers(2, 5) # Result will match e.g. ACACACACAC
    >>> constraint = AvoidPattern(pattern)
    """

    return SequencePattern(
        size=kmer_size * n_repeats,
        expression=r"([ATGC]{%d})\1{%d}" % (kmer_size, n_repeats-1),
        name="%d-repeats %d-mers" % (n_repeats, kmer_size),
        in_both_strands=False  # a kmer repeat one strand is also on the other
    )
