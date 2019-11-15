"""Implements the SequencePattern, DnaNotationPattern classes.

These classes are responsible for looking for a pattern in a sequence
(including overlapping patterns !), separating patterns with fixed size
and patterns with maximal size (from problem localization purposes).

The module also implements functions to specify common DNA patterns:
homopolymers, repeats, enzymatic restriction sites.


"""

import re
from ..biotools import reverse_complement
from ..Location import Location


class SequencePattern:
    """Pattern/ that will be looked for in a DNA sequence.

    Use this class for matching regular expression patterns, and
    DnaNotationPattern for matching explicit sequences or sequences using Ns
    etc.

    Examples
    --------
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

    registered_string_pattern_classes = []

    def __init__(
        self,
        expression,
        size=None,
        name=None,
        in_both_strands=True,
        lookahead="loop",
    ):
        if size is None:
            size = len(expression)
        self.expression = expression
        self.lookahead = lookahead
        if lookahead == "re":
            expression = "(?=(%s))" % expression
        if "(" not in expression:
            expression = "(%s)" % expression
        self.lookahead_expression = expression
        self.compiled_expression = re.compile(self.lookahead_expression)
        self.size = size
        self.name = name
        self.in_both_strands = in_both_strands

    def find_matches(self, sequence, location=None):
        """Return the locations where the sequence matches the expression.

        Parameters
        ----------

        sequence
          A string of "ATGC..."

        location
          Location indicating a segment to which to restrict
          the search. Only patterns entirely included in the segment will be
          returned

        Returns
        -------

        matches
          A list of the locations of matches, of the form
          ``[(start1, end1), (start2, end2),...]``.

        """

        if location is not None:
            subsequence = location.extract_sequence(sequence)
            return [
                (loc + location.start)
                if (location.strand != -1)
                else Location(
                    location.end - loc.end, location.end - loc.start, strand=-1
                )
                for loc in self.find_matches(subsequence)
            ]
        matches = self.find_matches_in_string(sequence)

        if self.in_both_strands:
            reverse = reverse_complement(sequence)
            L = len(sequence)
            matches += [
                (L - end, L - start, -1)
                for (start, end, strand) in self.find_matches_in_string(reverse)
            ]

        return [Location(start, end, strand) for start, end, strand in matches]

    def find_matches_in_string(self, sequence):
        if self.lookahead == "loop":
            matches = []
            position = 0
            while True:
                result = re.search(self.compiled_expression, sequence)
                if result is None:
                    return matches
                start, end = result.start(), result.end()
                matches.append((start + position, end + position, 1))
                sequence = sequence[start + 1 :]
                position += start + 1
        else:
            return [
                (match.start(), match.start() + len(match.groups()[0]), 1)
                for match in re.finditer(self.compiled_expression, sequence)
            ]

    def __str__(self):
        return self.expression + (
            "" if self.name is None else " (%s)" % self.name
        )

    @classmethod
    def from_string(cls, string):
        for myclass in cls.registered_string_pattern_classes:
            pattern = myclass.from_string(string)
            if pattern is not None:
                return pattern
        return SequencePattern(string)
