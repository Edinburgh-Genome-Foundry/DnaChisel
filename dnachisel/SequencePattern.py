"""Implements the SequencePattern, DnaNotationPattern classes.

These classes are responsible for looking for a pattern in a sequence
(including overlapping patterns !), separating patterns with fixed size
and patterns with maximal size (from problem localization purposes).

The module also implements functions to specify common DNA patterns:
homopolymers, repeats, enzymatic restriction sites.


"""

import re
import itertools
from .Location import Location
from .biotools import (dna_pattern_to_regexpr, is_palyndromic,
                       reverse_complement, NUCLEOTIDE_TO_REGEXPR,
                       IUPAC_NOTATION)
from Bio.Restriction.Restriction_Dictionary import rest_dict

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

    def __init__(self, expression, size=None, name=None, in_both_strands=True,
                 lookahead='loop'):
        if size is None:
            size = len(expression)
        self.expression = expression
        self.lookahead = lookahead
        if lookahead == 're':
            expression = '(?=(%s))' % expression
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
                (loc + location.start) if (location.strand != -1)
                else Location(location.end - loc.end,
                              location.end - loc.start,
                              strand=-1)
                for loc in self.find_matches(subsequence)
            ]
        matches = self.find_all_re_matches(sequence)

        if self.in_both_strands:
            reverse = reverse_complement(sequence)
            L = len(sequence)
            matches += [
                (L - end, L - start, -1)
                for (start, end, strand) in self.find_all_re_matches(reverse)
            ]

        return [
            Location(start, end, strand)
            for start, end, strand in matches
        ]
    def find_all_re_matches(self, sequence):
        if self.lookahead == 'loop':
            matches = []
            position = 0
            while True:
                result = re.search(self.compiled_expression, sequence)
                if result is None:
                    return matches
                start, end = result.start(), result.end()
                matches.append((start + position, end + position, 1))
                sequence = sequence[start + 1:]
                position += start + 1
        else:
            return [
                (match.start(), match.start() + len(match.groups()[0]), 1)
                for match in re.finditer(self.compiled_expression, sequence)
            ]

    def __str__(self):
        return self.expression + ("" if self.name is None else
                                  " (%s)" % self.name)
    
    @classmethod
    def from_string(cls, string):
        for myclass in cls.registered_string_pattern_classes:
            pattern = myclass.from_string(string)
            if pattern is not None:
                return pattern
        return SequencePattern(string)


class DnaNotationPattern(SequencePattern):
    """Class for patterns in plain DNA notation: ATTGCCA, GCNNKTA, etc.

    If the sequence is not palyndromic, the pattern will be looked for in
    both strands of sequences.
    """

    def __init__(self, sequence, name=None, in_both_strands='auto'):
        """Initialize"""
        if in_both_strands == 'auto':
            in_both_strands = not is_palyndromic(sequence)
        SequencePattern.__init__(
            self,
            size=len(sequence),
            expression=self.dna_sequence_to_regexpr(sequence),
            name=name,
            in_both_strands=in_both_strands
        )
        self.sequence = sequence

    @staticmethod
    def dna_sequence_to_regexpr(sequence):
        """Return a regular expression to find the pattern in a sequence."""

        regexpr = ''.join([
            NUCLEOTIDE_TO_REGEXPR[n]
            for n in sequence
        ])
        return regexpr

    def all_variants(self):
        """Return all ATGC sequence variants of a sequence"""
        return[
            "".join(nucleotides)
            for nucleotides in itertools.product(*[
                IUPAC_NOTATION[n] for n in self.sequence
            ])
        ]

    def __repr__(self):
        """Represent the pattern as PatternType(name) """
        return self.sequence + ("" if self.name is None else
                                " (%s)" % self.name)
    def __str__(self):
        """Represent the pattern as PatternType(name) """
        return self.sequence + ("" if self.name is None else
                                " (%s)" % self.name)
    @staticmethod
    def from_string(string):
        if set(string) <= set(NUCLEOTIDE_TO_REGEXPR.keys()):
            return DnaNotationPattern(string)

# DEFINITION OF COMMON PATTERNS

class EnzymeSitePattern(DnaNotationPattern):
    """Class to represent Enzyme site patterns

    Examples
    --------

    >>> enzyme_pattern = EnzymeSitePattern("BsaI")
    >>> constraint = AvoidPattern(enzyme_pattern)

    """

    def __init__(self, enzyme_name):
        self.enzyme_site = rest_dict[enzyme_name]["site"]
        DnaNotationPattern.__init__(self, self.enzyme_site, name=enzyme_name)
    
    @staticmethod
    def from_string(string):
        """Convert BsmBI_site to EnzymeSitePattern(BsmBI)"""
        match = re.match("(\S+)_site", string)
        if match is not None:
            enzyme_name = match.groups()[0]
            if enzyme_name in rest_dict:
                return EnzymeSitePattern(enzyme_name)

    def __str__(self):
        return "%s(%s)" % (self.name, self.enzyme_site)

        


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
        DnaNotationPattern.__init__(self, number * nucleotide,
                                    in_both_strands=True)
    @staticmethod
    def from_string(string):
        match = re.match("(\d+)x(\S)$", string)
        if match is not None:
            number, nucleotide = match.groups()
            return HomopolymerPattern(nucleotide, int(number))

    def __str__(self):
        return "%sx%s" % (self.number, self.nucleotide)


class RepeatedKmerPattern(SequencePattern):
    """Direct repeats like ATT-ATT, ATGC-ATGC-ATGC, etc.

    Shorthand string version: "3x4mer", "5x2mer", etc.

    Examples
    --------

    >>> homopolymer("A", 6) # returns DnaNotationPattern("AAAAAA")

    """
    def __init__(self, n_repeats, kmer_size):
        self.n_repeats = n_repeats
        self.kmer_size = kmer_size
        SequencePattern.__init__(
            self, size=kmer_size * n_repeats,
            expression=r"([ATGC]{%d})\1{%d}" % (kmer_size, n_repeats-1),
            name="%d-repeats %d-mers" % (n_repeats, kmer_size),
            in_both_strands=False,  # a repeat on a strand is also on the other
            lookahead='loop')

    @staticmethod
    def from_string(string):
        match = re.match("(\d+)x(\d+)mer$", string)
        if match is not None:
            n_repeats, kmer_size = match.groups()
            return RepeatedKmerPattern(int(n_repeats), int(kmer_size))

    def __str__(self):
        return "%sx%smer" % (self.n_repeats, self.kmer_size)

SequencePattern.registered_string_pattern_classes = [
    HomopolymerPattern,
    RepeatedKmerPattern,
    EnzymeSitePattern,
    DnaNotationPattern
]