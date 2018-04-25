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

    def __init__(self, expression, size=None, name=None, in_both_strands=True,
                 use_lookahead=True):
        if size is None:
            size = len(expression)
        self.expression = expression
        if use_lookahead:
            expression = '(?=(%s))' % expression
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
                loc + location.start
                for loc in self.find_matches(subsequence)
            ]
        matches = [
            (match.start(), match.start() + len(match.groups()[0]), 1)
            for match in re.finditer(self.compiled_expression, sequence)
        ]

        if self.in_both_strands:
            reverse = reverse_complement(sequence)
            L = len(sequence)
            matches += [
                (L - match.start() - len(match.groups()[0]),
                 L - match.start(), -1)
                for match in re.finditer(self.compiled_expression, reverse)
            ]

        return [
            Location(start, end, strand)
            for start, end, strand in matches
        ]

    def __str__(self):
        return self.expression + ("" if self.name is None else
                                  " (%s)" % self.name)


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
        # The ?= implements 'lookahead' which enables to find overlapping
        # patterns. It is followed by (.{S}) where S is the sequence size
        # to make sure the full group is selected.
        # return '(?=(%s))(.\{%d\})' % (regexpr, len(sequence))

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

# DEFINITION OF COMMON PATTERNS


def homopolymer_pattern(nucleotides, number):
    """Return a DnaNotationPattern with the sequence of a homopolymer.

    Examples
    --------

    >>> homopolymer("A", 6) # returns DnaNotationPattern("AAAAAA")

    """
    if len(nucleotides) == 1:
        return DnaNotationPattern(number * nucleotides)
    else:
        return SequencePattern(number * ("[%s]" % ("|".join(nucleotides))),
                               size=number)

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

    # FIXME: this regular expression does not support lookahead so we will
    # need the following special pattern.
    #     def myfindall(regex, seq):
    # ...    resultlist=[]
    # ...    pos=0
    # ...
    # ...    while True:
    # ...       result = regex.search(seq, pos)
    # ...       if result is None:
    # ...          break
    # ...       resultlist.append(seq[result.start():result.end()])
    # ...       pos = result.start()+1
    # ...    return resultlist

    return SequencePattern(
        size=kmer_size * n_repeats,
        expression=r"([ATGC]{%d})\1{%d}" % (kmer_size, n_repeats-1),
        name="%d-repeats %d-mers" % (n_repeats, kmer_size),
        in_both_strands=False,  # a kmer repeat one strand is also on the other
        use_lookahead=False
    )
