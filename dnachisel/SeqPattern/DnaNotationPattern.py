
import itertools
from ..biotools import (
    reverse_complement,
    NUCLEOTIDE_TO_REGEXPR,
    IUPAC_NOTATION,
)
from .SequencePattern import SequencePattern

class DnaNotationPattern(SequencePattern):
    """Class for patterns in plain DNA notation: ATTGCCA, GCNNKTA, etc.

    If the sequence is not palyndromic, the pattern will be looked for in
    both strands of sequences.
    """

    def __init__(self, sequence, name=None):
        """Initialize"""
        SequencePattern.__init__(
            self,
            size=len(sequence),
            expression=self.dna_sequence_to_regexpr(sequence),
            name=name,
            is_palyndromic=reverse_complement(sequence) == sequence,
        )
        self.sequence = sequence

    @staticmethod
    def dna_sequence_to_regexpr(sequence):
        """Return a regular expression to find the pattern in a sequence."""

        regexpr = "".join([NUCLEOTIDE_TO_REGEXPR[n] for n in sequence])
        return regexpr

    def all_variants(self):
        """Return all ATGC sequence variants of a sequence"""
        return [
            "".join(nucleotides)
            for nucleotides in itertools.product(
                *[IUPAC_NOTATION[n] for n in self.sequence]
            )
        ]

    def __repr__(self):
        """Represent the pattern as PatternType(name) """
        return self.sequence + (
            "" if self.name is None else " (%s)" % self.name
        )

    def __str__(self):
        """Represent the pattern as PatternType(name) """
        return self.sequence + (
            "" if self.name is None else " (%s)" % self.name
        )

    @staticmethod
    def from_string(string):
        if set(string) <= set(NUCLEOTIDE_TO_REGEXPR.keys()):
            return DnaNotationPattern(string)