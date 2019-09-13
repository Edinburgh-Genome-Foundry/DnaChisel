import itertools

from Bio.Seq import Seq
import numpy as np

from .biotables import (
    NUCLEOTIDE_TO_REGEXPR,
    COMPLEMENTS,
    CODONS_SEQUENCES,
    IUPAC_NOTATION,
)


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses BioPython for speed.
    """
    if len(dna_sequence) <= 30:
        return "".join([COMPLEMENTS[nuc] for nuc in dna_sequence])
    # This alternative has overhead but is really fast on long sequences
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.

    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def reverse_translate(protein_sequence, randomize_codons=False):
    """Return a DNA sequence which translates to the provided protein sequence

    Note: at the moment, the first valid codon found is used for each
    amino-acid (so it is deterministic but no codon-optimization is done).
    """
    if randomize_codons:
        random_indices = np.random.randint(0, 1000, len(protein_sequence))
        return "".join(
            [
                CODONS_SEQUENCES[aa][random_indice % len(CODONS_SEQUENCES[aa])]
                for aa, random_indice in zip(protein_sequence, random_indices)
            ]
        )
    return "".join([CODONS_SEQUENCES[aa][0] for aa in protein_sequence])


def translate(dna_sequence, translation_table="Bacterial"):
    """Translate the DNA sequence into an amino-acids sequence "MLKYQT...".
    If ``translation_table`` is the name or number of  NCBI genetic table,
    Biopython will be used. See here for options:

    http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25

    ``translation_table`` can also be a dictionnary of the form
    ``{"ATT": "M", "CTC": "X", etc.}`` for more exotic translation tables


    """
    if isinstance(translation_table, dict):
        return "".join(
            [
                translation_table[dna_sequence[i : i + 3]]
                for i in range(0, len(dna_sequence), 3)
            ]
        )
    else:
        return str(Seq(dna_sequence).translate(table=translation_table))


def dna_pattern_to_regexpr(dna_pattern):
    """Return a regular expression pattern for the provided DNA pattern

    For instance ``dna_pattern_to_regexpr('ATTNN')`` returns
    ``"ATT[A|T|G|C][A|T|G|C]"``.
    """
    return "".join([NUCLEOTIDE_TO_REGEXPR[n] for n in dna_pattern])


def all_iupac_variants(iupac_sequence):
    """Return all unambiguous possible versions of the given sequence."""
    return [
        "".join(nucleotides)
        for nucleotides in itertools.product(
            *[IUPAC_NOTATION[n] for n in iupac_sequence]
        )
    ]
