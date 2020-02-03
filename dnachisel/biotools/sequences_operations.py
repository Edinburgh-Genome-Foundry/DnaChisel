import itertools

from Bio.Seq import Seq
from Bio.Data import CodonTable
import numpy as np

from .biotables import (
    NUCLEOTIDE_TO_REGEXPR,
    COMPLEMENTS,
    IUPAC_NOTATION,
    get_backtranslation_table,
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


def reverse_translate(
    protein_sequence, randomize_codons=False, table="Standard"
):
    """Return a DNA sequence which translates to the provided protein sequence.

    Parameters
    ----------

    protein_sequence
      A sequence string of aminoacids, e.g. "MVKK..."

    table
      Genetic code table to use (e.g. 'Standard', 'Bacterial', etc.).
      See dnachisel.biotools.CODON_TABLE_NAMES for a list of available genetic
      code tables.
    
    randomize_codons
      If False, the first valid codon found is used for each, which can create
      biases (GC content, etc.), if True, each amino acid gets replaced by a
      randomly selected codon for this amino acid.
    """
    backtranslation_table = get_backtranslation_table(table_name=table)
    if randomize_codons:
        random_numbers = np.random.randint(0, 1000, len(protein_sequence))
        random_indices = [
            random_number % len(backtranslation_table[aa])
            for aa, random_number in zip(protein_sequence, random_numbers)
        ]
        return "".join(
            [
                backtranslation_table[aa][random_indice]
                for aa, random_indice in zip(protein_sequence, random_indices)
            ]
        )
    return "".join([backtranslation_table[aa][0] for aa in protein_sequence])


def translate(dna_sequence, table="Standard", assume_start_codon=False):
    """Translate the DNA sequence into an amino-acids sequence "MLKYQT...".
    If ``translation_table`` is the name or number of a NCBI genetic table,
    Biopython will be used. See here for options:

    http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25

    ``translation_table`` can also be a dictionnary of the form
    ``{"ATT": "M", "CTC": "X", etc.}`` for more exotic translation tables

    If assume_start_codon is True and the first codon is a start codon in the
    given genetic table, then the first amino acid will be M (methionine).

    """
    if isinstance(table, dict):
        return "".join(
            [
                table[dna_sequence[i : i + 3]]
                for i in range(0, len(dna_sequence), 3)
            ]
        )
    else:
        table = CodonTable.unambiguous_dna_by_name[table]
        if assume_start_codon and (dna_sequence[:3] in table.start_codons):
            return "M" + str(Seq(dna_sequence[3:]).translate(table=table))
        return str(Seq(dna_sequence).translate(table=table))


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
