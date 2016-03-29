
import biotables
import re
import numpy as np


def complement(sequence):
    return "".join([biotables.COMPLEMENTS[c] for c in sequence])

def reverse_complement(sequence):
    return complement(sequence)[::-1]

def reverse_translate(protein_sequence):
    return "".join([
        biotables.CODONS_SEQUENCES[aa][0]
        for aa in protein_sequence
    ])

def translate(dna_sequence):
    return "".join([
        biotables.CODON_TRANSLATIONS[dna_sequence[3*k:3*(k+1)]]
        for k in range(len(dna_sequence)/3)
    ])

def to_regexpr(self):
    """Return a regular expression pattern the sequence in ATGC sequences.

    For instance DNASequence('ATTNN').to_regexpr() returns
    "ATT[A|T|G|C][A|T|G|C]".
    """
    return "".join([
        biotables.NUCLEOTIDE_TO_REGEXPR[n]
        for n in self
    ])

def windows_overlap(window1, window2):
    """Return the overlap span between two windows.

    Parameters
    ----------

    window1, window2
      Each window is a couple of the form (start, end) indicating the range of
      a segment of integers.

    Returns
    -------

    None
      In case the two windows do not overlap.
    [start, end]
      The coordinates of the overlap segment if there is one.
    """
    start1, end1 = window1
    start2, end2 = window2

    if start2 < start1 :
        return windows_overlap(window2, window1)

    if start1 <= start2 <= end1:
        return [ start2, min(end1, end2)]
    else:
        return None

assert(windows_overlap([1,5], [3,7]) == [3,5])
assert(windows_overlap([3,7] ,[1,5]) == [3,5])
assert(windows_overlap([1,5] ,[1,5]) == [1,5])
assert(windows_overlap([1,5] ,[6,10]) is None)

def read_fasta(filename):
    """Read A sequence in a FASTA file with Biopython."""
    import Bio.SeqIO as seqio
    with open(filename) as f:
        return str(seqio.read(f, "fasta").seq)

def gc_percent(sequence, window_size = None):
    """Compute global or local GC content.

    Parameters
    ----------

    sequence
      An ATGC DNA sequence (upper case!)

    window_size
      If provided, the local GC content for the different sliding windows of
      this size is returned, else the global GC content is returned.

    Returns
    --------

      A number between 0 and 1 indication the proportion
      of GC content. If window_size is provided, returns
      a list of len(sequence)-window_size values indicating
      the local GC contents (sliding-window method). The i-th value
      indicates the GC content in the window [i, i+window_size]
    """
    # The code is a little cryptic but speed gain is 300x
    # compared with pure-python string operations

    arr = np.fromstring(sequence+"", dtype="uint8")
    arr_GCs = (arr == 71) | (arr == 67) # 67=C, 71=G

    if window_size is None:
        return 1.0 * arr_GCs.sum() / len(sequence)
    else:
        cs = np.cumsum(arr_GCs)
        a = cs[window_size-1:]
        b = np.hstack([[0],cs[:-window_size]])
        return 1.0*(a-b) / window_size
