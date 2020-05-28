"""Methods for finding indices where two sequences differ."""

import numpy as np


def sequences_differences_array(seq1, seq2):
    """Return an array [0, 0, 1, 0, ...] with 1s for sequence differences.

    seq1, seq2 should both be ATGC strings.
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            "Only use on same-size sequences (%d, %d)" % (len(seq1), len(seq2))
        )
    arr1 = np.frombuffer(seq1.encode(), dtype="uint8")
    arr2 = np.frombuffer(seq2.encode(), dtype="uint8")
    return arr1 != arr2


def sequences_differences(seq1, seq2):
    """Return the number of nucleotides that differ in the two sequences.

    seq1, seq2 should be strings of DNA sequences e.g. "ATGCTGTGC"
    """
    return sequences_differences_array(seq1, seq2).sum()


def sequences_differences_segments(seq1, seq2):
    """Return the list of segments on which sequence seq1 differs from seq2.

    The list is of the form [(start1, end1), (start2, end2), etc.]

    Parameters
    ----------

    seq1, seq2
      ATGC sequences to be compared
    """
    arr = 1 * sequences_differences_array(seq1, seq2)
    diffs = np.diff([0] + list(arr) + [0]).nonzero()[0]
    half = int(len(diffs) / 2)
    return [(diffs[2 * i], diffs[2 * i + 1]) for i in range(half)]
