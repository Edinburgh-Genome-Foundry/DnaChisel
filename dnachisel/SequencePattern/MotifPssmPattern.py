from Bio.Seq import Seq
from Bio import motifs
from Bio.Align.AlignInfo import PSSM
from .SequencePattern import SequencePattern
import numpy as np


class MotifPssmPattern(SequencePattern):
    """Motif pattern represented by a Position-Specific Scoring Matrix (PSSM).

    This class is better instantiated by creating a pattern from a set of
    sequence strings with ``MotifPssmPattern.from_sequences()``,  or reading
    pattern(s) for a file with ``MotifPssmPattern.from_file()``.

    Parameters
    ----------

    pssm
      A Bio.Align.AlignInfo.PSSM object.

    threshold
     locations of the sequence with a PSSM score above this value will be
     considered matches. For convenience, a relative_threshold can be
     given instead.

    relative_threshold
      Value between 0 and 1 from which the threshold will be auto-computed.
      0 means "match everything", 1 means "only match the one (or several)
      sequence(s) with the absolute highest possible score".
    """

    def __init__(
        self, pssm, threshold=None, relative_threshold=None,
    ):
        self.name = pssm.name
        self.pssm = pssm.pssm
        if relative_threshold is not None:
            mini, maxi = self.pssm.min, self.pssm.max
            threshold = mini + relative_threshold * (maxi - mini)
        self.threshold = threshold
        self.relative_threshold = relative_threshold
        self.size = pssm.length
        self.pssm_matrix = np.array([self.pssm[n] for n in "ATGC"])
        self.is_palyndromic = False

    @classmethod
    def apply_pseudocounts(cls, motif, pseudocounts):
        """Add pseudocounts to the motif's pssm matrix.

        Add nothing if pseudocounts is None, auto-compute pseudocounts if
        pseudocounts="jaspar", else just attribute the pseudocounts as-is
        to the motif.
        """
        if pseudocounts is not None:
            if pseudocounts == "jaspar":
                pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)
            motif.pseudocounts = pseudocounts

    def find_matches_in_string(self, sequence):
        """Find all positions where the PSSM score is above threshold."""

        # NOTE: Before, I made my PSSM searches with Biopython. It was looong!
        # Now I use Numpy and np.choice(), and I never looked back
        # sequence = Seq(sequence, alphabet=alphabet)
        # search = self.pssm.search(
        #     sequence, threshold=self.threshold, both=False
        # )
        indices = find_pssm_matches_with_numpy(
            pssm_matrix=self.pssm_matrix, sequence=sequence, threshold=self.threshold,
        )
        return [(i, i + self.size, 1) for i in indices]

    @classmethod
    def from_sequences(
        cls,
        sequences,
        name="unnamed",
        pseudocounts="jaspar",
        threshold=None,
        relative_threshold=None,
    ):
        """Return a PSSM pattern computed from same-length sequences.

        Parameters
        ----------

        sequences
          A list of same-length sequences.

        name
          Name to give to the pattern (will appear in reports etc.).

        pseudocounts
          Either a dict {"A": 0.01, "T": ...} or "jaspar" for automatic
          pseudocounts from the Biopython.motifs.jaspar module (recommended),
          or None for no pseudocounts at all (not recommended!).

        threshold
          locations of the sequence with a PSSM score above this value will be
          considered matches. For convenience, a relative_threshold can be
          given instead.

        relative_threshold
          Value between 0 and 1 from which the threshold will be auto-computed.
          0 means "match everything", 1 means "only match the one (or several)
          sequence(s) with the absolute highest possible score".
        """
        sequences = [Seq(s) for s in sequences]
        motif = motifs.create(sequences)
        cls.apply_pseudocounts(motif, pseudocounts)
        pssm = PSSM(motif.pssm)
        pssm.name = name
        return MotifPssmPattern(
            pssm=pssm, threshold=threshold, relative_threshold=relative_threshold,
        )

    @classmethod
    def list_from_file(
        cls,
        motifs_file,
        file_format,
        threshold=None,
        pseudocounts="jaspar",
        relative_threshold=None,
    ):
        """Return a list of PSSM patterns from a file in JASPAR, MEME, etc.

        Parameters
        ----------

        motifs_file
          Path to a motifs file, or file handle.

        file_format
          File format. one of "jaspar", "meme", "TRANSFAC".

        pseudocounts
          Either a dict {"A": 0.01, "T": ...} or "jaspar" for automatic
          pseudocounts from the Biopython.motifs.jaspar module (recommended),
          or None for no pseudocounts at all (not recommended!).

        threshold
          locations of the sequence with a PSSM score above this value will be
          considered matches. For convenience, a relative_threshold can be
          given instead.

        relative_threshold
          Value between 0 and 1 from which the threshold will be auto-computed.
          0 means "match everything", 1 means "only match the one (or several)
          sequence(s) with the absolute highest possible score".
        """
        if isinstance(motifs_file, str):
            with open("./jaspar.txt", "r") as f:
                motifs_list = motifs.parse(f, file_format)
        else:
            motifs_list = motifs.parse(motifs_file, file_format)
        if pseudocounts is not None:
            for motif in motifs_list:
                cls.apply_pseudocounts(motif, pseudocounts)

        return [
            MotifPssmPattern(
                pssm, threshold=threshold, relative_threshold=relative_threshold,
            )
            for pssm in motifs_list
        ]

    def __str__(self):
        if self.relative_threshold is not None:
            threshold = "%d%%" % (100 * self.relative_threshold)
        else:
            threshold = "%.2f" % self.threshold
        return "%s-PSSM(%s+)" % (self.name, threshold)

    def __repr__(self):
        if self.relative_threshold is not None:
            threshold = "%d%%" % (100 * self.relative_threshold)
        else:
            threshold = "%.2f" % self.threshold
        return "%s-PSSM(%s+)" % (self.name, threshold)


def find_pssm_matches_with_numpy(pssm_matrix, sequence, threshold):
    """Return every index in the +1 strand wit a PSSM score above threshold.

    My numpy-based implementation is 10 times faster than Biopython for some
    reason. Weird. Can someone else check?

    Parameters:
    -----------

    pssm
      A matrix whose rows give the frequency motif of ATGC (in this order).

    sequence
      A string representing a DNA sequence.

    threshold
      Every index with a score above this threshold will be returned.
    """
    nucleotide_to_index = dict(zip("ATGC", range(4)))
    len_pattern = len(pssm_matrix[0])

    # If sequence is small, use normal python to avoid numpy overhead

    if len(sequence) < 60:
        nucl_indices = [nucleotide_to_index[n] for n in sequence]
        return [
            i
            for i in range(len(sequence) - len_pattern)
            if np.choose(nucl_indices[i : len_pattern + i], pssm_matrix).sum()
            >= threshold
        ]

    # If sequence is large, use Numpy for speed. tested experimentally

    nucl_indices = np.array([nucleotide_to_index[n] for n in sequence], dtype="uint8")
    len_pattern = len(pssm_matrix[0])
    scores = np.array(
        [
            np.choose(nucl_indices[k : len_pattern + k], pssm_matrix).sum()
            for k in range(len(sequence) - len_pattern)
        ]
    )
    return np.nonzero(scores >= threshold)[0]
