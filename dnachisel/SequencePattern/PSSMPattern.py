from Bio.Seq import Seq
from Bio import motifs
from Bio.Align.AlignInfo import PSSM
from Bio.Alphabet import IUPAC
from .SequencePattern import SequencePattern
import numpy as np

alphabet = IUPAC.IUPACUnambiguousDNA()


def find_pssm_matches_with_numpy(pssm_matrix, sequence, threshold):
    """Return every index in the +1 strand wit a PSSM score above threshold.
    
    My numpy-based implementation is 10 times faster than Biopython for some
    reason. Weird. Can someone else check?

    Parameters:
    -----------

    pssm
      The .pssm attribute of of a Biopython PSSM object (I know).
    
    sequence
      An "ATGC" string
    
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

    nucl_indices = np.array(
        [nucleotide_to_index[n] for n in sequence], dtype="uint8"
    )
    len_pattern = len(pssm_matrix[0])
    scores = np.array(
        [
            np.choose(nucl_indices[k : len_pattern + k], pssm_matrix).sum()
            for k in range(len(sequence) - len_pattern)
        ]
    )
    return np.nonzero(scores >= threshold)[0]


class PSSMPattern(SequencePattern):
    def __init__(
        self,
        pssm,
        threshold=None,
        relative_threshold=None,
    ):
        self.name = pssm.name
        self.pssm = pssm.pssm
        if relative_threshold is not None:
            mini, maxi = self.pssm.min, self.pssm.max
            threshold = mini + relative_threshold * (maxi - mini)
        self.threshold = threshold
        self.relative_threshold = relative_threshold
        self.length = pssm.length
        self.pssm_matrix = np.array([self.pssm[n] for n in "ATGC"])
        self.is_palyndromic = False

    @classmethod
    def apply_pseudocounts(cls, motif, pseudocounts):
        if pseudocounts is not None:
            if pseudocounts == "jaspar":
                pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)
            motif.pseudocounts = pseudocounts

    def find_matches_in_string(self, sequence):
        sequence = Seq(sequence, alphabet=alphabet)
        # NOTE: Before, I made my PSSM searches with Biopython. It was looong!
        # Now I use Numpy and np.choice(), and I never looked back
        # search = self.pssm.search(
        #     sequence, threshold=self.threshold, both=False
        # )
        indices = find_pssm_matches_with_numpy(
            pssm_matrix=self.pssm_matrix,
            sequence=sequence,
            threshold=self.threshold,
        )
        return [(i, i + self.length, 1) for i in indices]

    @classmethod
    def from_sequences(
        cls,
        sequences,
        name="unnamed",
        threshold=None,
        relative_threshold=None,
        pseudocounts="jaspar",
    ):
        sequences = [Seq(s) for s in sequences]
        motif = motifs.create(sequences)
        cls.apply_pseudocounts(motif, pseudocounts)
        pssm = PSSM(motif.pssm)
        pssm.name = name
        return PSSMPattern(
            pssm=pssm,
            threshold=threshold,
            relative_threshold=relative_threshold,
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
        if isinstance(motifs_file, str):
            with open("./jaspar.txt", "r") as f:
                motifs_list = motifs.parse(f, format=file_format)
        else:
            motifs_list = motifs.parse(motifs_file, format=file_format)
        if pseudocounts is not None:
            for motif in motifs_list:
                cls.apply_pseudocounts(motif, pseudocounts)

        return [
            PSSMPattern(
                pssm,
                threshold=threshold,
                relative_threshold=relative_threshold,
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
