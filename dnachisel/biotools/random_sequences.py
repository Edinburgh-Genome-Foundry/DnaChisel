"""Methods to generate random sequences.

Also see utils.random_compatible_sequence() for random sequences verifying
user-provided constraints.
"""

import numpy as np


def random_dna_sequence(length, gc_share=None, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)
    if gc_share is not None:
        g_or_c = gc_share / 2.0
        not_g_or_c = (1 - gc_share) / 2.0
        probas = {"G": g_or_c, "C": g_or_c, "A": not_g_or_c, "T": not_g_or_c}
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def random_protein_sequence(length, seed=None):
    """Return a random protein sequence "MNQTW...YL*" of the specified length.

    Parameters
    ----------

    length
      Length of the protein sequence (in number of amino-acids). Note that the
      sequence will always start with ``"M"`` and end with a stop codon ``"*"``
      with (length-2) random amino-acids in the middle

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)

    aa_list = list("ACEDGFIHKLNQPSRTWVY")
    aa_choices = np.random.choice(aa_list, length - 2)
    return "M" + "".join(aa_choices) + "*"
