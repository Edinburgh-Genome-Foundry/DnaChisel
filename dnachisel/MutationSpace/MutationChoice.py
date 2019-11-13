from ..biotools import windows_overlap
import itertools
import numpy as np


class MutationChoice:
    """Represent a segment of a sequence with several possible variants.

    Parameters
    ----------

    segment
      A pair (start, end) indicating the range of nucleotides concerned. We
      are applying Python range, so

    variants
      A set of sequence variants, at the given position

    Examples
    --------

    >>> choice = MutationChoice((70, 73), {})


    """

    __slots__ = ["segment", "start", "end", "variants", "is_any_nucleotide"]

    def __init__(self, segment, variants, is_any_nucleotide=False):
        if isinstance(segment, int):
            segment = (segment, segment + 1)
        self.segment = segment
        self.start, self.end = segment
        self.variants = variants
        self.is_any_nucleotide = is_any_nucleotide
        # self.possible_subsequences = set(m.subsequence for m in mutations)

    def random_variant(self, sequence):
        """Return one of the variants, randomly."""
        subsequence = sequence[self.start : self.end]
        variants = [v for v in self.variants if v != subsequence]
        # the sorting of variants seems essential to ensure reproducibility
        # between sessions.
        # it does not slow down the global algorithm (or less than 3%)
        variants = sorted(variants)
        return variants[np.random.randint(len(variants))]

    def merge_with(self, others):
        """Merge this mutation choice with others to form a single choice

        Examples:
        ---------

        >>> ((2, 5), {'ATT', 'ATA'})

       merged with:

        >>> [
        >>>     ((0, 3), {'GTA', 'GCT', 'GTT'}),
        >>>     ((3, 4), {'A', 'T', 'G', 'C'}),
        >>>     ((4, 7), {'ATG', 'ACC', 'CTG'})
        >>> ]

        returns the only choices on the full interval which are compatible with
        at least one choice in each of the MutationChoices
        >>> (0, 7), {'GTATACC', 'GTATATG'}

        """
        others = sorted(others, key=lambda o: o.start)
        others_start = others[0].start
        final_segment = others_start, others[-1].end
        final_variants = set()
        for candidate in self.variants:
            slots = []
            for other in others:
                istart, iend = windows_overlap(other.segment, self.segment)
                slot = []
                for variant in other.variants:
                    subseq = variant[istart - other.start : iend - other.start]
                    subcandidate = candidate[
                        istart - self.start : iend - self.start
                    ]
                    if subseq == subcandidate:
                        slot.append(variant)
                slots.append(slot)
            for subseqs in itertools.product(*slots):
                seq = "".join(subseqs)
                matching_seq = seq[
                    self.start - others_start : self.end - others_start
                ]
                if matching_seq == candidate:
                    final_variants.add(seq)
        return MutationChoice(segment=final_segment, variants=final_variants)

    def extract_varying_region(self):
        """Return MutationChoices for the central varying region and 2 flanks.

        For instance:

        >>> choice = MutationChoice((5, 12), [
        >>>     'ATGCGTG',
        >>>     'AAAAATG',
        >>>     'AAATGTG',
        >>>     'ATGAATG',
        >>> ])
        >>> choice.extract_varying_region()

        Result :

        >>> [
        >>>     MutChoice(5-6 A),
        >>>     MutChoice(6-10 TGCG-AATG-TGAA-AAAA),
        >>>     MutChoice(10-12 TG)
        >>> ]

        """

        if len(self.variants) <= 1:
            return [self]
        variants = list(self.variants)
        reference = variants[0]
        start = -1
        end = len(reference)
        for i in range(len(reference)):
            for variant in variants[1:]:
                if variant[i] != reference[i]:
                    if start == -1:
                        start = i
                    end = i + 1
                    break
        result = []
        if start > 0:
            result.append(
                MutationChoice(
                    (self.start, self.start + start), set([reference[:start]])
                )
            )
        result.append(
            MutationChoice(
                (self.start + start, self.start + end),
                set([v[start:end] for v in variants]),
            )
        )
        if end < len(reference):
            result.append(
                MutationChoice(
                    (self.start + end, self.end),
                    set([v[end:] for v in variants]),
                )
            )
        return result

    def __repr__(self):
        """Represent."""
        subsequences = "-".join(self.variants)
        return "MutChoice(%d-%d %s)" % (self.start, self.end, subsequences)

    def __str__(self):
        """Represent."""
        subsequences = "-".join(self.variants)
        return "MutChoice(%d-%d %s)" % (self.start, self.end, subsequences)
