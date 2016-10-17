"""Collection of useful pre-defined objectives/constraints for DnaChisel."""

import copy
from collections import Counter, defaultdict

import numpy as np

from ..biotools.biotables import CODON_USAGE
from ..biotools.biotools import (gc_content, reverse_complement,
                                 sequences_differences, windows_overlap,
                                 blast_sequence, translate,
                                 subdivide_window)
from .Objective import (Objective, PatternObjective, TerminalObjective,
                        ObjectiveEvaluation, VoidObjective)


class CodonOptimize(Objective):
    """Objective to codon-optimize a coding sequence for a particular organism.

    Several codon-optimization policies exist. At the moment this Objective
    implements a method in which codons are replaced by the most frequent
    codon in the species.

    (as long as this doesn't break any Objective or lowers the global
    optimization objective)

    Supported organisms are ``E. coli``, ``S. cerevisiae``, ``H. Sapiens``,
    ``C. elegans``, ``D. melanogaster``, ``B. subtilis``.

    Parameters
    ----------

    organism
      Name of the organism to codon-optimize for. Supported organisms are
      ``E. coli``, ``S. cerevisiae``, ``H. Sapiens``, ``C. elegans``,
      ``D. melanogaster``, ``B. subtilis``.
      Note that the organism can be omited if a ``codon_usage_table`` is
      provided instead

    window
      Pair (start, end) indicating the position of the gene to codon-optimize.
      If not provided, the whole sequence is considered as the gene. Make
      sure the length of the sequence in the window is a multiple of 3.

    strand
      Either 1 if the gene is encoded on the (+) strand, or -1 for antisense.

    codon_usage_table
      A dict of the form ``{"TAC": 0.112, "CCT": 0.68}`` giving the RSCU table
      (relative usage of each codon). Only provide if no ``organism`` name
      was provided.

    Examples
    --------

    >>> objective = CodonOptimizationObjective(
    >>>     organism = "E. coli",
    >>>     window = (150, 300), # coordinates of a gene
    >>>     strand = -1
    >>> )


    """

    def __init__(self, organism=None, window=None, strand=1,
                 codon_usage_table=None, boost=1.0):
        Objective.__init__(self, boost=boost)
        self.window = window
        self.strand = strand
        self.organism = organism
        if organism is not None:
            codon_usage_table = CODON_USAGE[self.organism]
        if codon_usage_table is None:
            raise ValueError("Provide either an organism name or a codon "
                             "usage table")
        self.codon_usage_table = codon_usage_table

    def evaluate(self, canvas):

        window = (self.window if self.window is not None
                  else [0, len(canvas.sequence)])
        start, end = window
        subsequence = canvas.sequence[start:end]
        if self.strand == -1:
            subsequence = reverse_complement(subsequence)
        length = len(subsequence)
        if (length % 3):
            raise ValueError("CodonOptimizationObjective on a window/sequence"
                             "with size %d not multiple of 3)" % length)
        score = sum([
            self.codon_usage_table[subsequence[3 * i:3 * (i + 1)]]
            for i in range(length / 3)
        ])
        return ObjectiveEvaluation(
            self, canvas, score, windows=[[start, end]],
            message="Codon opt. on window %s scored %.02E" %
                    (str(window), score)
        )

    def __str__(self):
        return "CodonOptimize(%s, %s)" % (str(self.window), self.organism)


class MinimizeDifferences(Objective):
    """Objective to minimize the differences to a given sequence.


    This can be used to enforce "conservative" optimization, in which we try
    to minimize the changes from the original sequence

    Parameters
    ----------

    window
      Pair (start, end) indicating the segment of the sequence. If none
      provided, the whole sequence is considered.

    target_sequence
      The DNA sequence that the canvas' sequence (or subsequence) should equal.
      Can be omitted if ``original_sequence`` is provided instead

    original_sequence
      A DNA sequence (will generally be the canvas' sequence itself) with
      already the right sequence at the given ``window``. Only provide if
      you are not providing a ``target_sequence``


    Examples
    --------

    >>> from dnachisel import *
    >>> sequence = random_dna_sequence(length=10000)
    >>> # Fix the sequence's local gc content while minimizing changes.
    >>> canvas = DnaCanvas(
    >>>     sequence = sequence,
    >>>     Objectives = [GCContentObjective(0.3,0.6, gc_window=50)],
    >>>     objective = [MinimizeDifferencesObjective(
    >>>                     original_sequence=sequence)]
    >>> )
    >>> canvas.solve_all_Objectives_one_by_one()
    >>> canvas.maximize_all_objectives_one_by_one()

    """

    best_possible_score = 0

    def __init__(self, window=None, target_sequence=None,
                 original_sequence=None, boost=1.0):
        self.boost = boost
        self.window = window
        if target_sequence is None:
            target_sequence = (original_sequence if window is None else
                               original_sequence[window[0]:window[1]])
        self.reference_sequence = target_sequence

    def evaluate(self, canvas):
        window = (self.window if self.window is not None
                  else [0, len(canvas.sequence)])
        start, end = window
        subsequence = canvas.sequence[start: end]
        diffs = - sequences_differences(subsequence, self.reference_sequence)
        return ObjectiveEvaluation(
            self, canvas, score=-diffs, windows=[window],
            message="Found %d differences with target sequence" % diffs
        )

    def __str__(self):
        return "MinimizeDifferencesObj(%s, %s...)" % (
            "global" if self.window is None else str(self.window),
            self.sequence[:7]
        )


class DoNotModify(Objective):
    """Specify that some locations of the sequence should not be changed.

    ``DoNotModify`` Objectives are used to constrain the mutations space
    of DNA Canvas.

    Parameters
    ----------

    window
      Couple ``(start, end)`` indicating the position of the segment that
      must be left unchanged.
    """

    best_possible_score = 1

    def __init__(self, window=None, indices=None, boost=1.0):
        self.window = window
        self.indices = np.array(indices)
        self.boost = boost

    def evaluate(self, canvas):
        sequence = canvas.sequence
        original = canvas.original_sequence
        if (self.window is None) and (self.indices is None):
            return ObjectiveEvaluation(sequence == original)
        elif self.window is not None:
            start, end = self.window
            score = 1 if (sequence[start:end] == original[start:end]) else -1
            return ObjectiveEvaluation(self, canvas, score)
        else:
            sequence = np.fromstring(sequence, dtype="uint8")
            original = np.fromstring(original, dtype="uint8")
            if (sequence[self.indices] == original[self.indices]).min():
                score = 1
            else:
                score = -1

            return ObjectiveEvaluation(self, canvas, score)

    def localize(self, window):
        """Localize the DoNotModify to the overlap of its window and the new.
        """
        if self.window is not None:
            new_window = windows_overlap(self.window, window)
            if new_window is None:
                return VoidObjective(parent_objective=self)
            return self.copy_with_changes(window=new_window)
        else:
            start, end = window
            inds = self.indices
            new_indices = inds[(start <= inds) & (inds <= end)]
            return self.copy_with_changes(indices=new_indices)

    def __repr__(self):
        return "DoNotModify(%s)" % str(self.window)


class EnforcePattern(PatternObjective):
    """Enforce that the given pattern is present in the sequence.

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    window
      A couple (start, end) indicating the segment of DNA to which to restrict
      the search
    """

    def __init__(self, pattern, window=None, occurences=1, boost=1.0):
        PatternObjective.__init__(self, pattern, window)
        self.occurences = occurences
        self.boost = boost

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        windows = self.pattern.find_matches(canvas.sequence, window)
        score = -abs(len(windows) - self.occurences)

        if score == 0:
            message = "Passed. Pattern found at positions %s" % windows
        else:
            if self.occurences == 0:
                message = "Failed. Pattern not found."
            else:
                message = ("Failed. Pattern found %d times instead of %d"
                           " wanted at positions %s") % (len(windows),
                                                         self.occurences,
                                                         window)
        return ObjectiveEvaluation(
            self, canvas, score, message=message,
            windows=None if window is None else [window],
        )

    def __repr__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.window)


class AvoidPattern(PatternObjective):
    """Enforce that the given pattern is absent in the sequence.
    """

    def evaluate(self, canvas):
        windows = self.pattern.find_matches(canvas.sequence, self.window)
        score = -len(windows)
        if score == 0:
            message = "Passed. Pattern not found !"
        else:
            message = "Failed. Pattern found at positions %s" % windows
        return ObjectiveEvaluation(
            self, canvas, score, windows=windows, message=message
        )

    def __repr__(self):
        return "NoPattern(%s, %s)" % (self.pattern, self.window)


class EnforceTranslation(Objective):
    """Enforce that the DNA segment sequence translates to a specific
    amino-acid sequence.


    Parameters
    -----------

    window
      A pair (start, end) indicating the segment that is a coding sequence

    strand
      Set to 1 (default) if the gene is read in direct sense, -1 for antisense

    translation
      String representing the protein sequence that the DNA segment should
      translate to, eg. "MKY...LL*" ("*" stands for stop codon).
      This parameter can be omited if the ``sequence`` parameter is provided

    sequence
      A sequence of DNA that already encodes the right protein in the given
      ``window`` (will generally be equal to the sequence provided to
      the canvas if it already encodes the right protein).
      Can be provided instead of ``translation`` (the ``translation`` will be
      computed from this ``sequence``)

    Examples
    --------

    >>> from dnachisel import *
    >>> sequence = some_dna_sequence # with a gene in segment 150-300
    >>> Objective = EnforceTranslationObjective(
    >>>     window=(150,300),
    >>>     strand = 1,
    >>>     translation= translate(sequence[150:300]) # "MKKLYQ...YNL*"
    >>> )
    >>> # OR EQUIVALENT IF THE GENE ALREADY ENCODES THE RIGHT PROTEIN:
    >>> Objective = EnforceTranslationObjective(
    >>>     window=(150,300),
    >>>     strand = 1,
    >>>     sequence = sequence
    >>> )
    """

    best_possible_score = 1

    def __init__(self, window, sequence=None, translation=None, strand=1,
                 boost=1.0):

        self.boost = boost
        self.window = window
        if translation is None:
            start, end = window
            subsequence = sequence[start:end]
            if strand == -1:
                subsequence = reverse_complement(subsequence)
            translation = translate(subsequence)
        self.translation = translation
        self.strand = strand
        window_size = window[1] - window[0]
        if window_size != 3 * len(translation):
            raise ValueError(
                ("Window size (%d bp) incompatible with translation (%d aa)") %
                (window_size, len(translation))
            )

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        start, end = window
        subsequence = canvas.sequence[start:end]
        if self.strand == -1:
            subsequence = reverse_complement(subsequence)
        success = 1 if (translate(subsequence) == self.translation) else -1
        return ObjectiveEvaluation(self, canvas, success)

    def localized(self, window):
        """TODO: implement the localization of this one. At the moment
        it just voids the Objective in regions that do not overlap at all
        with its window."""
        if self.window is not None:
            overlap = windows_overlap(window, self.window)
            if overlap is None:
                return VoidObjective(parent_objective=self)
            else:
                return self
        return self

    def __repr__(self):
        return "EnforceTranslation(%s)" % str(self.window)


class EnforceGCContent(Objective):
    """Objective on the local or global proportion of G/C nucleotides.

    Examples
    --------

    >>> # Enforce global GC content between 40 and 70 percent.
    >>> Objective = GCContentObjective(0.4, 0.7)
    >>> # Enforce 30-80 percent local GC content over 50-nucleotides windows
    >>> Objective = GCContentObjective(0.3, 0.8, gc_window=50)


    Parameters
    ----------

    gc_min
      Minimal proportion of G-C (e.g. ``0.35``)

    gc_max
      Maximal proportion of G-C (e.g. ``0.75``)

    gc_window
      Length of the sliding window, in nucleotides, for local GC content.
      If not provided, the global GC content of the whole sequence is
      considered

    window
      Couple (start, end) indicating that the Objective only applies to a
      subsegment of the sequence. Make sure it is bigger than ``gc_window``
      if both parameters are provided

    """

    def __init__(self, gc_min=0, gc_max=1.0, gc_objective=None,
                 gc_window=None, window=None, boost=1.0):
        if gc_objective is not None:
            gc_min = gc_max = gc_objective
        self.gc_objective = gc_objective
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.gc_window = gc_window
        self.window = window
        self.boost = boost

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        wstart, wend = window
        sequence = canvas.sequence[wstart:wend]
        gc = gc_content(sequence, self.gc_window)
        breaches = (np.maximum(0, self.gc_min - gc) +
                    np.maximum(0, gc - self.gc_max))
        score = - (breaches.sum())
        breaches_starts = (breaches > 0).nonzero()[0]

        if len(breaches_starts) == 0:
            breaches_windows = []
        elif len(breaches_starts) == 1:
            if self.gc_window is not None:
                start = breaches_starts[0]
                breaches_windows = [[start, start + self.gc_window]]
            else:
                breaches_windows = [[wstart, wend]]
        else:
            breaches_windows = []
            current_start = breaches_starts[0]
            last_end = current_start + self.gc_window
            for i in breaches_starts[1:]:
                if (i > last_end + self.gc_window):
                    breaches_windows.append([
                        wstart + current_start, wstart + last_end]
                    )
                    current_start = i
                    last_end = i + self.gc_window

                else:
                    last_end = i + self.gc_window
            breaches_windows.append(
                [wstart + current_start, wstart + last_end])

        if breaches_windows == []:
            message = "Passed !"
        else:
            message = ("Failed: GC content out of bound on segments " +
                       ", ".join(["%s-%s" % (s[0], s[1])
                                  for s in breaches_windows]))
        return ObjectiveEvaluation(self, canvas, score, breaches_windows,
                                   message=message)

    def localized(self, window):
        """Localize the GC content evaluation

        For a window [start, end], the GC content evaluation will be restricted
        to [start - gc_window, end + gc_window]
        """
        if self.window is not None:
            new_window = windows_overlap(self.window, window)
            if new_window is None:
                return VoidObjective(parent_objective=self)
        else:
            start, end = window
            if self.gc_window is not None:
                new_window = [max(0, start - self.gc_window),
                              end + self.gc_window]
            else:
                new_window = None
        return self.copy_with_changes(window=new_window)

    def __repr__(self):
        return "GCContent(min %.02f, max %.02f, gc_win %s, window %s)" % (
            self.gc_min, self.gc_max, "global" if (self.gc_window is None) else
                                      self.gc_window, self.window
        )


class AvoidIDTHairpins(Objective):

    best_possible_score = 0

    def __init__(self, stem_size=20, hairpin_window=200, boost=1.0):

        self.stem_size = stem_size
        self.hairpin_window = hairpin_window
        self.boost = boost

    def evaluate(self, canvas):
        sequence = canvas.sequence
        reverse = reverse_complement(sequence)
        windows = []
        for i in range(len(sequence) - self.hairpin_window):
            word = sequence[i:i + self.stem_size]
            rest = reverse[-(i + self.hairpin_window):-(i + self.stem_size)]
            if word in rest:
                windows.append([i, i + self.hairpin_window])
        score = -len(windows)

        return ObjectiveEvaluation(self, canvas, score, windows=windows)

    def localized(self):
        # TODO: I'm pretty sure this can be localized
        return self

    def __repr__(self):
        return "NoHairpinsIDTObjective(size=%d, window=%d)" % \
            (self.stem_size, self.hairpin_window)


class TerminalGCContent(TerminalObjective):

    def __init__(self, gc_min, gc_max, window_size, boost=1.0):
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.window_size = window_size
        self.boost = boost

    def evaluate_end(self, sequence):
        return (self.gc_min < gc_content(sequence) < self.gc_max)

    def __repr__(self):
        return "Terminal(%.02f < gc < %.02f, window: %d)" % \
            (self.gc_min, self.gc_max, self.window_size)


class SequenceLengthBounds(Objective):
    """Checks that the sequence length is between bounds.

    Quite an uncommon objective as it can't really be solved or optimized.
    But practical at times.
    """

    def __init__(self, min_length=None, max_length=None):
        self.min_length = min_length
        self.max_length = max_length

    def evaluate(self, canvas):
        L, mini, maxi = len(canvas.sequence), self.min_length, self.max_length
        if mini is None:
            score = L <= maxi
        elif maxi is None:
            score = L >= mini
        else:
            score = (mini <= L <= maxi)
        return ObjectiveEvaluation(self, canvas, score - 1)

    def __repr__(self):
        return "Length(%d < L < %d)" % (self.min_length, self.max_length)


class AvoidBlastMatches(Objective):
    """Enforce that the given pattern is absent in the sequence.
    """

    def __init__(self, blast_db, word_size=4, perc_identity=100,
                 num_alignments=1000, num_threads=3, min_align_length=20,
                 window=None):
        self.blast_db = blast_db
        self.word_size = word_size
        self.perc_identity = perc_identity
        self.num_alignments = num_alignments
        self.num_threads = num_threads
        self.min_align_length = min_align_length
        self.window = window

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        wstart, wend = window
        sequence = canvas.sequence[wstart:wend]
        blast_record = blast_sequence(
            sequence, blast_db=self.blast_db,
            word_size=self.word_size,
            perc_identity=self.perc_identity,
            num_alignments=self.num_alignments,
            num_threads=self.num_threads
        )
        windows = sorted([
            sorted((hit.query_start + wstart, hit.query_end + wstart))
            for alignment in blast_record.alignments
            for hit in alignment.hsps
            if abs(hit.query_end - hit.query_start) >= self.min_align_length
        ])

        if windows == []:
            return ObjectiveEvaluation(self, canvas, score=1,
                                       message="Passed: no BLAST match found")

        return ObjectiveEvaluation(self, canvas, score=-len(windows),
                                   windows=windows,
                                   message="Failed - matches at %s" % windows)

    def localized(self, window):
        """Localize the evaluation """
        if self.window is not None:
            new_window = windows_overlap(self.window, window)
            if new_window is None:
                return VoidObjective(parent_objective=self)
        else:
            start, end = window
            radius = self.min_align_length
            new_window = [max(0, start - radius), end + radius]

        return self.copy_with_changes(window=new_window)

    def __repr__(self):
        return "NoBlastMatchesObjective%s(%s, %d+ bp, perc %d+)" % (
            self.window, self.blast_db, self.min_align_length,
            self.perc_identity
        )


class AvoidNonuniqueKmers(Objective):
    """
        from dnachisel import *
        sequence = random_dna_sequence(50000)
        canvas = DnaCanvas(
            sequence,
            constraints= [AvoidNonuniqueKmers(10,
                                    include_reverse_complement=True)]
        )
        print canvas.constraints_summary()
    """

    def __init__(self, length, window=None, include_reverse_complement=False):
        self.length = length
        self.window = window
        self.include_reverse_complement = include_reverse_complement

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        wstart, wend = window
        sequence = canvas.sequence[wstart:wend]
        rev_complement = reverse_complement(sequence)
        kmers_locations = defaultdict(lambda: [])
        for i in range(len(sequence) - self.length):
            start, end = i, i + self.length
            kmers_locations[sequence[start:end]].append((start, end))
        if self.include_reverse_complement:
            for i in range(len(sequence) - self.length):
                start, end = i, i + self.length
                kmers_locations[rev_complement[start:end]].append(
                    (len(sequence) - end, len(sequence) - start)
                )

        windows = sorted([
            min(positions_list, key=lambda p: p[0])
            for positions_list in kmers_locations.values()
            if len(positions_list) > 1
        ])

        if windows == []:
            return ObjectiveEvaluation(
                self, canvas, score=1,
                message="Passed: no nonunique %d-mer found." % self.length)

        return ObjectiveEvaluation(
            self, canvas, score=-len(windows),
            windows=windows,
            message="Failed, the following positions are the first occurences"
                    "of non-unique kmers %s" % windows)

    def __repr__(self):
        return "NoNonuniqueKmers(%d)" % (self.length)


#
# class GCContentObjective(Objective):
#     """Objective to obtain the desired global GC content.
#
#     Parameters
#     ----------
#
#     gc_content
#       GC content to target (e.g. ``0.4``)
#
#     window
#       Restrict the GC content observation to a DNA segment (start, end)
#
#     """
#
#     best_possible_score = 0
#
#     def __init__(self, gc_content, exponent=1.0, window=None,  boost=1.0,
#                  subdivision_window=None):
#         Objective.__init__(self, boost=boost)
#         self.gc_content = gc_content
#         self.exponent = exponent
#         self.window = window
#         self.subdivision_window = subdivision_window
#
#     def evaluate(self, canvas):
#         window = (self.window if self.window is not None
#                   else [0, len(canvas.sequence)])
#         start, end = window
#         subsequence = canvas.sequence[start: end]
#         gc = gc_content(subsequence)
#         score = -(abs(gc - self.gc_content) ** self.exponent)
#         if self.subdivision_window is None:
#             windows = [window]
#         else:
#             windows = subdivide_window(window, self.subdivision_window)
#
#         return ObjectiveEvaluation(
#             self, canvas, score=score, windows=windows,
#             message="scored %.02E. GC content is %.03f (%.03f wanted)" %
#                     (score, gc, self.gc_content))
#
#     def __str__(self):
#         return "GCContentObj(%.02f, %s)" % (
#             self.gc_content,
#             "global" if (self.window is None) else
#             ("window: %s" % str(self.window))
#         )
