import copy
from biotools import (gc_content, translate, reverse_complement,
                      windows_overlap, blast_sequence)
import numpy as np


class ConstraintEvaluation:
    """Store relevant infos about the evaluation of a constraint on a canvas

    Examples
    --------

    >>> evaluation_result = ConstraintEvaluation(
    >>>     constraint = constraint,
    >>>     canvas = canvas,
    >>>     score= evaluation_score, # float
    >>>     windows=[(w1_start, w1_end), (w2_start, w2_end)...],
    >>>     message = "Failed: constraint not met on windows (122,155),..."
    >>> )

    Parameters
    ----------

    constraint
      The Constraint object that was evaluated.

    canvas
      The canvas that the constraint was evaluated on.

    score
      The score associated to the evaluation. Note that for constraints, a
      positive score (0 or more) means that the constraint is passing, and it
      is expected that when it is negative, the distance to 0 reflects how far
      the constraint is from passing (this is used by the constraints solver).

    windows
      A list of couples (start, end) indicating the windows on which the
      constraints breaches were observed.
      This parameter can be left to None if the constraint breaches are not
      localized (for global constraints such as global GC content).

    message
      A message that will be returned by ``str(evaluation)``. It will notably
      be displayed by ``canvas.print_constraints_summaries``.

    """

    def __init__(self, constraint, canvas, score, windows=None, message=None):
        self.constraint = constraint
        self.canvas = canvas
        self.score = score
        self.passes = score >= 0
        self.windows = windows
        self.message = message

    def __repr__(self):
        return (
            self.message if (self.message is not None) else
            str((
                "Passes" if self.passes else "Fails",
                "score : %.02f" % self.score,
                self.windows
            ))
        )


class Constraint:
    """General class to define constraints.

    New types of constraints are defined by subclassing ``Constraint`` and
    providing a custom ``evaluate`` and ``localized`` methods.

    """

    def __init__(self, window=None):
        self.window = window

    def evaluate(self, canvas):
        """Evaluate the constraint on the provided canvas.

        Returns a ``ConstraintEvaluation`` object indicating whether the
        constraint passed, the score of the evaluation, the windows on which
        to focus searches in case the constraint failed and a string message
        summarizing the evaluation (see ``ConstraintEvaluation``).
        """
        pass

    def localized(self, window):
        """Return a modified version of the constraint for the case where
        sequence modifications are only performed inside the provided window.

        For instance if a constraint concerns local GC content, and we are
        only making local mutations to destroy a restriction site, then we only
        need to check the local GC content around the restriction site after
        each mutation (and not compute it for the whole sequence), so
        ``GCContentConstraint.localised(window)`` will return a constraint
        that only looks for GC content around the provided window.

        If a constraint concerns a DNA segment that is completely disjoint from
        the provided window, this must return a ``VoidConstraint``.

        Must return an object of class ``Constraint``.
        """
        return self

    def copy_with_changes(self, **kwargs):
        new_constraint = copy.deepcopy(self)
        new_constraint.__dict__.update(kwargs)
        return new_constraint


class VoidConstraint(Constraint):
    """Void Constraints are a special case of constraints that always pass.

    Void Constraints are generally obtained when a constraint is "made void"
    by a localization. For instance, if we are optimizing the segment (10,50)
    of a DNA segment, the constraint EnforceTranslation([300,500]) does not
    apply as it concerns a Gene that is in a completely different segment.
    Therefore the localized version of EnforceTranslation will be void.
    """

    def __init__(self, parent_constraint=None):
        self.parent_constraint = parent_constraint
        self.message = ("Pass (not relevant in this context)"
                        % parent_constraint)

    def evaluate(self, canvas):
        """The evaluation of VoidConstraints always passes with score=1.0
        It returns a message indicating that the parent constraint was voided
        """
        return ConstraintEvaluation(self, canvas, score=1.0,
                                    message=self.message,
                                    windows=None)

    def __repr__(self):
        return "Void %s" % self.parent_constraint


class PatternConstraint(Constraint):
    """Class for constraints such as presence or absence of a pattern."""

    def __init__(self, pattern, window=None):
        self.pattern = pattern
        self.window = window

    def localized(self, window):
        pattern_size = self.pattern.size
        if self.window is not None:
            overlap = windows_overlap(self.window, window)
            if overlap is None:
                return VoidConstraint(parent_constraint=self)
            else:
                start, end = window
                ostart, oend = overlap
                if ostart == start:
                    new_window = [start, min(end, oend + pattern_size)]
                    #new_window = [start, oend + pattern_size]
                else:  # oend = end
                    new_window = [max(start, ostart - pattern_size), end]
                    #new_window = [ostart - pattern_size, end]
        else:
            start, end = window
            margin = pattern_size - 1
            new_window = [max(0, start - margin), end + margin]

        return self.copy_with_changes(window=new_window)


class EnforcePatternConstraint(PatternConstraint):
    """Enforce that the given pattern is present in the sequence.

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    window
      A couple (start, end) indicating the segment of DNA to which to restrict
      the search
    """

    def __init__(self, pattern, window=None, occurences=1):
        PatternConstraint.__init__(self, pattern, window)
        self.occurences = occurences

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
        return ConstraintEvaluation(
            self, canvas, score, message=message,
            windows=None if window is None else [window],
        )

    def __repr__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.window)


class NoPatternConstraint(PatternConstraint):
    """Enforce that the given pattern is absent in the sequence.
    """

    def evaluate(self, canvas):
        windows = self.pattern.find_matches(canvas.sequence, self.window)
        score = -len(windows)
        if score == 0:
            message = "Passed. Pattern not found !"
        else:
            message = "Failed. Pattern found at positions %s" % windows
        return ConstraintEvaluation(
            self, canvas, score, windows=windows, message=message
        )

    def __repr__(self):
        return "NoPattern(%s, %s)" % (self.pattern, self.window)


class EnforceTranslationConstraint(Constraint):
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
    >>> constraint = EnforceTranslationConstraint(
    >>>     window=(150,300),
    >>>     strand = 1,
    >>>     translation= translate(sequence[150:300]) # "MKKLYQ...YNL*"
    >>> )
    >>> # OR EQUIVALENT IF THE GENE ALREADY ENCODES THE RIGHT PROTEIN:
    >>> constraint = EnforceTranslationConstraint(
    >>>     window=(150,300),
    >>>     strand = 1,
    >>>     sequence = sequence
    >>> )



    """

    def __init__(self, window, sequence=None, translation=None, strand=1):
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
        return ConstraintEvaluation(self, canvas, success)

    def localized(self, window):
        """TODO: implement the localization of this one. At the moment
        it just voids the constraint in regions that do not overlap at all
        with its window."""
        if self.window is not None:
            overlap = windows_overlap(window, self.window)
            if overlap is None:
                return VoidConstraint(parent_constraint=self)
            else:
                return self
        return self

    def __repr__(self):
        return "EnforceTranslation(%s)" % str(self.window)


class GCContentConstraint(Constraint):
    """Constraint on the local or global proportion of G/C nucleotides.

    Examples
    --------

    >>> # Enforce global GC content between 40 and 70 percent.
    >>> constraint = GCContentConstraint(0.4, 0.7)
    >>> # Enforce 30-80 percent local GC content over 50-nucleotides windows
    >>> constraint = GCContentConstraint(0.3, 0.8, gc_window=50)


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
      Couple (start, end) indicating that the constraint only applies to a
      subsegment of the sequence. Make sure it is bigger than ``gc_window``
      if both parameters are provided

    """

    def __init__(self, gc_min, gc_max, gc_window=None, window=None):
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.gc_window = gc_window
        self.window = window

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        wstart, wend = window
        sequence = canvas.sequence[wstart:wend]
        gc = gc_content(sequence, self.gc_window)
        breaches = np.maximum(0, self.gc_min - gc) + \
            np.maximum(0, gc - self.gc_max)
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
        return ConstraintEvaluation(self, canvas, score, breaches_windows,
                                    message=message)

    def localized(self, window):
        """Localize the GC content evaluation

        For a window [start, end], the GC content evaluation will be restricted
        to [start - gc_window, end + gc_window]
        """
        if self.window is not None:
            new_window = windows_overlap(self.window, window)
            if new_window is None:
                return VoidConstraint(parent_constraint=self)
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


class DoNotModifyConstraint(Constraint):
    """Specify that a segment of the sequence should not be changed.

    ``DoNotModify`` Constraints are used to constrain the mutations space
    of DNA Canvas.

    Parameters
    ----------

    window
      Couple ``(start, end)`` indicating the position of the segment that
      must be left unchanged.
    """

    def __init__(self, window):
        self.window = window

    def evaluate(self, canvas):
        sequence = canvas.sequence
        original = canvas.original_sequence
        if self.window is None:
            return ConstraintEvaluation(sequence == original)
        else:
            start, end = self.window
            score = 1 if (sequence[start:end] == original[start:end]) else -1
            return ConstraintEvaluation(self, canvas, score)

    def localize(self, window):
        """Localize the DoNotModify to the overlap of its window and the new.
        """
        new_window = windows_overlap(self.window, window)
        if new_window is None:
            return VoidConstraint(parent_constraint=self)
        return self.copy_with_changes(window=new_window)

    def __repr__(self):
        return "DoNotModify(%s)" % str(self.window)


class NoHairpinsIDTConstraint(Constraint):

    def __init__(self, stem_size=20, hairpin_window=200):

        self.stem_size = stem_size
        self.hairpin_window = hairpin_window

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

        return ConstraintEvaluation(self, canvas, score, windows=windows)

    def __repr__(self):
        return "NoHairpinsIDTConstraint(size=%d, window=%d)" % \
            (self.stem_size, self.hairpin_window)


class TerminalConstraint(Constraint):

    def evaluate(self, canvas):
        sequence = canvas.sequence
        L = len(sequence)
        wsize = self.window_size
        ends_sequences = [
            ((0, wsize), sequence[:wsize]),
            ((L - wsize, L), sequence[-wsize:])
        ]
        windows = []
        for window, sequence in ends_sequences:
            if not self.evaluate_end(sequence):
                windows.append(window)

        if windows == []:
            message = "Passed (no breach at the ends)"
        else:
            message = "Failed: breaches at ends %s" % str(windows)

        return ConstraintEvaluation(self, canvas, score=len(windows),
                                    windows=windows, message=message)


class TerminalGCContentConstraint(TerminalConstraint):

    def __init__(self, gc_min, gc_max, window_size):
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.window_size = window_size

    def evaluate_end(self, sequence):
        return (self.gc_min < gc_content(sequence) < self.gc_max)

    def __repr__(self):
        return "Terminal(%.02f < gc < %.02f, window: %d)" % \
            (self.gc_min, self.gc_max, self.window_size)


class SequenceLengthConstraint(Constraint):

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
        return ConstraintEvaluation(self, canvas, score - 1)

    def __repr__(self):
        return "Length(%d < L < %d)" % (self.min_length, self.max_length)


class NoBlastMatchConstraint(Constraint):
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
            if abs(hit.query_end - hit.query_start) > self.min_align_length
        ])

        if windows == []:
            return ConstraintEvaluation(self, canvas, score=1,
                                        message="Passed: no BLAST match found")

        return ConstraintEvaluation(self, canvas, score=-len(windows),
                                    windows=windows,
                                    message="Failed - matches at %s" % windows)

    def localized(self, window):
        """Localize the evaluation

        For a window [start, end], the GC content evaluation will be restricted
        to [start - r, end + r]
        """
        if self.window is not None:
            new_window = windows_overlap(self.window, window)
            if new_window is None:
                return VoidConstraint(parent_constraint=self)
        else:
            start, end = window
            new_window = [max(0, start - self.min_align_length),
                          end + self.min_align_length]

        return self.copy_with_changes(window=new_window)

    def __repr__(self):
        return "NoBlastMatchesConstraint%s(%s, %d+ bp, perc %d+)" % (
            self.window, self.blast_db, self.min_align_length,
            self.perc_identity
        )
