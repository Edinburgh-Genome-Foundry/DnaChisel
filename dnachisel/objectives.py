from biotables import CODON_USAGE
from biotools import gc_percent, reverse_complement, sequences_differences
import numpy as np


class ObjectiveEvaluation:

    def __init__(self, objective, canvas, score, windows=None, message=None):
        self.objective = objective
        self.canvas = canvas
        self.score = score
        self.windows = windows
        self.message = message

    def __str__(self):
        return (
            self.message if (self.message is not None) else
            "score: %.02E, windows: %s" % (self.score, str(self.windows))
        )


class Objective:

    def __init__(self, boost=1. - 0):
        self.boost = boost

    def localized(self, window):
        return self


class CodonOptimizationObjective(Objective):

    def __init__(self, organism, window=None, strand=1, boost=1.0):
        Objective.__init__(self, boost=boost)
        self.window = window
        self.strand = strand
        self.organism = organism

    def evaluate(self, canvas):

        usage = CODON_USAGE[self.organism]
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
            usage[subsequence[3 * i:3 * (i + 1)]]
            for i in range(length / 3)
        ])
        return ObjectiveEvaluation(
            self, canvas, score, windows=[[start, end]],
            message="Codon opt. on window %s scored %.02E" %
                    (str(window), score)
        )

    def __str__(self):
        return "CodonOptimize(%s, %s)" % (str(self.window), self.organism)


class GCContentObjective(Objective):

    def __init__(self, gc_percent, exponent=1.0, window=None,  boost=1.0):
        Objective.__init__(self, boost=boost)
        self.gc_percent = gc_percent
        self.exponent = exponent
        self.window = window

    def evaluate(self, canvas):
        window = (self.window if self.window is not None
                  else [0, len(canvas.sequence)])
        start, end = window
        subsequence = canvas.sequence[start: end]
        gc = gc_percent(subsequence)
        score = -(abs(gc - self.gc_percent) ** self.exponent)
        return ObjectiveEvaluation(self, canvas, score=score, windows=[window],
                                   message="scored %.02E. GC content is %.03f (%.03f wanted)" %
                                   (score, gc, self.gc_percent))

    def __str__(self):
        return "GCContentObj(%.02f, %s)" % (
            self.gc_percent,
            "global" if (self.window is None) else
            ("window: %s" % str(self.window))
        )


class MinimizeDifferencesObjective(Objective):

    def __init__(self, window=None, sequence=None, original_sequence=None,
                 boost=1.0):
        Objective.__init__(self, boost=boost)
        self.window = window
        if sequence is None:
            sequence = original_sequence[window[0]:window[1]]
        self.sequence = sequence

    def evaluate(self, canvas):
        window = (self.window if self.window is not None
                  else [0, len(canvas.sequence)])
        start, end = window
        subsequence = canvas.sequence[start: end]
        diffs = - sequences_differences(subsequence, self.sequence)
        return ObjectiveEvaluation(
            self, canvas, score=-diffs, windows=[window],
            message="Found %d differences with target sequence" % diffs
        )

    def __str__(self):
        return "MinimizeDifferencesObj(%s, %s...)" % (
            "global" if self.window is None else str(self.window),
            self.sequence[:7]
        )
