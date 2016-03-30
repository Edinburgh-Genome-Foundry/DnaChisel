import copy
from biotools import gc_percent, translate, reverse_complement, windows_overlap
import numpy as np



class ConstraintEvaluation:

    def __init__(self, constraint, canvas, score, windows=None, message=None):
        self.constraint = constraint
        self.canvas = canvas
        self.score = score
        self.passes = score >= 0
        self.windows = windows
        self.message = message

    def __str__(self):
        return (
            self.message if (self.message is not None) else
            str((
                "Passes" if self.passes else "Fails",
                "score : %.02f" % self.score,
                self.windows
            ))
        )


class Constraint:

    def __init__(self, window=None):
        self.window = window

    def evaluate(self, sequence):
        pass

    def localized(self, window):
        return self

    def copy_with_changes(self, **kwargs):
        new_constraint = copy.deepcopy(self)
        new_constraint.__dict__.update(kwargs)
        return new_constraint


class PatternConstraint(Constraint):
    """Class for constraints such as presence or absence of a pattern."""

    def __init__(self, pattern, window=None):
        self.pattern = pattern
        self.window = window

    def localized(self, window):
        start, end = window
        pattern_size = self.pattern.size
        return self.copy_with_changes(window=[max(0,start - pattern_size),
                                              end + pattern_size])


class EnforcePatternConstraint(PatternConstraint):

    def __init__(self, pattern, window=None, occurences=1):
        PatternConstraint.__init__(self, pattern, window)
        self.occurences = occurences

    def evaluate(self, canvas):
        window = self.window
        if window is None:
            window = (0, len(canvas.sequence))
        windows = self.pattern.find_matches(canvas.sequence, window)
        score = -abs(len(windows) - self.occurences)
        return ConstraintEvaluation(self, canvas, score, windows=[window])

    def evaluation_message(self, evaluation):
        if evaluation.passes:
            return "Passed. Pattern found at positions %s" % (
                evaluation.windows
            )
        else:
            if self.occurences == 0:
                return "Failed. Pattern not found."
            else:
                return ("Failed. Pattern found %d times instead of %d wanted,"
                        " at positions %s") % (len(evaluation.windows),
                                               self.occurences,
                                               pattern_window)

    def __str__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.window)


class NoPatternConstraint(PatternConstraint):

    def evaluate(self, canvas):
        windows = self.pattern.find_matches(canvas.sequence, self.window)
        score = -len(windows)
        if score == 0:
            message= "Passed. Pattern not found !"
        else:
            message= "Failed. Pattern found at positions %s" % windows
        return ConstraintEvaluation(
            self, canvas, score, windows=windows, message=message
        )

    def __str__(self):
        return "NoPattern(%s, %s)" % (self.pattern, self.window)


class EnforceTranslationConstraint(Constraint):

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

    def __str__(self):
        return "EnforceTranslation(%s)" % str(self.window)


class GCContentConstraint(Constraint):

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
        gc = gc_percent(sequence, self.gc_window)
        breaches = np.maximum(0, self.gc_min - gc) + \
            np.maximum(0, gc - self.gc_max)
        score = - (breaches.sum())
        breaches_starts = (breaches > 0).nonzero()[0]


        if len(breaches_starts) == 0:
            breaches_windows = []
        elif len(breaches_starts) == 1:
            start = breaches_starts[0]
            breaches_windows = [[start, start + self.gc_window]]
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
            message= ("Failed: GC content out of bound on segments " +
                      ", ".join(["%s-%s" % (s[0], s[1])
                                for s in breaches_windows]))
        return ConstraintEvaluation(self, canvas, score, breaches_windows,
                                    message = message)

    def localized(self, window):
        if self.window is not None:
            new_window = windows_overlap(self.window, window)
        else:
            start, end = window
            if self.gc_window is not None:
                new_window = [max(0, start - self.gc_window),
                              end + self.gc_window]
            else:
                new_window = None
        return self.copy_with_changes(window=new_window)

    def __str__(self):
        return "GCContent(min %.02f, max %.02f, gc_win %s, window %s)" % (
            self.gc_min, self.gc_max, "global" if (self.gc_window is None) else
                                      self.gc_window, self.window
        )

class DoNotModifyConstraint(Constraint):

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

    def __str__(self):
        return "DoNotModify(%s)" % str(self.window)
