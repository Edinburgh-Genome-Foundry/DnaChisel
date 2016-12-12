import copy

import numpy as np

from ..biotools.biotables import CODON_USAGE
from ..biotools import (gc_content, reverse_complement, windows_overlap,
                        sequences_differences)



class ObjectiveEvaluation:
    """Store relevant infos about the evaluation of an objective on a canvas

    Examples
    --------

    >>> evaluation_result = ConstraintEvaluation(
    >>>     objective=objective,
    >>>     canvas = canvas,
    >>>     score= evaluation_score, # float
    >>>     windows=[(w1_start, w1_end), (w2_start, w2_end)...],
    >>>     message = "Score: 42 (found 42 sites)"
    >>> )

    Parameters
    ----------

    objective
      The Objective that was evaluated.

    canvas
      The canvas that the objective was evaluated on.

    score
      The score associated to the evaluation.

    windows
      A list of couples (start, end) indicating the windows on which the
      the optimization shoul be localized to improve the objective.

    message
      A message that will be returned by ``str(evaluation)``. It will notably
      be displayed by ``canvas.print_objectives_summaries``.

    """

    def __init__(self, objective, canvas, score, windows=None, message=None):
        self.objective = objective
        self.canvas = canvas
        self.score = score
        self.passes = score >= 0
        self.windows = windows
        self.message = message

    def __repr__(self):
        return (
            str((
                "Passes" if self.passes else "Fails",
                "score : %.02f" % self.score,
                self.windows
            )) +
            (self.message if (self.message is not None) else "")
        )

    def __str__(self):
        return (
            self.message if (self.message is not None) else
            "score: %.02E, windows: %s" % (self.score, str(self.windows))
        )


class Objective:
    """General class to define objectives to optimize.

    Note that all objective have a ``boost`` attribute that is a multiplicator
    that will be used when computing the global objective score of a canvas
    with ``canvas.all_objectives_score()``.

    New types of objectives are defined by subclassing ``Objective`` and
    providing a custom ``evaluate`` and ``localized`` methods.

    """

    best_possible_score = None
    can_be_solved_locally = False

    def __init__(self, evaluate, boost=1.0):
        self.boost = boost
        self.evaluate = evaluate

    # def evaluate(self, canvas):
    #     """Evaluate the objective on the provided canvas.
    #
    #     Returns a ``ObjectiveEvaluation`` object indicating whether the
    #     objective passed, the score , the windows on which to focus searches
    #     in case the objective failed and a string message
    #     summarizing the evaluation (see ``ObjectiveEvaluation``).
    #     """
    #     pass

    def localized(self, window):
        """Return a modified version of the objective for the case where
        sequence modifications are only performed inside the provided window.

        For instance if an objective concerns local GC content, and we are
        only making local mutations to destroy a restriction site, then we only
        need to check the local GC content around the restriction site after
        each mutation (and not compute it for the whole sequence), so
        ``EnforceGCContent.localized(window)`` will return an objective
        that only looks for GC content around the provided window.

        If an objective concerns a DNA segment that is completely disjoint from
        the provided window, this must return a ``VoidConstraint``.

        Must return an object of class ``Constraint``.
        """
        return self

    def copy_with_changes(self, **kwargs):
        new_objective = copy.deepcopy(self)
        new_objective.__dict__.update(kwargs)
        return new_objective

    def initialize_problem(self, problem, role="constraint"):
        return self



class VoidObjective(Objective):
    """Void Objectives are a special case of Objectives that always pass.

    Void Objectives are generally obtained when a Objective is "made void"
    by a localization. For instance, if we are optimizing the segment (10,50)
    of a DNA segment, the Objective EnforceTranslation([300,500]) does not
    apply as it concerns a Gene that is in a completely different segment.
    Therefore the localized version of EnforceTranslation will be void.
    """

    def __init__(self, parent_objective=None, boost=1.0):
        self.parent_objective = parent_objective
        self.message = ("Pass (not relevant in this context)"
                        % parent_objective)
        self.boost = boost

    def evaluate(self, canvas):
        """The evaluation of VoidObjectives always passes with score=1.0
        It returns a message indicating that the parent Objective was voided
        """
        return ObjectiveEvaluation(self, canvas, score=1.0,
                                   message=self.message,
                                   windows=None)

    def __repr__(self):
        return "Voided %s" % repr(self.parent_objective)

class PatternObjective(Objective):
    """Class for Objectives such as presence or absence of a pattern.

    The particularity of the PatternObjectives is that they will either infer
    or ask for the length of the associated pattern and use this to localize
    the objective efficiently when performing local optimization or solving.

    Parameters
    ----------

    pattern


    """
    can_be_solved_locally = False

    def __init__(self, pattern, window=None, boost=1.0):
        self.pattern = pattern
        self.window = window

    def localized(self, window):
        """Localize the pattern to the given window. Taking into account the
        objective's own window, and the size of the pattern."""
        pattern_size = self.pattern.size
        if self.window is not None:
            overlap = windows_overlap(self.window, window)
            if overlap is None:
                return VoidObjective(parent_objective=self)
            else:
                start, end = window
                ostart, oend = overlap
                if ostart == start:
                    new_window = [start, min(end, oend + pattern_size)]
                else:
                    new_window = [max(start, ostart - pattern_size), end]
        else:
            start, end = window
            margin = pattern_size - 1
            new_window = [max(0, start - margin), end + margin]

        return self.copy_with_changes(window=new_window)

class TerminalObjective(Objective):
    """Objectives that apply in the same way to both ends of the sequence.

    These are particularly useful for modeling constraints from providers
    who have terminal-ends constraints.

    Subclasses of these objectives should have a `window_size` and a
    `evaluate_end` method"""

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

        return ObjectiveEvaluation(self, canvas, score=len(windows),
                                   windows=windows, message=message)
