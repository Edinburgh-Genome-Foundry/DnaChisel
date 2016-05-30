from biotables import CODON_USAGE
from biotools import gc_content, reverse_complement, sequences_differences
import numpy as np


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
      The canvas that the constraint was evaluated on.

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
        self.windows = windows
        self.message = message

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

    max_possible_score = None

    def __init__(self, boost=1. - 0):
        self.boost = boost

    def evaluate(self):
        """"""
        pass

    def localized(self, window):
        return self


class CodonOptimizationObjective(Objective):
    """Objective to codon-optimize a coding sequence for a particular organism.

    Several codon-optimization policies exist. At the moment this Objective
    implements a method in which codons are replaced by the most frequent
    codon in the species.

    (as long as this doesn't break any constraint or lowers the global
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
                 codon_usage_table=None, boost=1.0,):
        Objective.__init__(self, boost=boost)
        self.window = window
        self.strand = strand
        self.organism = organism
        if organism is not None:
            codon_usage_table =  CODON_USAGE[self.organism]
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

def subdivide_window(window, max_span):
    start, end = window
    inds = list(range(start, end, max_span))+[end]
    return zip(inds, inds[1:])

class GCContentObjective(Objective):
    """Objective to obtain the desired global GC content.

    Parameters
    ----------

    gc_content
      GC content to target (e.g. ``0.4``)

    window
      Restrict the GC content observation to a DNA segment (start, end)

    """

    max_possible_score = 0

    def __init__(self, gc_content, exponent=1.0, window=None,  boost=1.0,
                 subdivision_window=None):
        Objective.__init__(self, boost=boost)
        self.gc_content = gc_content
        self.exponent = exponent
        self.window = window
        self.subdivision_window = subdivision_window

    def evaluate(self, canvas):
        window = (self.window if self.window is not None
                  else [0, len(canvas.sequence)])
        start, end = window
        subsequence = canvas.sequence[start: end]
        gc = gc_content(subsequence)
        score = -(abs(gc - self.gc_content) ** self.exponent)
        if self.subdivision_window is None:
            windows = [window]
        else:
            windows = subdivide_window(window, self.subdivision_window)

        return ObjectiveEvaluation(
            self, canvas, score=score, windows=windows,
            message="scored %.02E. GC content is %.03f (%.03f wanted)" %
                    (score, gc, self.gc_content))

    def __str__(self):
        return "GCContentObj(%.02f, %s)" % (
            self.gc_content,
            "global" if (self.window is None) else
            ("window: %s" % str(self.window))
        )


class MinimizeDifferencesObjective(Objective):
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
    >>>     constraints = [GCContentConstraint(0.3,0.6, gc_window=50)],
    >>>     objective = [MinimizeDifferencesObjective(
    >>>                     original_sequence=sequence)]
    >>> )
    >>> canvas.solve_all_constraints_one_by_one()
    >>> canvas.maximize_all_objectives_one_by_one()

    """

    max_possible_score = 0

    def __init__(self, window=None, target_sequence=None,
                 original_sequence=None, boost=1.0):
        Objective.__init__(self, boost=boost)
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
