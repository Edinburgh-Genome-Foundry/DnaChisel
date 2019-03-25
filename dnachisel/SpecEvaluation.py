# -*- coding: utf-8 -*-
"""Classes for Objective Evaluations."""


import textwrap
import itertools
try:
    import matplotlib.cm as cm
    MATPLOTLIB_AVAILABLE = True
except:
    MATPLOTLIB_AVAILABLE = False

from Bio.SeqFeature import SeqFeature

from .Location import Location

def colors_cycle(lightness_factor=1.0, color_shift=0):
    if MATPLOTLIB_AVAILABLE:
        cycle = itertools.cycle([
            cm.Paired(color_shift + 0.21 * i % 1.0)
            for i in range(30)
        ])
        return (
            '#%02x%02x%02x' % tuple([int(255 * c * lightness_factor)
                                    for c in rgb_tuple[:3]])
            for rgb_tuple in cycle
        )
    else:
        return itertools.cycle([
            "#f1cccc", "#f1e5cc", "#e3f1cc", "#ccf1e3", "#ccd7f1", "#e0ccf1",
            "f1cce7"
        ])
    

def score_as_text(score):
    raw = str(int(score) if (int(score) == score) else score)
    as_float = "%.02f" % score
    as_eng = "%.02E." % score
    return min([raw, as_float, as_eng], key=len).rjust(9)

class SpecEvaluation:
    """Store relevant infos about the evaluation of an objective on a problem.

    Examples
    --------

    >>> evaluation_result = ConstraintEvaluation(
    >>>     specification=specification,
    >>>     problem = problem,
    >>>     score= evaluation_score, # float
    >>>     locations=[Locations list...],
    >>>     message = "Score: 42 (found 42 sites)"
    >>> )

    Parameters
    ----------

    objective
      The Objective that was evaluated.

    problem
      The problem that the objective was evaluated on.

    score
      The score associated to the evaluation.

    locations
      A list of couples (start, end) indicating the locations on which the
      the optimization shoul be localized to improve the objective.

    message
      A message that will be returned by ``str(evaluation)``. It will notably
      be displayed by ``problem.print_objectives_summaries``.

    """

    def __init__(self, specification, problem, score, locations=None,
                 message=None, data=None):
        """Initialize."""
        self.specification = specification
        self.problem = problem
        self.score = score
        self.passes = score >= 0
        self.is_optimal = (score == specification.best_possible_score)
        self.locations = locations
        self.message = self.default_message if message is None else message
        self.data = {} if data is None else data



    @property
    def default_message(self):
        """Return the default message for console/reports."""
        return "Score: %s. Locations: %s" % (score_as_text(self.score),
                                             self.locations)
    @property
    def score_as_text(self):
        return score_as_text(self.score)

    def to_text(self, role='constraint', wrapped=True):
        """Return a string representation of the evaluation.

        Parameters
        ----------
        role
          either 'objective' or 'constraint'

        wrapped
          Whether to wrap-and-indent the result

        Returns
        -------
        message
          A long string detailing the evaluations's status and regions of
          breach or suboptimality.

        """
        message = self.message
        if wrapped:
            indents = 6 if (role == 'constraint') else 11
            indent = indents*' ' + "│ "
            message = "\n".join(
                textwrap.wrap(message, width=80, initial_indent=indent,
                              subsequent_indent=indent)
            )

        if role == "objective":
            return "{optimal}{score} ┍ {spec} \n{message}".format(
                self=self, optimal="✔" if self.is_optimal else " ",
                spec=self.specification.label(with_location=True),
                score=score_as_text(self.score), message=message
            )
        else:
            return "{passes} ┍ {spec}\n{message}".format(
                self=self, passes="✔PASS" if self.passes else " FAIL",
                spec=self.specification.label(with_location=True),
                message=message)

    def locations_to_biopython_features(self, feature_type="misc_feature",
                                        color="red", label_prefix="",
                                        merge_overlapping=False):
        """Return a list of locations (of breach/suboptimality) as annotations.

        Parameters
        ----------
        feature_type
          Genbank type of the annotations

        color
          Color property attached to the annotations

        label_prefix
          The locations will be labelled of the form
          "prefix NameOfSpecification()"
        """
        locations = self.locations
        if merge_overlapping:
            locations = Location.merge_overlapping_locations(locations)
        return [
            SeqFeature(location.to_biopython_location(), type=feature_type,
                       qualifiers=dict(
                           label=label_prefix + " " + str(self.specification),
                           color=color
            ))
            for location in locations
        ]


class SpecEvaluations:
    """Base class for the collective handling of lists of SpecEvaluations.

    See ProblemObjectivesEvaluations and ProblemConstraintsEvaluations for
    the useful subclasses.

    Parameters
    ----------
    evaluations
      list of SpecEvaluations

    problem
      (optional) problem on which the evaluations were carried.

    """
    color_lightness = 1.0
    color_shift = 0

    def __init__(self, evaluations, problem=None):
        """Initialize."""
        self.evaluations = evaluations
        self.problem = problem

    def __iter__(self):
        """Iterate over evaluations."""
        return self.evaluations.__iter__()

    def __len__(self):
        """Return the number of evaluations."""
        return len(self.evaluations)

    def all_evaluations_pass(self):
        """Return whether all evaluations pass."""
        return all([ev.passes for ev in self.evaluations])

    def scores_sum(self, as_text=False):
        """Return the sum of all evaluations scores.

        Scores are multiplied by their respective boost factor.
        """
        result = sum([
            ev.specification.boost * ev.score
            for ev in self.evaluations
        ])
        if as_text:
            result = score_as_text(result)
        return result


    def filter(self, eval_filter):
        """Create a new instance with a subset of the evaluations.

        ``eval_filter`` is either a function ``(evaluation) => True/False``
        or one of "passing", "failing", "optimal", "suboptimal", to obtain the
        corresponding constraints.
        """
        if isinstance(eval_filter, str):
            eval_filter = {
                "passing": lambda e: e.passes,
                "failing": lambda e: not e.passes,
                "optimal": lambda e: e.is_optimal,
                "suboptimal": lambda e: not e.is_optimal
            }[eval_filter]
        return self.__class__(evaluations=[
            ev for ev in self.evaluations if eval_filter(ev)
        ], problem=self.problem)

    def to_text(self):
        """Return a long representation of the evaluations."""
        return "\n".join(["===> %s" % self.text_summary_message()] + [
            e.to_text(role=self.specifications_role)
            for e in self.evaluations
        ]) + "\n\n"

    def evaluations_with_locations(self):
        """Return the list of all evaluations whose location is not None."""
        return [
            ev for ev in self.evaluations
            if ev.locations is not None
        ]

    def success_and_failures_as_features(self, feature_type="misc_feature"):
        """Return all evaluations as Biopython features.

        With color property depending on whether the evaluation is passing,
        failing, optimal, or suboptimal, the color being determined by
        ``self.success_failure_color``.

        """
        return [
            ev.specification.to_biopython_feature(
                feature_type=feature_type,
                color=self.success_failure_color(ev),
                passes='true' if ev.passes else 'false',
                is_optimal='true' if ev.is_optimal else 'false',
            )
            for ev in self.evaluations
            if ev.specification.__dict__.get('location', False)
        ]

    def locations_as_features(self, features_type="misc_feature",
                              with_specifications=True, label_prefix="From",
                              colors="cycle", merge_overlapping=False):
        """Return all locations from all evaluations as biopython features.

        Parameters
        ----------
        features_type
          The Genbank feature type

        with_specifications
          If True, features are added to the list to indicate the scope of the
          different specifications. If false, only the specification breaches
          are returned.

        label_prefix
          Each breach may be labeled "prefix NameOfSpec(props)", for instance,
          "From AvoidPattern(100-200)", to indicate where the breach belongs.

        colors
          Either a list of colors (one for each specification) or "cycle"
          for cycling through predefined colors. The colors are applied to all
          breaches.

        """
        if colors == "cycle":
            cycle = colors_cycle(lightness_factor=self.color_lightness,
                                 color_shift=self.color_shift)
            colors = [next(cycle) for ev in self.evaluations]
        features = [
            location.to_biopython_feature(
                feature_type="misc_feature",
                specification=label_prefix + " " + str(ev.specification),
                color=color
            )
            for (ev, color) in zip(self.evaluations_with_locations(), colors)
            for location in (
                Location.merge_overlapping_locations(ev.locations)
                if merge_overlapping else ev.locations
            )
        ]
        if with_specifications:
            features += [
                ev.specification.to_biopython_feature(
                    feature_type="misc_feature",
                    label=str(ev.specification),
                    role=self.specifications_role,
                    color=color
                )
                for ev, color in zip(self.evaluations, colors)
                if ev.specification.__dict__.get('location', False)
            ]
        return features


class ProblemConstraintsEvaluations(SpecEvaluations):
    """Special multi-evaluation class for all constraints of a same problem.

    See submethod ``.from_problem``

    """

    specifications_role = "constraint"

    @staticmethod
    def from_problem(problem):
        """Create an instance by evaluating all constraints in the problem.

        The ``problem`` is a DnaChisel DnaOptimizationProblem.

        """
        return ProblemConstraintsEvaluations([
            specification.evaluate(problem)
            for specification in problem.constraints
        ], problem=problem)

    def success_failure_color(self, evaluation):
        """Return color #60f979 if evaluation.passes else #f96c60."""
        return "#60f979" if evaluation.passes else "#f96c60"

    def text_summary_message(self):
        """Return a global SUCCESS or FAILURE message for all evaluations."""
        failed = [e for e in self.evaluations if not e.passes]
        if failed == []:
            return "SUCCESS - all constraints evaluations pass"
        else:
            return "FAILURE: %d constraints evaluations failed" % len(failed)

class ProblemObjectivesEvaluations(SpecEvaluations):
    """Special multi-evaluation class for all objectives of a same problem.

    See submethod ``.from_problem``

    """
    color_lightness = 0.8
    color_shift = 0.14

    specifications_role = "objective"

    @staticmethod
    def from_problem(problem):
        """Create an instance by evaluating all objectives in the problem.

        The ``problem`` is a DnaChisel DnaOptimizationProblem.

        """
        return ProblemObjectivesEvaluations([
            specification.evaluate(problem)
            for specification in problem.objectives
        ], problem=problem)

    def success_failure_color(self, evaluation):
        """Return color #cbf960 if evaluation is optimal else #f9a260."""
        return "#cbf960" if evaluation.is_optimal else "#f9a260"

    def text_summary_message(self):
        """Return a TOTAL SCORE message."""
        if len(self.evaluations) == 0:
            return "No specifications"
        return "TOTAL OBJECTIVES SCORE: " + self.scores_sum(as_text=True)
