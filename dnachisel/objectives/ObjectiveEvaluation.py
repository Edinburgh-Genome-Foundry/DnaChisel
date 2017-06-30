"""Classes for Objective Evaluations."""

from Bio.SeqFeature import SeqFeature
from ..plotting_tools import colors_cycle


class ObjectiveEvaluation:
    """Store relevant infos about the evaluation of an objective on a problem.

    Examples
    --------

    >>> evaluation_result = ConstraintEvaluation(
    >>>     objective=objective,
    >>>     problem = problem,
    >>>     score= evaluation_score, # float
    >>>     locations=[(w1_start, w1_end), (w2_start, w2_end)...],
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

    def __init__(self, objective, problem, score, locations=None,
                 message=None):
        """Initialize."""
        self.objective = objective
        self.problem = problem
        self.score = score
        self.passes = score >= 0
        self.is_optimal = (score == objective.best_possible_score)
        self.locations = locations
        self.message = self.default_message if message is None else message

    @property
    def default_message(self):
        """Return the default message for console/reports."""
        return "Score: %.02E. Locations: %s" % (self.score, self.locations)

    def to_text(self, role=None):
        """Return a string representation of the evaluation."""
        if role == "objective":
            return ("{optimal} Scored {self.score:.02E} | {self.objective} | "
                    "{self.message}").format(
                self=self, optimal="OPTIMAL  | " if self.is_optimal else ""
            )
        else:
            return "{passes} | {self.objective} | {self.message}".format(
                self=self, passes="PASS" if self.passes else "FAIL")

    def locations_to_biopython_features(self, feature_type="misc_feature",
                                        color="red", label_prefix=None):
        return [
            SeqFeature(location.to_biopython_location(), type=feature_type,
                       qualifiers=dict(
                           label=label_prefix + " " + str(self.objective),
                           color=color
            ))
            for location in self.locations
        ]


class ObjectiveEvaluations:

    def __init__(self, evaluations, problem=None):
        self.evaluations = evaluations
        self.problem = problem

    def __iter__(self):
        return self.evaluations.__iter__()

    def __len__(self):
        return len(self.evaluations)

    def all_evaluations_pass(self):
        return all([ev.passes for ev in self.evaluations])

    def scores_sum(self):
        return sum([
            ev.objective.boost * ev.score
            for ev in self.evaluations
        ])

    def filter(self, eval_filter):
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
        return "\n".join(["===> %s" % self.text_summary_message()] + [
            e.to_text(role=self.objectives_role)
            for e in self.evaluations
        ]) + "\n\n"

    def evaluations_with_locations(self):
        return [
            ev for ev in self.evaluations
            if ev.locations is not None
        ]

    def success_and_failures_as_features(self, feature_type="misc_feature"):
        return [
            ev.objective.to_biopython_feature(
                feature_type=feature_type,
                color=self.success_failure_color(ev),
                passes='true' if ev.passes else 'false',
                is_optimal='true' if ev.is_optimal else 'false',
            )
            for ev in self.evaluations
            if ev.objective.__dict__.get('location', False)
        ]

    def locations_as_features(self, features_type="misc_feature",
                              with_objectives=True, label_prefix="From",
                              colors="cycle"):
        if colors == "cycle":
            cycle = colors_cycle()
            colors = [next(cycle) for ev in self.evaluations]

        features = [
            location.to_biopython_feature(
                feature_type="misc_feature",
                objective=label_prefix + " " + str(ev.objective),
                color=color
            )
            for (ev, color) in zip(self.evaluations_with_locations(), colors)
            for location in ev.locations
        ]
        if with_objectives:
            features += [
                ev.objective.to_biopython_feature(
                    feature_type="misc_feature",
                    label=str(ev.objective),
                    role=self.objectives_role,
                    color=color
                )
                for ev, color in zip(self.evaluations, colors)
                if ev.objective.__dict__.get('location', False)
            ]
        return features


class ProblemConstraintsEvaluations(ObjectiveEvaluations):
    objectives_role = "constraint"

    @staticmethod
    def from_problem(problem):
        return ProblemConstraintsEvaluations([
            objective.evaluate(problem)
            for objective in problem.constraints
        ], problem=problem)

    def success_failure_color(self, evaluation):
        return "#60f979" if evaluation.passes else "#f96c60"

    def text_summary_message(self):
        failed = [e for e in self.evaluations if not e.passes]
        if failed == []:
            return "SUCCESS - all constraints evaluations pass"
        else:
            return "FAILURE: %d constraints evaluations failed" % len(failed)


class ProblemObjectivesEvaluations(ObjectiveEvaluations):
    objectives_role = "objective"

    @staticmethod
    def from_problem(problem):
        return ProblemObjectivesEvaluations([
            objective.evaluate(problem)
            for objective in problem.objectives
        ], problem=problem)

    def success_failure_color(self, evaluation):
        return "#cbf960" if evaluation.is_optimal else "#f9a260"

    def text_summary_message(self):
        return "TOTAL OBJECTIVES SCORE: %.02f" % self.scores_sum()
