# -*- coding: utf-8 -*-
"""Implements the SpecEvaluation class

SpecEvaluation is a class describing the result of an the evaluation of a
*DnaOptimizationProblem* by a *Specification*. It contains a score, a message,
a list of *Locations* of sub-optimal regions.

Several evaluations can be grouped using classes
*ProblemConstraintsEvaluations* and *ProblemObjectivesEvaluations*, which
implement methods for printing or exporting as Genbank a set of evaluations.
"""

import textwrap

from Bio.SeqFeature import SeqFeature
from ...biotools import score_to_formatted_string
from ...Location import Location


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

    def __init__(
        self,
        specification,
        problem,
        score,
        locations=None,
        message=None,
        data=None,
    ):
        """Initialize."""
        self.specification = specification
        self.problem = problem
        self.score = score
        self.passes = score >= 0
        self.is_optimal = score == specification.best_possible_score
        self.locations = locations
        self.message = self.default_message if message is None else message
        self.data = {} if data is None else data

    @property
    def default_message(self):
        """Return the default message for console/reports."""
        return "Score: %s. Locations: %s" % (
            score_to_formatted_string(self.score),
            self.locations,
        )

    @property
    def score_to_formatted_string(self):
        return score_to_formatted_string(self.score)

    def to_text(self, role="constraint", wrapped=True, max_message_length=500):
        """Return a string representation of the evaluation.

        Example output for a constraint:

        >>> FAIL ┍ UniquifyAllKmers[10-1000](k:9)
        >>>      │ Score:        -2. Locations: [232-241, 233-242]

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
        if len(message) > max_message_length:
            half = int(max_message_length / 2)
            message = message[:half] + " ... " + message[-half:]
        if wrapped:
            indents = 6 if (role == "constraint") else 11
            indent = indents * " " + "│ "
            message = "\n".join(
                textwrap.wrap(
                    message,
                    width=80,
                    initial_indent=indent,
                    subsequent_indent=indent,
                )
            )

        if role == "objective":
            return "{optimal}{score} ┍ {spec} \n{message}".format(
                self=self,
                optimal="✔" if self.is_optimal else " ",
                spec=self.specification.label(with_location=True),
                score=score_to_formatted_string(self.score),
                message=message,
            )
        else:
            return "{passes} ┍ {spec}\n{message}".format(
                self=self,
                passes="✔PASS" if self.passes else " FAIL",
                spec=self.specification.label(with_location=True),
                message=message,
            )

    def locations_to_biopython_features(
        self,
        feature_type="misc_feature",
        color="red",
        label_prefix="",
        merge_overlapping=False,
    ):
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

        merge_overlapping
          If true, then overlapping locations (0-5, 2-9) will be merged into a
          single one (0-9).
        """
        locations = self.locations
        if merge_overlapping:
            locations = Location.merge_overlapping_locations(locations)
        return [
            SeqFeature(
                location.to_biopython_location(),
                type=feature_type,
                qualifiers=dict(
                    label=label_prefix + " " + str(self.specification),
                    color=color,
                ),
            )
            for location in locations
        ]
