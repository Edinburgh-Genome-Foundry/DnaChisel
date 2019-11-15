from .SpecEvaluations import SpecEvaluations


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
        return ProblemObjectivesEvaluations(
            [
                specification.evaluate(problem)
                for specification in problem.objectives
            ],
            problem=problem,
        )

    def success_failure_color(self, evaluation):
        """Return color #cbf960 if evaluation is optimal else #f9a260."""
        return "#cbf960" if evaluation.is_optimal else "#f9a260"

    def text_summary_message(self):
        """Return a TOTAL SCORE message."""
        if len(self.evaluations) == 0:
            return "No specifications"
        return "TOTAL OBJECTIVES SCORE: " + self.scores_sum(as_text=True)
