from .SpecEvaluation import SpecEvaluation
from .SpecEvaluations import SpecEvaluations


class ProblemConstraintsEvaluations(SpecEvaluations):
    """Special multi-evaluation class for all constraints of a same problem.

    See submethod ``.from_problem``

    """

    specifications_role = "constraint"

    @staticmethod
    def from_problem(problem, autopass_constraints=True):
        """Create an instance by evaluating all constraints in the problem.

        The ``problem`` is a DnaChisel DnaOptimizationProblem.

        """

        def evaluate(constraint):
            if (
                autopass_constraints
                and constraint.enforced_by_nucleotide_restrictions
            ):
                return SpecEvaluation(
                    constraint,
                    problem,
                    score=1,
                    locations=[],
                    message="Enforced by nucleotides restrictions",
                )
            else:
                return constraint.evaluate(problem)

        return ProblemConstraintsEvaluations(
            [evaluate(constraint) for constraint in problem.constraints],
            problem=problem,
        )

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
