"""Attempt at a special subclass of DnaOptimizationProblem for circular
sequences. Needs more docs. See example in the examples folder.
"""
from ..reports.optimization_reports import (
    write_optimization_report,
    write_no_solution_report,
)
from .DnaOptimizationProblem import DnaOptimizationProblem
from .NoSolutionError import NoSolutionError
from ..Location import Location


class CircularViewProblem(DnaOptimizationProblem):
    """Class representing an optimization problem as the concatenation of 3
    times the same sequence and specifications, in order to model the
    circularity of the DNA sequence"""

    def _replace_sequence(self, new_sequence):
        L = len(new_sequence) // 3

        def return_the_loony(a, b, c):
            """Return the element of a, b, c, which is unlike the other 2"""
            if a == b:
                return c
            elif a == c:
                return b
            else:
                return a

        self.sequence = 3 * "".join(
            [
                return_the_loony(
                    new_sequence[i],
                    new_sequence[i + L],
                    new_sequence[i + 2 * L],
                )
                for i in range(L)
            ]
        )


class CircularDnaOptimizationProblem(DnaOptimizationProblem):
    """Class for solving circular DNA optimization problems.

    The high-level interface is the same as for DnaOptimizationProblem:
    initialization, resolve_constraints(), and optimize() work the same way.
    """

    def _circularized_specs(self, specs, central_specs_only=False):
        L = len(self.sequence)
        new_specs = []
        for spec in specs:
            loc = spec.location
            if (loc.start, loc.end) == (0, L):
                new_location = Location(0, 3 * L, strand=loc.strand)
                new_specs.append(
                    spec.copy_with_changes(
                        location=new_location, derived_from=spec
                    )
                )
            else:
                new_location = loc + L
                new_specs += [spec.shifted(i) for i in [0, L, 2 * L]]
        central_loc = Location(L, 2 * L)
        if central_specs_only:
            new_specs = [
                spec
                for spec in new_specs
                if spec.location.overlap_region(central_loc) is not None
            ]
        return new_specs

    def _circularized_view(
        self,
        with_objectives=False,
        with_constraints=False,
        central_specs_only=False,
    ):

        return CircularViewProblem(
            sequence=3 * self.sequence,
            constraints=self._circularized_specs(
                self.constraints, central_specs_only=central_specs_only
            )
            if with_constraints
            else [],
            objectives=self._circularized_specs(
                self.objectives, central_specs_only=central_specs_only
            )
            if with_objectives
            else [],
            logger=self.logger,
        )

    def _recentered_evaluations(self, evals):
        L = len(self.sequence)
        central_loc = Location(L, 2 * L)
        for e in list(evals.evaluations):
            e.specification = e.specification.derived_from
            if e.locations is not None:
                e.locations = [
                    (loc - L)
                    for loc in e.locations
                    if loc.overlap_region(central_loc) is not None
                ]
                e.message = e.default_message
        return evals

    def constraints_evaluations(self, autopass=True):
        """Return the evaluation of constraints.

        The "autopass_constraints" enables to just assume that constraints
        enforced by the mutation space are verified.

        """
        circularized = self._circularized_view(
            with_constraints=True, central_specs_only=False
        )
        evals = circularized.constraints_evaluations(autopass=autopass)
        return self._recentered_evaluations(evals)
    
    def all_constraints_pass(self, autopass=True):
        """Return whether the current problem sequence passes all constraints.
        """
        evals = self.constraints_evaluations(autopass=autopass)
        return evals.all_evaluations_pass()

    def objectives_evaluations(self):
        circularized = self._circularized_view(
            with_objectives=True, central_specs_only=False
        )
        evals = circularized.objectives_evaluations()
        return self._recentered_evaluations(evals)

    def resolve_constraints(self, final_check=True):
        problem = self._circularized_view(with_constraints=True)
        L = len(self.sequence)
        central_loc = Location(L, 2 * L)
        for c in problem.constraints:
            if c.location.overlap_region(central_loc) is not None:
                problem.resolve_constraint(c)

        self.sequence = problem.sequence[L : 2 * L]
        if final_check:
            self.perform_final_constraints_check()

    def optimize(self):
        problem = self._circularized_view(
            with_constraints=True, with_objectives=True
        )
        problem.optimize()
        L = len(self.sequence)
        self.sequence = problem.sequence[L : 2 * L]

    def optimize_with_report(
        self,
        target,
        project_name="My project",
        file_path=None,
        file_content=None,
    ):
        """Resolve constraints, optimize objectives, write a multi-file report.

        WARNING: in case of optimization failure, the report generated will
        show a "pseudo-circular" sequence formed by concatenating the sequence
        with itself three times.

        TODO: fix the warning above, at some point?

        The report's content may vary depending on the optimization's success.

        Parameters
        ----------

        target
          Either a path to a folder that will containt the report, or a path to
          a zip archive, or "@memory" to return raw data of a zip archive
          containing the report.

        project_name
          Project name to write on PDF reports

        Returns
        -------

        (success, message, zip_data)
          Triplet where success is True/False, message is a one-line string
          summary indication whether some clash was found, or some solution, or
          maybe no solution was found because the random searches were too
          short
        """
        self.logger(message="Solving constraints")
        try:
            self.resolve_constraints()
        except NoSolutionError as error:
            problem = self._circularized_view(
                with_constraints=True, with_objectives=True
            )
            self.logger(message="No solution found: making report")
            data = write_no_solution_report(
                target,
                problem,
                error,
                file_path=file_path,
                file_content=file_content,
            )
            start, end, s = error.location.to_tuple()
            message = "No solution found in zone [%d, %d]: %s." % (
                start,
                end,
                str(error),
            )
            return False, message, data
        self.logger(message="Now optimizing the sequence")
        self.optimize()
        self.logger(message="Success! Generating report.")
        data = write_optimization_report(
            target,
            self,
            project_name=project_name,
            file_path=file_path,
            file_content=file_content,
        )
        return True, "Optimization successful.", data
