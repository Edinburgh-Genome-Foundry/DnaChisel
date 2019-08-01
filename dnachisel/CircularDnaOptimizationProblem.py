"""Attempt at a special subclass of DnaOptimizationProblem for circular
sequences. Needs more docs. See example in the examples folder.
"""

from .DnaOptimizationProblem import DnaOptimizationProblem
from .Location import Location


class CircularViewProblem(DnaOptimizationProblem):
    """Class representing an optimization problem as the concatenation of 3
    times the same sequence and specifications, in order to model the
    circularity of the DNA sequence"""
    def change_sequence(self, new_sequence):
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

    def constraints_evaluations(self):
        circularized = self._circularized_view(
            with_constraints=True, central_specs_only=False
        )
        evals = circularized.constraints_evaluations()
        return self._recentered_evaluations(evals)

    def objectives_evaluations(self):
        circularized = self._circularized_view(
            self, with_objectives=True, central_specs_only=False
        )
        evals = circularized.objectives_evaluations()
        return self._recentered_evaluations(evals)

    def resolve_constraints(self, final_check=True):
        problem = self._circularized_view(with_constraints=True)
        problem.resolve_constraints(final_check=final_check)
        L = len(self.sequence)
        self.sequence = problem.sequence[L : 2 * L]

    def optimize(self):
        problem = self._circularized_view(
            self, with_constraints=True, with_objectives=True
        )
        problem.optimize()
        L = len(self.sequence)
        self.sequence = problem.sequence[L : 2 * L]

    def optimize_with_report(self, target, project_name="My project"):
        problem = self._circularized_view(
            with_constraints=True, with_objectives=True
        )
        return problem.optimize_with_report(target, project_name=project_name)

