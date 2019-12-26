from ...Location import Location
from ...Specification.SpecEvaluation import (
    ProblemConstraintsEvaluations,
)
from ..NoSolutionError import NoSolutionError


class ConstraintsSolverMixin:
    @property
    def constraints_before(self):
        """"""
        if self._constraints_before is None:
            sequence = self.sequence
            self.sequence = self.sequence_before
            self._constraints_before = self.constraints_evaluations()
            self.sequence = sequence
        return self._constraints_before

    def constraints_evaluations(self, autopass=True):
        """Return a list of the evaluations of each constraint of the problem.

        The "autopass" enables to just assume that constraints
        enforced by the mutation space are verified.
        """
        return ProblemConstraintsEvaluations.from_problem(
            self, autopass_constraints=autopass
        )

    def all_constraints_pass(self, autopass=True):
        """Return whether the current problem sequence passes all constraints.
        """
        if len(self.constraints) == 0:
            return True
        return all(
            c.evaluate(self).passes
            for c in self.constraints
            if (not autopass) or (not c.enforced_by_nucleotide_restrictions)
        )

    def constraints_text_summary(self, failed_only=False, autopass=True):
        evals = self.constraints_evaluations(autopass=autopass)
        if failed_only:
            evals = evals.filter("failing")
        return evals.to_text()

    def get_focus_constraint(self):
        focus_constraints = [c for c in self.constraints if c.is_focus]
        if len(focus_constraints) == 1:
            focus = focus_constraints[0]
            other_constraints = [c for c in self.constraints if not c.is_focus]
            return focus, other_constraints
        return None, None

    def resolve_constraints_by_exhaustive_search(self):
        """Solve all constraints by exploring the whole search space.

        This method iterates over ``self.iter_mutations_space()`` (space of
        all sequences that could be reached through successive mutations) and
        stops when it finds a sequence which meets all the constraints of the
        problem.
        """
        focus_constraint, other_constraints = self.get_focus_constraint()
        sequence_before = self.sequence
        all_variants = self.mutation_space.all_variants(self.sequence)
        space_size = int(self.mutation_space.space_size)
        self.logger(mutation__total=space_size)
        for variant in self.logger.iter_bar(mutation=all_variants):
            self.sequence = variant
            if focus_constraint is not None:
                if focus_constraint.evaluate(self).passes:
                    if all(c.evaluate(self).passes for c in other_constraints):
                        self.logger(mutation__index=space_size)
                        return
            elif self.all_constraints_pass():
                self.logger(mutation__index=space_size)
                return
        self.sequence = sequence_before
        raise NoSolutionError(
            "Exhaustive search failed to satisfy all constraints.",
            problem=self,
        )

    def resolve_constraints_by_random_mutations(self):
        """Solve all constraints by successive sets of random mutations.

        This method modifies the problem sequence by applying a number
        ``mutations_per_iteration`` of random mutations. The constraints are
        then evaluated on the new sequence. If all constraints pass, the new
        sequence becomes the problem's new sequence.
        If not all constraints pass, the sum of all scores from failing
        constraints is considered. If this score is superior to the score of
        the previous sequence, the new sequence becomes the problem's new
        sequence.

        This operation is repeated `max_iter` times at most, after which
        a ``NoSolutionError`` is thrown if no solution was found.
        """

        focus_constraint, other_constraints = self.get_focus_constraint()
        if focus_constraint is not None:
            self.resolve_single_constraint_by_random_mutations(
                focus_constraint, other_constraints
            )
            return

        evaluations = self.constraints_evaluations()
        score = sum([e.score for e in evaluations if not e.passes])
        iters = range(self.max_random_iters)
        for i in self.logger.iter_bar(mutation=iters):

            if all(e.passes for e in evaluations):
                self.logger(mutation__index=iters)
                return
            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=self.mutations_per_iteration,
                sequence=self.sequence,
            )

            evaluations = self.constraints_evaluations()
            new_score = sum([e.score for e in evaluations if not e.passes])

            if new_score > score:
                score = new_score
            else:
                self.sequence = previous_sequence
        raise NoSolutionError(
            "Random search did not find a solution in the given number of "
            "attempts. Try to increase the number of attempts with:\n\n"
            "problem.max_random_iters = 5000 # or even 10000, 20000, etc.\n\n"
            "If the problem persists, you may be in presence of a complex or "
            "unsolvable problem.",
            problem=self,
        )

    def resolve_single_constraint_by_random_mutations(
        self, constraint, other_constraints
    ):
        evaluation = constraint.evaluation
        score = evaluation.score
        iters = range(self.max_random_iters)
        for i in self.logger.iter_bar(mutation=iters):
            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=self.mutations_per_iteration,
                sequence=self.sequence,
            )
            new_evaluation = constraint.evaluate(self)
            if new_evaluation.score > score:
                if all(c.evaluate(self).passes for c in other_constraints):
                    score = new_evaluation.score
                    if new_evaluation.passes:
                        self.logger(mutation__index=iters)
                        return
                else:
                    self.sequence = previous_sequence
            else:
                self.sequence = previous_sequence

        raise NoSolutionError(
            "Random search did not find a solution in the given number of "
            "attempts. Try to increase the number of attempts with:\n\n"
            "problem.max_random_iters = 5000 # or even 10000, 20000, etc.\n\n"
            "If the problem persists, you may be in presence of a complex or "
            "unsolvable problem.",
            problem=self,
        )

    def resolve_constraints_locally(self):
        """Perform a local search, either stochastic or exhaustive.
        """
        if self.mutation_space.space_size < self.randomization_threshold:
            self.resolve_constraints_by_exhaustive_search()
        else:
            self.resolve_constraints_by_random_mutations()

    def resolve_constraint(self, constraint):
        """Resolve a constraint through successive localizations."""

        # EVALUATE THE CONSTRAINT, FIND BREACHING LOCATIONS

        evaluation = constraint.evaluate(self)
        if evaluation.passes:
            return

        locations = sorted(evaluation.locations)
        iterator = self.logger.iter_bar(
            location=locations, bar_message=lambda loc: str(loc)
        )

        # FOR EACH LOCATION, CREATE A LOCAL PROBLEM AND RESOLVE LOCALLY.

        for i, location in enumerate(iterator):

            # SEVERAL "EXTENSIONS" OF THE LOCAL ZONE WILL BE TESTED IN TURN
            # IN CASE THE LOCAL SEQUENCE IS FROZEN DUE TO NUCLEOTIDE INTER-
            # DEPENDENCIES (CODONS, ETC.)

            for extension in self.local_extensions:
                new_location = location.extended(extension)
                mutation_space = self.mutation_space.localized(new_location)

                if mutation_space.space_size == 0:

                    # If the sequence is frozen at this location, either
                    # "continue" (go straight to the next, larger extension)
                    # or if we are already in the largest extension, return
                    # an error with data that will be used by the report
                    # generator.

                    if extension != self.local_extensions[-1]:
                        continue
                    else:
                        error = NoSolutionError(
                            location=new_location,
                            problem=self,
                            message="Constraint breach in region that cannot "
                            "be mutated.",
                        )
                        error.location = new_location
                        error.constraint = constraint
                        error.message = "While solving %s in %s:\n\n%s" % (
                            constraint,
                            new_location,
                            str(error),
                        )
                        self.logger(
                            location__index=len(locations),
                            location__message="Cold exit",
                        )
                        raise error
                new_location = Location(*mutation_space.choices_span)

                # This blocks solves the problem of overlapping breaches,
                # which can make the local optimization impossible.
                # If the next constraint breach overlaps with the current
                # location, localize the constraint with a with_righthand=False
                # flag, which will be used by the constraints ".localized"
                # method to only consider the right-hand side.

                if (i < (len(locations) - 1)) and (
                    locations[i + 1].overlap_region(new_location)
                ):
                    this_local_constraint = constraint.localized(
                        new_location, with_righthand=False, problem=self
                    )
                else:
                    this_local_constraint = constraint.localized(
                        new_location, problem=self
                    )
                evaluation = this_local_constraint.evaluate(self)

                # MAYBE THE LOCAL BREACH WAS ALREADY RESOLVED AS A SIDE EFFECT
                # OF SOLVING PREVIOUS BREACHES. IN THAT CASE, PASS.

                if evaluation.passes:
                    continue

                # ELSE, CREATE A NEW LOCAL PROBLEM WITH LOCALIZED CONSTRAINTS

                this_local_constraint.is_focus = True
                this_local_constraint.evaluation = evaluation

                localized_constraints = [
                    cst.localized(new_location, problem=self)
                    for cst in self.constraints
                    if cst != constraint
                    and not cst.enforced_by_nucleotide_restrictions
                ]
                passing_localized_constraints = [
                    cst
                    for cst in localized_constraints
                    if cst is not None and cst.evaluate(self).passes
                ]
                local_problem = self.__class__(
                    sequence=self.sequence,
                    constraints=(
                        [this_local_constraint] + passing_localized_constraints
                    ),
                    mutation_space=mutation_space,
                )
                local_problem.randomization_threshold = (
                    self.randomization_threshold
                )
                local_problem.max_random_iters = self.max_random_iters
                local_problem.mutations_per_iteration = (
                    self.mutations_per_iteration
                )

                # STORE THE LOCAL PROBLEM IN THE LOGGER.
                # This is useful for troubleshooting.

                self.logger.store(
                    problem=self,
                    local_problem=local_problem,
                    location=location,
                )

                # RESOLVE THE LOCAL PROBLEM. RETURN AN ERROR IF IT FAILS.

                try:
                    if hasattr(constraint, "resolution_heuristic"):
                        constraint.resolution_heuristic(local_problem)
                    else:
                        local_problem.resolve_constraints_locally()
                    self._replace_sequence(local_problem.sequence)
                    break
                except NoSolutionError as error:
                    if extension == self.local_extensions[-1]:
                        error.location = new_location
                        error.constraint = constraint
                        error.message = "While solving %s in %s:\n\n%s" % (
                            constraint,
                            new_location,
                            str(error),
                        )
                        self.logger(
                            location__index=len(locations),
                            location__message="Cold exit",
                        )
                        raise error
                    else:
                        continue

    def resolve_constraints(self, final_check=True, cst_filter=None):
        """Solve a particular constraint using local, targeted searches.

        Parameters
        ----------

        constraint
          The ``Specification`` object for which the sequence should be solved

        final_check
          If True, a final check of that all constraints pass will be run at
          the end of the process, when constraints have been resolved one by
          one, to check that the solving of one constraint didn't undo the
          solving of another.

        cst_filter
          An optional filter to only resolve a subset function (constraint => True/False)

        """
        constraints = [
            c
            for c in self.constraints
            if not c.enforced_by_nucleotide_restrictions
            and ((cst_filter is None) or cst_filter(c))
        ]
        if len(constraints) == 0:
            return
        constraints = sorted(constraints, key=lambda c: -c.priority)
        for constraint in self.logger.iter_bar(
            constraint=constraints, bar_message=lambda c: str(c)
        ):
            try:
                self.resolve_constraint(constraint=constraint)
            except NoSolutionError as error:
                self.logger(constraint__index=len(constraints))
                raise error
        if final_check:
            self.perform_final_constraints_check()

    def perform_final_constraints_check(self):
        for cst in self.constraints:
            if not cst.evaluate(self).passes:
                raise NoSolutionError(
                    "The solving of all constraints failed to solve"
                    " all constraints, as some appear unsolved at the end"
                    " of the optimization. This is an unintended behavior,"
                    " likely due to a complex problem. Try running the"
                    " solver on the same sequence again, or report the"
                    " error to the maintainers:\n\n"
                    + self.constraints_text_summary(
                        failed_only=True, autopass=False
                    ),
                    problem=self,
                )
