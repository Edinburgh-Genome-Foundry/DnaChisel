from ...Location import Location
from ...Specification.SpecEvaluation import ProblemObjectivesEvaluations
from ..NoSolutionError import NoSolutionError


class ObjectivesMaximizerMixin:
    @property
    def objectives_before(self):
        if self._objectives_before is None:
            sequence = self.sequence
            self.sequence = self.sequence_before
            self._objectives_before = self.objectives_evaluations()
            self.sequence = sequence
        return self._objectives_before

    def objectives_evaluations(self):
        """Return a list of the evaluation of each objective of the problem"""
        return ProblemObjectivesEvaluations.from_problem(self)

    def objective_scores_sum(self):
        return self.objectives_evaluations().scores_sum()

    def objectives_text_summary(self):
        return self.objectives_evaluations().to_text()

    def optimize_by_exhaustive_search(self):
        """
        """
        if not self.all_constraints_pass():
            summary = self.constraints_text_summary(failed_only=True)
            raise NoSolutionError(
                summary
                + "Optimization can only be done when all constraints are "
                "verified.",
                self,
            )

        if all(
            [obj.best_possible_score is not None for obj in self.objectives]
        ):
            best_possible_score = sum(
                [obj.best_possible_score for obj in self.objectives]
            )
        else:
            best_possible_score = None

        current_best_score = self.objective_scores_sum()
        current_best_sequence = self.sequence
        all_variants = self.mutation_space.all_variants(self.sequence)
        space_size = int(self.mutation_space.space_size)
        self.logger(mutation__total=space_size)
        for variant in self.logger.iter_bar(mutation=all_variants):
            self.sequence = variant
            if self.all_constraints_pass():
                score = self.objective_scores_sum()
                if score > current_best_score:
                    current_best_score = score
                    current_best_sequence = self.sequence
                    if (best_possible_score is not None) and (
                        current_best_score >= best_possible_score
                    ):
                        self.logger(mutation__index=space_size)
                        break
        self.sequence = current_best_sequence

    def optimize_by_random_mutations(self):
        """
        """
        if not self.all_constraints_pass():
            summary = self.constraints_text_summary()
            raise ValueError(
                summary + "Optimization can only be done when all"
                " constraints are verified"
            )
        score = self.objective_scores_sum()

        if all(
            [obj.best_possible_score is not None for obj in self.objectives]
        ):
            best_possible_score = sum(
                [
                    obj.best_possible_score * obj.boost
                    for obj in self.objectives
                ]
            )
        else:
            best_possible_score = None
        iters = self.max_random_iters
        stagnating_iterations = 0
        for iteration in self.logger.iter_bar(mutation=range(iters)):
            if (best_possible_score is not None) and (
                score >= best_possible_score
            ):
                self.logger(mutation__index=iters)
                break
            if (self.optimization_stagnation_tolerance is not None) and (
                stagnating_iterations > self.optimization_stagnation_tolerance
            ):
                break

            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=self.mutations_per_iteration,
                sequence=self.sequence,
            )
            if self.all_constraints_pass():
                new_score = self.objective_scores_sum()
                if new_score > score:
                    score = new_score
                    stagnating_iterations = 0
                else:
                    self.sequence = previous_sequence
            else:
                self.sequence = previous_sequence
            stagnating_iterations += 1

    def optimize_objective(self, objective):
        """Optimize the total objective score, focusing on a single objective.

        This method will attempt to increase the global objective score by
        focusing on a single objective. First the locations of under-optimal
        subsequences for this objective are identified, then these locations
        are optimized one after the other, left to right.

        For each location, a local problem is created and the optimization uses
        either a custom optimization algorithm, an exhaustive search, or a
        random search, to optimize the local problem
        """
        # EVALUATE OBJECTIVE. RETURN IF THERE IS NOTHING TO BE DONE.
        evaluation = objective.evaluate(self)
        locations = evaluation.locations
        if (objective.best_possible_score is not None) and (
            evaluation.score == objective.best_possible_score
        ):
            return

        # FOR EACH LOCATION, CREATE AND OPTIMIZE A LOCAL PROBLEM.

        for location in self.logger.iter_bar(
            location=locations, bar_message=lambda l: str(l)
        ):
            # Localize the mutation space by freezing any nucleotide outside of
            # it
            mutation_space = self.mutation_space.localized(location)
            if mutation_space.space_size == 0:
                continue

            # Update the location so it matches the span of the mutation_space
            # the resulting location will be equal or smaller to the original
            # location.
            location = Location(*mutation_space.choices_span)
            localized_constraints = [
                cst.localized(location, problem=self)
                for cst in self.constraints
            ]
            localized_constraints = [
                cst for cst in localized_constraints if cst is not None
            ]
            localized_objectives = [
                obj.localized(location, problem=self)
                for obj in self.objectives
                if obj.boost != 0
            ]
            localized_objectives = [
                obj for obj in localized_objectives if obj is not None
            ]
            local_problem = self.__class__(
                sequence=self.sequence,
                constraints=localized_constraints,
                mutation_space=mutation_space,
                objectives=localized_objectives,
            )
            self.logger.store(
                problem=self, local_problem=local_problem, location=location
            )
            local_problem.randomization_threshold = (
                self.randomization_threshold
            )
            local_problem.max_random_iters = self.max_random_iters
            local_problem.optimization_stagnation_tolerance = (
                self.optimization_stagnation_tolerance
            )
            local_problem.mutations_per_iteration = (
                self.mutations_per_iteration
            )

            # OPTIMIZE THE LOCAL PROBLEM

            if hasattr(objective, "optimization_heuristic"):
                # Some specifications implement their own optimization method.
                objective.optimization_heuristic(local_problem)
            else:
                # Run an exhaustive or random search depending on the size
                # of the mutation space.
                space_size = local_problem.mutation_space.space_size
                exhaustive_search = space_size < self.randomization_threshold
                if exhaustive_search:
                    local_problem.optimize_by_exhaustive_search()
                else:
                    local_problem.optimize_by_random_mutations()

            # UPDATE THE PROBLEM's SEQUENCE

            self.sequence = local_problem.sequence

    def optimize(self):
        """Maximize the total score by optimizing each objective in turn."""

        objectives = [
            obj
            for obj in self.objectives
            if not obj.optimize_passively and obj.boost != 0
        ]
        if len(objectives) == 0:
            return
        for objective in self.logger.iter_bar(
            objective=objectives, bar_message=lambda o: str(o)
        ):
            self.optimize_objective(objective=objective)
