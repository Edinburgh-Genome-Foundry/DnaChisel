"""Define the DnaOptimizationProblem class.

DnaOptimizationProblem is where the whole problem is defined: sequence,
constraints, objectives.
"""

import itertools as itt

import numpy as np

from .biotools.biotools import (sequence_to_biopython_record,
                                find_specification_in_feature,
                                sequences_differences_segments)

from .specifications import (Specification, AvoidChanges, EnforcePattern,
                             ProblemObjectivesEvaluations,
                             ProblemConstraintsEvaluations,
                             DEFAULT_SPECIFICATIONS_DICT)

from .Location import Location
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from .Mutation import MutationSpace


from tqdm import tqdm


class NoSolutionError(Exception):
    """Exception returned when a DnaOptimizationProblem aborts.
    This means that the constraints are found to be unsatisfiable.
    """

    def __init__(self, message, problem, location=None):
        """Initialize."""
        Exception.__init__(self, message)
        self.problem = problem
        self.location = location


class NoSolutionFoundError(Exception):
    """Exception returned when a DnaOptimizationProblem solving fails to find
    a solution.
    This means that the search was not long enough or that there was
    no solution, because in some complex way, some constraints
    are incompatible.
    """
    pass


class DnaOptimizationProblem:
    """Problem specifications: sequence, constraints, optimization objectives.

    The original constraints, objectives, and original sequence of the problem
    are stored in the DNA Canvas. This class also has methods to display
    reports on the constraints and objectives, as well as solving the
    constraints and objectives.

    Examples
    --------

    >>> from dnachisel import *
    >>> canvas = DnaOptimizationProblem(
    >>>     sequence = "ATGCGTGTGTGC...",
    >>>     constraints = [constraint1, constraint2, ...],
    >>>     objectives = [objective1, objective2, ...]
    >>> )
    >>> canvas.resolve_constraints_one_by_one()
    >>> canvas.maximize_all_objectives_one_by_one()
    >>> print(canvas.constraints_text_summary())
    >>> print(canvas.objectives_summary())


    Parameters
    ----------

    sequence
      A string of ATGC characters (they must be upper case!), e.g. "ATTGTGTA"

    constraints
      A list of objects of type ``Specification``.

    objectives
      A list of objects of type ``Specification`` specifying what must be
      optimized in the problem. Note that each objective has a float ``boost``
      parameter. The larger the boost, the more the objective is taken into
      account during the optimization.

    Attributes
    ----------

    sequence
      The sequence

    constraints
      The list of constraints

    objectives
      The list of objectives

    possible_mutations
      A dictionnary indicating the possible mutations. Note that this is only
      computed the first time that canvas.possible_mutations is invoked.

    Notes
    -----

    The dictionnary ``self.possible_mutations`` is of the form
    ``{location1 : list1, location2: list2...}``
    where ``location`` is either a single index (e.g. 10) indicating the
    position of a nucleotide to be muted, or a couple ``(start, end)``
    indicating a whole segment whose sub-sequence should be replaced.
    The ``list`` s are lists of possible sequences to replace each location,
    e.g. for the mutation of a whole codon ``(3,6): ["ATT", "ACT", "AGT"]``.
    """

    # location_extension = 2

    def __init__(self, sequence, constraints=None, objectives=None,
                 progress_logger=None, mutation_space=None):
        """Initialize"""
        if isinstance(sequence, SeqRecord):
            self.record = sequence
            self.sequence = str(sequence.seq)
        else:
            self.record = None
            self.sequence = sequence
        self.constraints = [] if constraints is None else constraints
        self.objectives = [] if objectives is None else objectives
        if progress_logger is None:
            def progress_logger(*a, **k):
                pass
        self.progress_logger = progress_logger
        self.mutation_space = mutation_space
        self.initialize()

    def initialize(self):
        """Variables initialization before solving."""

        self.progress_logger(message="Initializing...")

        self.constraints = [
            constraint.initialize_on_problem(self, role="constraint")
            for constraint in self.constraints
        ]
        self.objectives = [
            objective.initialize_on_problem(self, role="objective")
            for objective in self.objectives
        ]

        self.sequence_before = self.sequence
        self._constraints_before = None
        self._objectives_before = None
        if self.mutation_space is None:
            self.mutation_space = MutationSpace.from_optimization_problem(self)
            self.sequence = self.mutation_space.constrain_sequence(
                self.sequence)

    def constraints_before(self):
        if self._constraints_before is None:
            sequence = self.sequence
            self.sequence = self.sequence_before
            self._constraints_before = self.constraints_evaluations()
            self.sequence = sequence
        return self._constraints_before

    def objectives_before(self):
        if self._objectives_before is None:
            sequence = self.sequence
            self.sequence = self.sequence_before
            self._obejctives_before = self.objectives_evaluations()
            self.sequence = sequence
        return self._objectives_before

    def constraints_evaluations(self):
        """Return a list of the evaluations of each constraint of the canvas.
        """
        return ProblemConstraintsEvaluations.from_problem(self)

    def all_constraints_pass(self):
        """Return True iff the current problem sequence passes all constraints.
        """
        return self.constraints_evaluations().all_evaluations_pass()

    def constraints_text_summary(self, failed_only=False):
        evals = self.constraints_evaluations()
        if failed_only:
            evals = evals.filter("failing")
        return evals.to_text()

    def objectives_evaluations(self):
        """Return a list of the evaluation of each objective of the canvas"""
        return ProblemObjectivesEvaluations.from_problem(self)

    def all_objectives_scores_sum(self):
        return self.objectives_evaluations().scores_sum()

    def objectives_text_summary(self):
        return self.objectives_evaluations().to_text()

    def extract_subsequence(self, location):
        """Return the subsequence (a string) at the given location).

         The ``location`` can be either an index (an integer) indicating the
         position of a single nucleotide, or a list/couple ``(start, end)``
         indicating a whole sub-segment.
        """
        if np.isscalar(location):
            return self.sequence[location]
        else:
            start, end = location
            return self.sequence[start:end]

    def resolve_constraints_by_exhaustive_search(self,
                                                 progress_bar=False, raise_exception_on_failure=True):
        """Solve all constraints by exploring the whole search space.

        This method iterates over ``self.iter_mutations_space()`` (space of
        all sequences that could be reached through successive mutations) and
        stops when it finds a sequence which meets all the constraints of the
        canvas.
        """
        sequence_before = self.sequence
        all_mutations = self.mutation_space.iterate_mutations(self.sequence)
        if progress_bar:
            all_mutations = tqdm(all_mutations, desc="Mutation", leave=False,
                                 total=self.mutation_space.space_size)
        for mutated_sequence in all_mutations:
            self.sequence = mutated_sequence
            if self.all_constraints_pass():
                return
        self.sequence = sequence_before

        if raise_exception_on_failure:
            raise NoSolutionError(
                "Exhaustive search failed to satisfy all constraints.",
                problem=self)

    def resolve_constraints_by_random_mutations(
            self, max_iter=1000, n_mutations=1, progress_bar=False,
            raise_exception_on_failure=False):
        """Solve all constraints by successive sets of random mutations.

        This method modifies the canvas sequence by applying a number
        ``n_mutations`` of random mutations. The constraints are then evaluated
        on the new sequence. If all constraints pass, the new sequence becomes
        the canvas's new sequence.
        If not all constraints pass, the sum of all scores from failing
        constraints is considered. If this score is superior to the score of
        the previous sequence, the new sequence becomes the canvas's new
        sequence.

        This operation is repeated `max_iter` times at most, after which
        a ``NoSolutionError`` is thrown.


        """

        evaluations = self.constraints_evaluations()
        score = sum([
            e.score
            for e in evaluations
            if not e.passes
        ])
        range_iter = range(max_iter)
        if progress_bar:
            range_iter = tqdm(range_iter, desc="Random Mutation", leave=False)
        for iteration in range_iter:
            if all(e.passes for e in evaluations):
                return
            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=n_mutations, sequence=self.sequence)

            evaluations = self.constraints_evaluations()
            new_score = sum([
                e.score
                for e in evaluations
                if not e.passes
            ])

            if new_score > score:
                score = new_score
            else:
                self.sequence = previous_sequence
        if raise_exception_on_failure:
            raise NoSolutionError(
                "Random search hit max_iterations without finding a solution.",
                problem=self
            )

    def resolve_constraints(self, constraints='all',
                            randomization_threshold=10000,
                            max_random_iters=1000,
                            progress_bars=0,
                            evaluation=None,
                            n_mutations=1,
                            final_check=True):
        """Solve a particular constraint using local, targeted searches.

        Parameters
        ----------

        constraint
          The ``Specification`` object for which the sequence should be solved

        randomization_threshold
          Local problems with a search space size under this threshold will be
          solved using deterministic, exhaustive search of the search space
          (see ``resolve_constraints_by_exhaustive_search``)
          When the space size is above this threshold, local searches will use
          a randomized search algorithm
          (see ``resolve_constraints_by_random_mutations``).

        max_random_iters
          Maximal number of iterations when performing a randomized search
          (see ``resolve_constraints_by_random_mutations``).

        """
        if constraints == 'all':
            constraints = self.constraints

        if isinstance(constraints, list):
            iter_constraints = constraints
            if progress_bars > 0:
                iter_constraints = tqdm(constraints, desc="Constraint",
                                        leave=False)
            self.progress_logger(n_constraints=len(self.constraints),
                                 constraint_ind=0)
            for i, constraint in enumerate(iter_constraints):
                self.resolve_constraints(
                    constraints=constraint,
                    randomization_threshold=randomization_threshold,
                    max_random_iters=max_random_iters,
                    progress_bars=progress_bars - 1, n_mutations=n_mutations
                )
                self.progress_logger(constraint_ind=i + 1)
            if final_check:
                for cst in constraints:
                    if not cst.evaluate(self).passes:
                        raise NoSolutionFoundError(
                            'The solving of all constraints failed to solve'
                            ' all constraints, constraint %s is'
                            ' still failing. This is an unintended behavior,'
                            ' likely due to a complex problem. Try running the'
                            ' solver on the same sequence again, or report the'
                            ' error to the maintainers' % cst)

            return

        constraint = constraints
        evaluation = constraint.evaluate(self)
        if evaluation.passes:
            return

        locations = evaluation.locations
        self.progress_logger(n_locations=len(locations), location_ind=0)
        if progress_bars > 0:
            locations = tqdm(locations, desc="Window", leave=False)

        for i, location in enumerate(locations):
            mutation_space = self.mutation_space.localized(location)
            location = Location(*mutation_space.choices_span)
            localized_constraints = [
                _constraint.localized(location)
                for _constraint in self.constraints
            ]
            passing_localized_constraints = [
                _constraint
                for _constraint in localized_constraints
                if _constraint.evaluate(self).passes
            ]
            local_problem = DnaOptimizationProblem(
                sequence=self.sequence,
                constraints=([constraint.localized(location)] +
                             passing_localized_constraints),
                mutation_space=mutation_space
            )

            space_size = local_problem.mutation_space.space_size
            exhaustive_search = space_size < randomization_threshold
            try:
                if exhaustive_search:
                    local_problem.resolve_constraints_by_exhaustive_search(
                        progress_bar=progress_bars > 1)
                    self.sequence = local_problem.sequence
                else:
                    local_problem.resolve_constraints_by_random_mutations(
                        max_iter=max_random_iters, n_mutations=n_mutations,
                        progress_bar=progress_bars > 1)
                    self.sequence = local_problem.sequence
            except NoSolutionError as error:
                error.location = location
                raise error
            self.progress_logger(location_ind=i + 1)

    # SPECIFICATIONS

    def optimize_by_exhaustive_search(self, progress_bar=False):
        """
        """
        if not self.all_constraints_pass():
            summary = self.constraints_text_summary(failed_only=True)
            raise NoSolutionError(
                summary +
                "Optimization can only be done when all constraints are"
                "verified."
            )

        if all([obj.best_possible_score is not None
                for obj in self.objectives]):
            best_possible_score = sum([obj.best_possible_score
                                       for obj in self.objectives])
        else:
            best_possible_score = None

        current_best_score = self.all_objectives_scores_sum()
        current_best_sequence = self.sequence
        all_mutations = self.mutation_space.iterate_mutations(self.sequence)
        if progress_bar:
            all_mutations = tqdm(all_mutations, desc="Mutation", leave=False,
                                 total=self.mutation_space.space_size)
        for mutated_sequence in all_mutations:
            self.sequence = mutated_sequence
            if self.all_constraints_pass():
                score = self.all_objectives_scores_sum()
                if score > current_best_score:
                    current_best_score = score
                    current_best_sequence = self.sequence
                    if ((best_possible_score is not None) and
                            (current_best_score >= best_possible_score)):
                        break
            self.sequence = self.sequence_before
        self.sequence = current_best_sequence
        assert self.all_constraints_pass()

    def optimize_by_random_mutations(self, max_iter=1000, n_mutations=1,
                                     progress_bar=False):
        """
        """

        if not self.all_constraints_pass():
            summary = self.constraints_text_summary()
            raise ValueError(summary + "Optimization can only be done when all"
                             " constraints are verified")
        score = self.all_objectives_scores_sum()

        if all([obj.best_possible_score is not None
                for obj in self.objectives]):
            best_possible_score = sum([obj.best_possible_score * obj.boost
                                       for obj in self.objectives])
        else:
            best_possible_score = None

        range_iters = range(max_iter)
        if progress_bar:
            range_iters = tqdm(range_iters, desc="Random mutation",
                               leave=False)

        for iteration in range_iters:
            if ((best_possible_score is not None) and
                    (score >= best_possible_score)):
                break

            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=n_mutations, sequence=self.sequence)
            if self.all_constraints_pass():
                new_score = self.all_objectives_scores_sum()
                if new_score > score:
                    score = new_score
                else:
                    self.sequence = previous_sequence
            else:
                self.sequence = previous_sequence
        #  assert self.all_constraints_pass()

    def optimize(self, objectives='all', n_mutations=1,
                 randomization_threshold=10000,
                 max_random_iters=1000,
                 optimize_independently=False,
                 progress_bars=False):
        """Maximize the objective via local, targeted mutations."""

        if objectives == 'all':
            objectives = [
                obj for obj in self.objectives
                if not obj.is_passive_objective
            ]

        if isinstance(objectives, list):

            iter_objectives = objectives
            if progress_bars > 0:
                iter_objectives = tqdm(objectives, desc="Objective",
                                       leave=False)
            self.progress_logger(n_objectives=len(objectives),
                                 objective_ind=0)
            for i, objective in enumerate(iter_objectives):
                self.optimize(
                    objectives=objective,
                    randomization_threshold=randomization_threshold,
                    max_random_iters=max_random_iters,
                    progress_bars=(progress_bars - 1),
                    n_mutations=n_mutations,
                    optimize_independently=optimize_independently
                )
                self.progress_logger(objective_ind=i + 1)

            return

        objective = objectives
        evaluation = objective.evaluate(self)
        locations = evaluation.locations
        if ((objective.best_possible_score is not None) and
                (evaluation.score == objective.best_possible_score)):
            return
        if locations is None:
            raise ValueError(
                ("With %s:" % objective) +
                "max_objective_by_localization requires either that"
                " locations be provided or that the objective evaluation"
                " returns locations."
            )

        if progress_bars > 0:
            locations = tqdm(locations, desc="Window", leave=False)

        for location in locations:
            mutation_space = self.mutation_space.localized(location)
            location = Location(*mutation_space.choices_span)
            localized_constraints = [
                _constraint.localized(location)
                for _constraint in self.constraints
            ]
            local_problem = DnaOptimizationProblem(
                sequence=self.sequence,
                constraints=localized_constraints,
                mutation_space=mutation_space,
                objectives=[
                    _objective.localized(location)
                    for _objective in self.objectives
                ]
            )
            space_size = local_problem.mutation_space.space_size
            exhaustive_search = space_size < randomization_threshold
            if exhaustive_search:
                local_problem.optimize_by_exhaustive_search(
                    progress_bar=progress_bars > 1)
            else:
                local_problem.optimize_by_random_mutations(
                    max_iter=max_random_iters, n_mutations=n_mutations,
                    progress_bar=progress_bars > 1)
            self.sequence = local_problem.sequence

    def include_pattern_by_successive_tries(self, pattern, location=None,
                                            paste_pattern_in=True,
                                            aim_at_location_centre=True):
        if location is None:
            location = Location(0, len(self.sequence))
        constraint = EnforcePattern(pattern=pattern, location=location)
        self.constraints.append(constraint)

        indices = list(range(location.start, location.end - pattern.size))
        if aim_at_location_centre:
            center = 0.5 * (location.start + location.end - pattern.size)
            indices = sorted(indices, key=lambda i: abs(i - center))

        for i in indices:
            sequence_before = self.sequence
            if paste_pattern_in:
                location = [i, i + pattern.size]
                self.mutate_sequence([(location, pattern.pattern)])
            try:
                self.resolve_constraints_one_by_one()
                return i
            except NoSolutionError:
                self.sequence = sequence_before
        self.constraints.pop()  # remove the pattern constraint
        raise NoSolutionError("Failed to insert the pattern")

    @staticmethod
    def from_record(record, specifications_dict="default"):
        if isinstance(record, str):
            if record.lower().endswith((".fa", ".fasta")):
                record = SeqIO.read(record, 'fasta')
            elif record.lower().endswith((".gb", ".gbk")):
                record = SeqIO.read(record, 'genbank')
            else:
                raise ValueError("Record is either a Biopython record or a "
                                 "file path ending in .gb, .gbk, .fa, .fasta.")
        if specifications_dict == "default":
            specifications_dict = DEFAULT_SPECIFICATIONS_DICT
        parameters = dict(
            sequence=record,
            constraints=[],
            objectives=[]
        )
        for feature in record.features:
            if feature.type != "misc_feature":
                continue
            if find_specification_in_feature(feature) is None:
                continue
            role, objective = Specification.from_biopython_feature(
                feature, specifications_dict)
            parameters[role + "s"].append(objective)

        return DnaOptimizationProblem(**parameters)

    def to_record(self, filepath=None, features_type="misc_feature",
                  with_original_features=True,
                  with_original_objective_features=False,
                  with_constraints=True,
                  with_objectives=True,
                  colors_dict=None):
        record = sequence_to_biopython_record(self.sequence)

        record.features = []
        if with_constraints:
            record.features += [
                cst.to_biopython_feature(role="constraint",
                                         feature_type=features_type,
                                         colors_dict=colors_dict)
                for cst in self.constraints
                if cst.__dict__.get('location', False)
            ]
        if with_objectives:
            record.features += [
                obj.to_biopython_feature(role="objective",
                                         feature_type=features_type,
                                         colors_dict=colors_dict)
                for obj in self.objectives
            ]
        if with_original_features and (self.record is not None):
            record.features += [
                f for f in self.record.features
                if with_original_objective_features or
                not find_specification_in_feature(f)
            ]

        if filepath is not None:
            SeqIO.write(record, filepath, "genbank")
        else:
            return record

    def sequence_edits_as_features(self, feature_type="misc_feature"):

        segments = sequences_differences_segments(self.sequence,
                                                  self.sequence_before)
        return [
            Location(start, end).to_biopython_feature(
                label="was " + self.sequence_before[start: end],
                is_edit="true")
            for start, end in segments
        ]
