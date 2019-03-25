"""Define the DnaOptimizationProblem class.

DnaOptimizationProblem is where the whole problem is defined: sequence,
constraints, objectives.
"""

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from .biotools.biotools import (sequence_to_biopython_record,
                                find_specification_in_feature,
                                sequences_differences_array,
                                sequences_differences_segments)
from .Specification import Specification, SpecificationsSet
from .SpecEvaluation import (ProblemObjectivesEvaluations,
                             ProblemConstraintsEvaluations)
from .Location import Location
from .MutationSpace import MutationSpace
from .reports.optimization_reports import (write_optimization_report,
                                           write_no_solution_report,
                                           PDF_REPORTS_AVAILABLE)
from proglog import default_bar_logger

DEFAULT_SPECIFICATIONS_DICT = {} # completed at library initialization


class NoSolutionError(Exception):
    """Exception returned when a DnaOptimizationProblem aborts.
    This means that the constraints are found to be unsatisfiable.
    """

    def __init__(self, message, problem, constraint=None, location=None):
        """Initialize."""
        Exception.__init__(self, message)
        self.message = message
        self.problem = problem
        self.constraint = constraint
        self.location = location

    def __str__(self):
        return self.message

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

    logger
      Either None for no logger, 'bar' for a tqdm progress bar logger, or
      any ProgLog progress bar logger.

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

    randomization_threshold = 10000
    max_random_iters = 1000
    n_mutations = 2
    local_extensions = (0, 5)

    def __init__(self, sequence, constraints=None, objectives=None,
                 logger='bar', mutation_space=None):
        """Initialize"""
        if isinstance(sequence, SeqRecord):
            self.record = sequence
            self.sequence = str(sequence.seq).upper()
        else:
            self.record = None
            self.sequence = sequence.upper()
        self.constraints = [] if constraints is None else list(constraints)
        self.objectives = [] if objectives is None else list(objectives)
        self.logger = default_bar_logger(
            logger, bars=('objective', 'constraint', 'location'),
            ignored_bars=('mutation',), min_time_interval=0.2)
        self.mutation_space = mutation_space
        self.initialize()

    def initialize(self):
        """Variables initialization before solving."""

        #Uncompress SpecificationSets into
        for specs in (self.constraints, self.objectives):
            specsets = [
                spec for spec in specs
                if isinstance(spec, SpecificationsSet)
            ]
            specs_in_sets = [
                spec
                for specset in specsets
                for spec in specset.specifications.values()
            ]
            for specset in specsets:
                specs.remove(specset)
            specs.extend(specs_in_sets)


        self.constraints = [
            constraint.initialize_on_problem(self, role="constraint")
            for constraint in self.constraints
        ]
        self.objectives = [
            objective.initialize_on_problem(self, role="objective")
            for objective in  self.objectives
        ]

        self.sequence_before = self.sequence
        self._constraints_before = None
        self._objectives_before = None
        if self.mutation_space is None:
            self.mutation_space = MutationSpace.from_optimization_problem(self)
            self.sequence = self.mutation_space.constrain_sequence(
                self.sequence)

    @property
    def constraints_before(self):
        if self._constraints_before is None:
            sequence = self.sequence
            self.sequence = self.sequence_before
            self._constraints_before = self.constraints_evaluations()
            self.sequence = sequence
        return self._constraints_before

    @property
    def objectives_before(self):
        if self._objectives_before is None:
            sequence = self.sequence
            self.sequence = self.sequence_before
            self._objectives_before = self.objectives_evaluations()
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

    def objective_scores_sum(self):
        return self.objectives_evaluations().scores_sum()

    def objectives_text_summary(self):
        return self.objectives_evaluations().to_text()

    def resolve_constraints_by_exhaustive_search(self):
        """Solve all constraints by exploring the whole search space.

        This method iterates over ``self.iter_mutations_space()`` (space of
        all sequences that could be reached through successive mutations) and
        stops when it finds a sequence which meets all the constraints of the
        canvas.
        """
        sequence_before = self.sequence
        all_variants = self.mutation_space.all_variants(self.sequence)
        space_size = int(self.mutation_space.space_size)
        self.logger(mutation__total=space_size)
        for variant in self.logger.iter_bar(mutation=all_variants):
            self.sequence = variant
            if self.all_constraints_pass():
                self.logger(mutation__index=space_size)
                return
        self.sequence = sequence_before
        raise NoSolutionError(
            "Exhaustive search failed to satisfy all constraints.",
            problem=self
        )

    def resolve_constraints_by_random_mutations(self):
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

        iters = range(3*self.max_random_iters)
        for i in self.logger.iter_bar(mutation=iters):

            if all(e.passes for e in evaluations):
                self.logger(mutation__index=iters)
                return
            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=self.n_mutations, sequence=self.sequence)

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
        raise NoSolutionError(
            "Random search did not find a solution in the given number of "
            "attempts. Try to increase the number of attempts with:\n\n"
            "problem.max_random_iters = 5000 # or even 10000, 20000, etc.\n\n"
            "If the problem persists, you may be in presence of a complex or "
            "unsolvable problem.",
            problem=self
        )

    def resolve_constraints_locally(self):
        """Orient the local search towards a stochastic or exhaustive search.
        """
        if self.mutation_space.space_size < self.randomization_threshold:
            self.resolve_constraints_by_exhaustive_search()
        else:
            self.resolve_constraints_by_random_mutations()

    def resolve_constraint(self, constraint):
        """Resolve a constraint through successive localizations."""
        evaluation = constraint.evaluate(self)
        if evaluation.passes:
            return

        locations = sorted(evaluation.locations)
        iterator = self.logger.iter_bar(location=locations,
                                        bar_message=lambda loc: str(loc))
        for i, location in enumerate(iterator):
            for extension in self.local_extensions:
                new_location = location.extended(extension)
                mutation_space = self.mutation_space.localized(new_location)
                if mutation_space.space_size == 0:
                    if extension == self.local_extensions[-1]:
                        error = NoSolutionError(
                            location=new_location,
                            problem=self,
                            message='Constraint breach in region that cannot '
                                    'be mutated.'
                        )
                        error.location = new_location
                        error.constraint = constraint
                        error.message = "While solving %s in %s:\n\n%s" % (
                            constraint, new_location, str(error)
                        )
                        self.logger(location__index=len(locations),
                                    location__message='Cold exit')
                        raise error
                    else:
                        continue
                new_location = Location(*mutation_space.choices_span)

                # This blocks solves the problem of overlapping breaches,
                # which can make the local optimization impossible.
                if (i < (len(locations) - 1)) and (
                   locations[i + 1].overlap_region(new_location)):
                    this_local_constraint = constraint.localized(
                        new_location, with_righthand=False, problem=self)
                else:
                    this_local_constraint = constraint.localized(
                        new_location, problem=self)

                if this_local_constraint.evaluate(self).passes:
                    continue
                localized_constraints = [
                    _constraint.localized(new_location, problem=self)
                    for _constraint in self.constraints
                    if _constraint != constraint
                    if not _constraint.enforced_by_nucleotide_restrictions
                ]
                localized_constraints = [cst for cst in localized_constraints
                                         if cst is not None]
                passing_localized_constraints = [
                    _constraint
                    for _constraint in localized_constraints
                    if _constraint.evaluate(self).passes
                ]
                local_problem = self.__class__(
                    sequence=self.sequence,
                    constraints=([this_local_constraint] +
                                 passing_localized_constraints),
                    mutation_space=mutation_space
                )
                self.logger.store(problem=self,
                                  local_problem=local_problem,
                                  location=location)
                local_problem.randomization_threshold = \
                    self.randomization_threshold
                local_problem.max_random_iters = self.max_random_iters
                local_problem.n_mutations = self.n_mutations
                try:
                    if hasattr(constraint, 'resolution_heuristic'):
                        constraint.resolution_heuristic(local_problem)
                    else:
                        local_problem.resolve_constraints_locally()
                    self.sequence = local_problem.sequence
                    break
                except NoSolutionError as error:
                    if extension == self.local_extensions[-1]:
                        error.location = new_location
                        error.constraint = constraint
                        error.message = "While solving %s in %s:\n\n%s" % (
                            constraint, new_location, str(error)
                        )
                        self.logger(location__index=len(locations),
                                    location__message='Cold exit')
                        raise error
                    else:
                        continue

    def resolve_constraints(self, final_check=True):
        """Solve a particular constraint using local, targeted searches.

        Parameters
        ----------

        constraint
          The ``Specification`` object for which the sequence should be solved

        """
        constraints = [
            c for c in self.constraints
            if not c.enforced_by_nucleotide_restrictions
        ]
        if len(constraints) == 0:
            return
        constraints = sorted(constraints, key=lambda c: -c.priority)
        for constraint in self.logger.iter_bar(constraint=constraints,
                                               bar_message=lambda c: str(c)):
            try:
                self.resolve_constraint(constraint=constraint)
            except NoSolutionError as error:
                self.logger(constraint__index=len(constraints))
                raise error
        if final_check:
            for cst in self.constraints:
                if not cst.evaluate(self).passes:
                    raise NoSolutionError(
                        'The solving of all constraints failed to solve'
                        ' all constraints. This is an unintended behavior,'
                        ' likely due to a complex problem. Try running the'
                        ' solver on the same sequence again, or report the'
                        ' error to the maintainers:\n\n' +
                        self.constraints_text_summary(failed_only=True),
                        problem=self)



    # SPECIFICATIONS

    def optimize_by_exhaustive_search(self):
        """
        """
        if not self.all_constraints_pass():
            summary = self.constraints_text_summary(failed_only=True)
            raise NoSolutionError(
                summary +
                "Optimization can only be done when all constraints are "
                "verified.",
                self
            )

        if all([obj.best_possible_score is not None
                for obj in self.objectives]):
            best_possible_score = sum([obj.best_possible_score
                                       for obj in self.objectives])
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
                    if ((best_possible_score is not None) and
                            (current_best_score >= best_possible_score)):
                        self.logger(mutation__index=space_size)
                        break
            self.sequence = self.sequence_before
        self.sequence = current_best_sequence
        assert self.all_constraints_pass()

    def optimize_by_random_mutations(self):
        """
        """

        if not self.all_constraints_pass():
            summary = self.constraints_text_summary()
            raise ValueError(summary + "Optimization can only be done when all"
                             " constraints are verified")
        score = self.objective_scores_sum()

        if all([obj.best_possible_score is not None
                for obj in self.objectives]):
            best_possible_score = sum([obj.best_possible_score * obj.boost
                                       for obj in self.objectives])
        else:
            best_possible_score = None
        iters = self.max_random_iters
        for iteration in self.logger.iter_bar(mutation=range(iters)):
            if ((best_possible_score is not None) and
                    (score >= best_possible_score)):
                self.logger(mutation__index=iters)
                break

            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                n_mutations=self.n_mutations, sequence=self.sequence)
            if self.all_constraints_pass():
                new_score = self.objective_scores_sum()
                if new_score > score:
                    score = new_score
                else:
                    self.sequence = previous_sequence
            else:
                self.sequence = previous_sequence
        #  assert self.all_constraints_pass()

    def optimize_objective(self, objective):

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

        for location in self.logger.iter_bar(location=locations,
                                             bar_message=lambda l: str(l)):
            mutation_space = self.mutation_space.localized(location)
            if mutation_space.space_size == 0:
                continue
            location = Location(*mutation_space.choices_span)
            localized_constraints = [
                _constraint.localized(location, problem=self)
                for _constraint in self.constraints
            ]
            localized_constraints = [cst for cst in localized_constraints
                                     if cst is not None]
            localized_objectives = [
                _objective.localized(location, problem=self)
                for _objective in self.objectives
            ]
            localized_objectives = [obj for obj in localized_objectives
                                    if obj is not None]
            local_problem = DnaOptimizationProblem(
                sequence=self.sequence,
                constraints=localized_constraints,
                mutation_space=mutation_space,
                objectives=localized_objectives
            )
            self.logger.store(problem=self,
                              local_problem=local_problem,
                              location=location)
            local_problem.randomization_threshold = self.randomization_threshold
            local_problem.max_random_iters = self.max_random_iters
            local_problem.n_mutations = self.n_mutations
            if hasattr(objective, 'optimization_heuristic'):
                objective.optimization_heuristic(local_problem)
            else:
                space_size = local_problem.mutation_space.space_size
                exhaustive_search = space_size < self.randomization_threshold
                if exhaustive_search:
                    local_problem.optimize_by_exhaustive_search()
                else:
                    local_problem.optimize_by_random_mutations()
            self.sequence = local_problem.sequence

    def optimize(self):
        """Maximize the objective via local, targeted mutations."""

        objectives = [
            obj for obj in self.objectives
            if not obj.optimize_passively
        ]
        if len(objectives) == 0:
            return
        for objective in self.logger.iter_bar(objective=objectives,
                                              bar_message=lambda o: str(o)):
            self.optimize_objective(objective=objective)

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
            role, spec = Specification.from_biopython_feature(
                feature, specifications_dict)
            parameters[role + "s"].append(spec)

        return DnaOptimizationProblem(**parameters)

    def to_record(self, filepath=None, features_type="misc_feature",
                  with_original_features=True,
                  with_original_spec_features=False,
                  with_constraints=True,
                  with_objectives=True,
                  with_sequence_edits=False,
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
                if with_original_spec_features or
                not find_specification_in_feature(f)
            ]
        if with_sequence_edits:
            record.features += self.sequence_edits_as_features()


        if filepath is not None:
            SeqIO.write(record, filepath, "genbank")
        else:
            return record

    def sequence_edits_as_array(self):
        return sequences_differences_array(self.sequence, self.sequence_before)


    def sequence_edits_as_features(self, feature_type="misc_feature"):

        segments = sequences_differences_segments(self.sequence,
                                                  self.sequence_before)
        return [
            Location(start, end).to_biopython_feature(
                label="was " + self.sequence_before[start: end],
                is_edit="true")
            for start, end in segments
        ]
    
    def optimize_with_report(self, target, project_name="My project"):
        """Resolve constraints, optimize objectives, write a multi-file report.

        The report's content may vary depending on the optimization's success 

        Parameters
        ----------
        target
        Either a path to a folder that will containt the report, or a path to
        a zip archive, or "@memory" to return raw data of a zip archive
        containing the report.

        project_name
        Project name to write on PDF reports

        **solver_parameters
        Extra keyword arguments passed to ``problem.resolve_constraints()``

        Returns
        -------
        (success, message, zip_data)
        Triplet where success is True/False, message is a one-line string
        summary indication whether some clash was found, or some solution, or
        maybe no solution was found because the random searches were too short
        """
        self.logger(message="Solving constraints")
        try:
            self.resolve_constraints()
        except NoSolutionError as error:
            self.logger(message="No solution found: making report")
            data = write_no_solution_report(target, self, error)
            start, end, s = error.location.to_tuple()
            message = ("No solution found in zone [%d, %d]: %s." %
                    (start, end, str(error)))
            return False, message, data
        self.logger(message="Now optimizing the sequence")
        self.optimize()
        self.logger(message="Success! Generating report.")
        data = write_optimization_report(
            target, self, project_name=project_name)
        return True, "Optimization successful.", data