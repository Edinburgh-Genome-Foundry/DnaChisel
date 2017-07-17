"""Define the DnaOptimizationProblem class.

DnaOptimizationProblem is where the whole problem is defined: sequence,
constraints, objectives.
"""

import itertools as itt

import numpy as np

from .biotools.biotools import (reverse_complement,
                                sequence_to_biopython_record,
                                find_objective_in_feature,
                                sequences_differences_segments)

from .objectives import (Objective, DoNotModify, EnforcePattern,
                         EnforceTranslation,
                         ProblemObjectivesEvaluations,
                         ProblemConstraintsEvaluations,
                         DEFAULT_OBJECTIVES_DICT)

from .Location import Location
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


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
    >>> canvas.solve_all_constraints_one_by_one()
    >>> canvas.maximize_all_objectives_one_by_one()
    >>> print(canvas.constraints_text_summary())
    >>> print(canvas.objectives_summary())


    Parameters
    ----------

    sequence
      A string of ATGC characters (they must be upper case!), e.g. "ATTGTGTA"

    constraints
      A list of objects of type ``Objective``.

    objectives
      A list of objects of type ``Objective`` specifying what must be optimized
      in the problem. Note that each objective has a float ``boost`` parameter.
      The larger the boost, the more the objective is taken into account during
      the optimization.

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

    def __init__(self, sequence, constraints=None, objectives=None,
                 progress_logger=None):
        if isinstance(sequence, SeqRecord):
            self.record = sequence
            self.sequence = str(sequence.seq)
        else:
            self.record = None
            self.sequence = sequence
        self.possible_mutations_dict = None
        self.constraints = [] if constraints is None else constraints
        self.objectives = [] if objectives is None else objectives
        if progress_logger is None:
            def progress_logger(*a, **k):
                pass
        self.progress_logger = progress_logger
        self.initialize()

    def initialize(self):
        """Variables initialization before solving."""

        self.progress_logger(message="Initializing...")

        self.constraints = [
            constraint.initialize_problem(self, role="constraint")
            for constraint in self.constraints
        ]
        self.objectives = [
            objective.initialize_problem(self, role="objective")
            for objective in self.objectives
        ]

        self.sequence_before = self.sequence
        self.constraints_before = self.constraints_evaluations()
        self.objectives_before = self.objectives_evaluations()


    def constraints_evaluations(self):
        """Return a list of the evaluations of each constraint of the canvas.
        """
        return ProblemConstraintsEvaluations.from_problem(self)

    def all_constraints_pass(self):
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

    # MUTATIONS

    @property
    def possible_mutations(self):
        if self.possible_mutations_dict is None:
            self.compute_possible_mutations()
        return self.possible_mutations_dict

    def compute_possible_mutations(self):
        """Compute all possible mutations that can be applied to the sequence.

        The result of the computations is stored in ``self.possible_mutations``
        (see ``DnaOptimizationProblem`` documentation).

        The possible mutations are constrained by the ``constraints`` of the
        DnaOptimizationProblem with respect to the following rules:

        - ``DoNotModify``  constraints disable mutations for the nucleotides of
          the concerned segments.
        - ``EnforceTranlation`` constraints ensure that on the concerned
          segments only codons that translate to the imposed amino-acid will
          be considered, so a triplet of nucleotides that should code for
          Asparagin will see its choices down to ["AAT", "AAC"], instead of
          the 64 possible combinations of free triplets.

        """
        self.possible_mutations_dict = {}
        unibase_mutable = np.ones(len(self.sequence))
        for constraint in self.constraints:
            if isinstance(constraint, DoNotModify):
                if constraint.location is not None:
                    unibase_mutable[constraint.location.start:
                                    constraint.location.end] = 0
                else:
                    unibase_mutable[constraint.indices] = 0
        for constraint in self.constraints:
            if (isinstance(constraint, EnforceTranslation) and
                (constraint.codons_sequences is not None)):
                start, end = constraint.location.start, constraint.location.end
                for i, aa in enumerate(constraint.translation):
                    if constraint.location.strand == 1:
                        cstart, cstop = start + 3 * i, start + 3 * (i + 1)
                        seq_codon = self.sequence[cstart:cstop]
                    else:
                        cstart, cstop = end - 3 * (i + 1), end - 3 * i
                        seq_codon = reverse_complement(
                            self.sequence[cstart:cstop])
                    possible_codons = constraint.codons_sequences[aa][:]
                    local_immutable_unibases = (
                        unibase_mutable[cstart:cstop] == 0
                    ).nonzero()[0]

                    def array_subsequence(seq, inds):
                        return np.array([ord(seq[i]) for i in inds])

                    def codon_is_compatible(codon):
                        a1 = array_subsequence(seq_codon,
                                               local_immutable_unibases)
                        a2 = array_subsequence(codon,
                                               local_immutable_unibases)
                        return (a1 == a2).min() == 1

                    if len(local_immutable_unibases):
                        reachable_possible_codons = [
                            codon
                            for codon in possible_codons
                            if codon_is_compatible(codon)
                        ]
                        if reachable_possible_codons == []:

                            raise NoSolutionError(
                                "An EnforceTranslation constraint seems to"
                                " clash with a DoNotTouch constraint: %s" %
                                constraint,
                                problem=self,
                                location=Location(cstart, cstop,
                                                  constraint.location.strand)
                            )
                        possible_codons = reachable_possible_codons
                    if constraint.location.strand == -1:
                        possible_codons = [
                            reverse_complement(possible_codon)
                            for possible_codon in possible_codons
                        ]
                    unibase_mutable[cstart:cstop] = 0

                    if possible_codons not in [[], [seq_codon]]:
                        self.possible_mutations_dict[(cstart, cstop)] = \
                            possible_codons
        # print unibase_mutable
        for i in unibase_mutable.nonzero()[0]:
            self.possible_mutations_dict[i] = ["A", "T", "G", "C"]

    def mutation_space_size(self):
        """Return the total number of possible sequence variants.

        The result is a float.
        """
        return np.prod([
            1.0 * len(v)
            for v in self.possible_mutations.values()
        ])

    def iter_mutations_space(self):
        return itt.product(*[
            [(k, seq) for seq in values]
            for k, values in self.possible_mutations.items()
        ])

    def get_random_mutations(self, n_mutations=1):
        """Pick a random set of possible mutations.

        Returns a list ``[(location1, new_sequence1), ...]`` where location is
        either an index or a couple ``(start, end)`` and ``new_sequence`` is
        a DNA string like ``ATGC``, indicating that the canvas' sequence
        should be modified through mutations of the form
        ``self.sequence[location1] = new_sequence1``.
        """
        locs = list(self.possible_mutations.keys())
        if n_mutations == 1:
            indices = [np.random.randint(0, len(locs), 1)[0]]
        else:
            indices = np.random.choice(range(0, len(locs)), n_mutations,
                                       replace=False)
        mutations = []
        for index in indices:
            location = locs[index]
            subsequence = self.extract_subsequence(location)
            choices = self.possible_mutations[location]
            if subsequence in choices:
                choices.remove(subsequence)
            if choices == []:
                mutations.append(None)
            else:
                choice = np.random.choice(choices)
                mutations.append((location, choice))
        return mutations

    def mutate_sequence(self, mutations):
        """Modify the canvas's sequence (inplace) through mutations.

        ``mutations`` must be a list ``[(location1, new_sequence1), ...]``
        where location is either an index or a couple ``(start, end)`` and
        ``new_sequence`` is a DNA string like ``ATGC``, indicating that the
        canvas' sequence should be modified through mutations of the form
        ``self.sequence[location1] = new_sequence1``.
        """
        sequence_buffer = np.fromstring(self.sequence, dtype=np.uint8)
        for mutation in mutations:
            if mutation is not None:
                ind, seq = mutation
                if np.isscalar(ind):
                    sequence_buffer[ind] = ord(str(seq))
                else:
                    start, end = ind
                    sequence_buffer[start:end] = np.fromstring(str(seq),
                                                               dtype=np.uint8)
        self.sequence = sequence_buffer.tostring().decode("utf8")

    # CONSTRAINTS

    def solve_all_constraints_by_exhaustive_search(self, verbose=False,
        progress_bar=False, raise_exception_on_failure=True):
        """Solve all constraints by exploring the whole search space.

        This method iterates over ``self.iter_mutations_space()`` (space of
        all sequences that could be reached through successive mutations) and
        stops when it finds a sequence which meets all the constraints of the
        canvas.
        """
        mutation_space = self.iter_mutations_space()
        if progress_bar:
            mutation_space = tqdm(mutation_space, desc="Mutation", leave=False)
        for mutations in mutation_space:
            self.mutate_sequence(mutations)
            if verbose:
                print(self.constraints_text_summary())
            if self.all_constraints_pass():
                return
            else:
                self.sequence = self.sequence_before

        if raise_exception_on_failure:
            raise NoSolutionError(
                "Exhaustive search failed to satisfy all constraints.",
                problem=self)

    def solve_all_constraints_by_random_mutations(
        self, max_iter=1000, n_mutations=1, verbose=False, progress_bar=False,
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
            if score == 0:
                return
            mutations = self.get_random_mutations(n_mutations)
            if verbose:
                print(self.constraints_text_summary())
            previous_sequence = self.sequence
            self.mutate_sequence(mutations)

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
                problem=problem
            )

    def solve_constraint_by_localization(self, constraint,
                                         randomization_threshold=10000,
                                         max_random_iters=1000, verbose=False,
                                         progress_bars=0,
                                         evaluation=None,
                                         n_mutations=1,
                                         consider_other_constraints=True):
        """Solve a particular constraint using local, targeted searches.

        Parameters
        ----------

        constraint
          The ``Objective`` object for which the sequence should be solved

        randomization_threshold
          Local problems with a search space size under this threshold will be
          solved using deterministic, exhaustive search of the search space
          (see ``solve_all_constraints_by_exhaustive_search``)
          When the space size is above this threshold, local searches will use
          a randomized search algorithm
          (see ``solve_all_constraints_by_random_mutations``).

        max_random_iters
          Maximal number of iterations when performing a randomized search
          (see ``solve_all_constraints_by_random_mutations``).

        verbose
          If True, each step of each search will print in the console the
          evaluation of each constraint.

        """
        if evaluation is None:
            evaluation = constraint.evaluate(self)

        if evaluation.passes:
            return

        if evaluation.locations is not None:
            locations = evaluation.locations
            self.progress_logger(n_locations=len(locations))
            if progress_bars > 0:
                locations = tqdm(locations, desc="Window", leave=False)

            for i, location in enumerate(locations):
                if verbose:
                    print(location)
                do_not_modify_location = location.extended(
                    5, upper_limit=len(self.sequence))
                if consider_other_constraints:
                    localized_constraints = [
                        _constraint.localized(do_not_modify_location)
                        for _constraint in self.constraints
                    ]
                    passing_localized_constraints = [
                        _constraint
                        for _constraint in localized_constraints
                        if _constraint.evaluate(self).passes
                    ]
                else:
                    passing_localized_constraints = []
                localized_problem = DnaOptimizationProblem(
                    sequence=self.sequence,
                    constraints=[
                        DoNotModify(
                            location=Location(0, do_not_modify_location.start)
                        ),
                        DoNotModify(
                            location=Location(do_not_modify_location.end,
                                              len(self.sequence))
                        )
                    ] + [
                        constraint.localized(do_not_modify_location)
                    ] + passing_localized_constraints
                )
                try:
                    if (localized_problem.mutation_space_size() <
                            randomization_threshold):
                        localized_problem.solve_all_constraints_by_exhaustive_search(
                            verbose=verbose, progress_bar=progress_bars > 1)
                        self.sequence = localized_problem.sequence
                    else:
                        localized_problem.solve_all_constraints_by_random_mutations(
                            max_iter=max_random_iters, n_mutations=n_mutations,
                            verbose=verbose, progress_bar=progress_bars > 1)
                        self.sequence = localized_problem.sequence
                except NoSolutionError as error:
                    error.location = location
                    raise error
                self.progress_logger(location_ind=i+1)

    def solve_all_constraints_one_by_one(self, max_loops=2,
                                         randomization_threshold=10000,
                                         max_random_iters=1000, verbose=False,
                                         progress_bars=0,
                                         n_mutations=1,
                                         solve_independently=False):
        """Solve each of the canvas' constraints in turn, using local, targeted
        searches.

        Parameters
        ----------

        max_loops
          Number of times that the constraints will be considered one after the
          other. The function may stop sooner, as soon as all constraints pass.
          If after all these loops some constraints are still not passing, a
          ``NoSolutionError`` is raised.

        randomization_threshold
          Local problems with a search space size under this threshold will be
          solved using deterministic, exhaustive search of the search space
          (see ``solve_all_constraints_by_exhaustive_search``)
          When the space size is above this threshold, local searches will use
          a randomized search algorithm
          (see ``solve_all_constraints_by_random_mutations``).

        max_random_iters
          Maximal number of iterations when performing a randomized search
          (see ``solve_all_constraints_by_random_mutations``).

        verbose
          If True, each step of each search will print in the console the
          evaluation of each constraints.

        """

        range_loops = range(max_loops)
        if progress_bars > 0:
            range_loops = tqdm(range_loops, desc="Loop", leave=False)
        self.progress_logger(n_iterations=max_loops)
        for iteration in range_loops:

            evaluations = self.constraints_evaluations()
            failed_evaluations = [
                evaluation
                for evaluation in evaluations
                if not evaluation.passes
            ]
            if failed_evaluations == []:
                return
            self.progress_logger(n_evaluations=len(failed_evaluations))
            if progress_bars > 1:
                failed_evaluations = tqdm(failed_evaluations, leave=False,
                                          desc="Failing constraint")
            for i, evaluation in enumerate(failed_evaluations):
                self.solve_constraint_by_localization(
                    evaluation.objective, randomization_threshold,
                    max_random_iters, verbose=verbose,
                    progress_bars=progress_bars - 2,
                    evaluation=evaluation,
                    n_mutations=n_mutations,
                    consider_other_constraints=not solve_independently
                )
                self.progress_logger(evaluation_ind=i + 1)
            self.progress_logger(iteration_ind=iteration + 1)
        if not self.all_constraints_pass():
            raise NoSolutionFoundError(
                "One-by-one could not solve all constraints before max_loops."
            )

    # OBJECTIVES

    def maximize_objectives_by_exhaustive_search(self, verbose=False,
                                                 progress_bar=False):
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
        mutation_space = self.iter_mutations_space()
        if progress_bar:
            mutation_space = tqdm(mutation_space, desc="Mutation", leave=False)
        for mutations in mutation_space:
            self.mutate_sequence(mutations)
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

    def maximize_objectives_by_random_mutations(self, max_iter=1000,
                                                n_mutations=1,
                                                verbose=False,
                                                progress_bar=False):
        """
        """
        if not self.all_constraints_pass():
            summary = self.constraints_text_summary()
            raise ValueError(summary + "Optimization can only be done when all"
                             " constraints are verified")
        mutations_locs = list(self.possible_mutations.keys())
        score = self.all_objectives_scores_sum()
        range_iters = range(max_iter)

        if all([obj.best_possible_score is not None
                for obj in self.objectives]):
            best_possible_score = sum([obj.best_possible_score
                                       for obj in self.objectives])
        else:
            best_possible_score = None

        if progress_bar:
            range_iters = tqdm(range_iters, desc="Random mutation",
                               leave=False)
        for iteration in range_iters:
            random_mutations_inds = np.random.randint(
                0, len(mutations_locs), n_mutations)
            mutations = [
                (mutations_locs[ind],
                 np.random.choice(
                    self.possible_mutations[mutations_locs[ind]], 1
                )[0]
                )
                for ind in random_mutations_inds
            ]
            if verbose:
                print(self.constraints_text_summary())
            previous_sequence = self.sequence
            self.mutate_sequence(mutations)
            if self.all_constraints_pass():
                new_score = self.all_objectives_scores_sum()
                if new_score > score:
                    score = new_score
                else:
                    self.sequence = previous_sequence
            else:
                self.sequence = previous_sequence

            if ((best_possible_score is not None) and
                (score >= best_possible_score)):
                break
        assert self.all_constraints_pass()

    def maximize_objective_by_localization(self, objective, locations=None,
                                           randomization_threshold=10000,
                                           max_random_iters=1000,
                                           verbose=False,
                                           progress_bars=False,
                                           n_mutations=1,
                                           optimize_independently=False):
        """Maximize the objective via local, targeted mutations."""
        if locations is None:
            evaluation = objective.evaluate(self)
            locations = evaluation.locations
            if ((objective.best_possible_score is not None) and
                (evaluation.score == objective.best_possible_score)):
                return
            if locations is None:
                raise ValueError(("With %s:" % objective) +
                    "max_objective_by_localization requires either that"
                    " locations be provided or that the objective evaluation"
                    " returns locations."
                )

        if progress_bars > 0:
            locations = tqdm(locations, desc="Window", leave=False)
        for location in locations:
            if verbose:
                print(location)
            do_not_modify_location = location.extended(
                5, upper_limit=len(self.sequence))
            if optimize_independently:
                objectives = [objective.localized(do_not_modify_location)]
            else:
                objectives = [
                    _objective.localized(do_not_modify_location)
                    for _objective in self.objectives
                ]

            localized_problem = DnaOptimizationProblem(
                sequence=self.sequence,
                constraints=[
                    _constraint.localized(do_not_modify_location)
                    for _constraint in self.constraints
                ] + [
                        DoNotModify(
                            location=Location(0, do_not_modify_location.start)
                        ),
                        DoNotModify(
                            location=Location(do_not_modify_location.end,
                                              len(self.sequence))
                        )
                ],
                objectives=objectives
            )

            if (localized_problem.mutation_space_size() <
                    randomization_threshold):
                localized_problem.maximize_objectives_by_exhaustive_search(
                    verbose=verbose, progress_bar=progress_bars > 1)
            else:
                localized_problem.maximize_objectives_by_random_mutations(
                    max_iter=max_random_iters, n_mutations=n_mutations,
                    verbose=verbose, progress_bar=progress_bars > 1)
            self.sequence = localized_problem.sequence

    def maximize_all_objectives_one_by_one(self, n_loops=1,
                                           randomization_threshold=10000,
                                           max_random_iters=1000,
                                           verbose=False,
                                           progress_bars=0,
                                           optimize_independently=False,
                                           n_mutations=1):

        range_loops = range(n_loops)
        if progress_bars > 0:
            range_loops = tqdm(range_loops, desc="Loop", leave=False)

        for iteration in range_loops:
            objectives = self.objectives
            if progress_bars > 1:
                objectives = tqdm(objectives, desc="Objective", leave=False)
            for objective in objectives:
                self.maximize_objective_by_localization(
                    objective,
                    randomization_threshold=randomization_threshold,
                    max_random_iters=max_random_iters,
                    verbose=verbose,
                    progress_bars=(progress_bars - 2),
                    n_mutations=n_mutations,
                    optimize_independently=optimize_independently
                )

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
                self.solve_all_constraints_one_by_one()
                return i # success
            except NoSolutionError:
                self.sequence = sequence_before
        self.constraints.pop()  # remove the pattern constraint
        raise NoSolutionError("Failed to insert the pattern")

    @staticmethod
    def from_record(record, objectives_dict="default"):
        if isinstance(record, str):
            if record.lower().endswith((".fa", ".fasta")):
                record = SeqIO.read(record, 'fasta')
            elif record.lower().endswith((".gb", ".gbk")):
                record = SeqIO.read(record, 'genbank')
            else:
                raise ValueError("Record is either a Biopython record or a "
                                 "file path ending in .gb, .gbk, .fa, .fasta.")
        if objectives_dict == "default":
            objectives_dict = DEFAULT_OBJECTIVES_DICT
        parameters = dict(
            sequence=record,
            constraints=[],
            objectives=[]
        )
        for feature in record.features:
            if feature.type != "misc_feature":
                continue
            if find_objective_in_feature(feature) is None:
                continue
            role, objective = Objective.from_biopython_feature(feature,
                                                               objectives_dict)
            parameters[role+"s"].append(objective)

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
                not find_objective_in_feature(f)
            ]

        if filepath is not None:
            SeqIO.write(record, filepath, "genbank")
        else:
            return record

    def sequence_edits_as_features(self, feature_type="misc_feature"):

        segments = sequences_differences_segments(self.sequence,
                                                  self.sequence_before)
        return [
            Location(start, end).to_biopython_feature(label="edit",
                                                      is_edit="true")
            for start, end in segments
        ]
