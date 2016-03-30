"""Define DNACanvas.

DNACanvas is where the whole problem is defined: sequence,
constraints, objectives.
"""

from copy import deepcopy, copy
import ctypes
import itertools as itt

import numpy as np

from biotools import translate, reverse_translate, gc_percent, read_fasta, reverse_complement
import biotables
import constraints as cst


class DNACanvas:

    def __init__(self, sequence, constraints=None, objectives=None):

        self.sequence = sequence
        self.original_sequence = sequence
        self.constraints = [] if constraints is None else constraints
        self.objectives = [] if objectives is None else objectives

        self.compute_possible_mutations()

    # MUTATIONS

    def compute_possible_mutations(self):
        self.possible_mutations = {}
        unibase_mutable = np.ones(len(self.sequence))
        for constraint in self.constraints:
            if isinstance(constraint, cst.DoNotModifyConstraint):
                start, end = constraint.window
                unibase_mutable[start:end] = 0
        for constraint in self.constraints:
            if isinstance(constraint, cst.EnforceTranslationConstraint):
                start, end = constraint.window
                for i, aa in enumerate(constraint.translation):
                    if constraint.strand == 1:
                        cstart, cstop = start + 3 * i, start + 3 * (i + 1)
                        seq_codon = self.sequence[cstart:cstop]
                    else:
                        cstart, cstop = end - 3 * (i+1), end - 3 * i
                        seq_codon = reverse_complement(
                                        self.sequence[cstart:cstop])
                    possible_codons = biotables.CODONS_SEQUENCES[aa][:]
                    local_immutable_unibases = (
                        unibase_mutable[cstart:cstop] == 0
                    ).nonzero()[0]

                    def array_subsequence(seq, inds):
                        return np.array([seq[i] for i in inds])
                    if len(local_immutable_unibases):
                        reachable_possible_codons = [
                            codon
                            for codon in possible_codons
                            if all(
                                array_subsequence(
                                    seq_codon,
                                    local_immutable_unibases
                                ) ==
                                array_subsequence(
                                    codon,
                                    local_immutable_unibases
                                )
                            )
                        ]
                        if reachable_possible_codons == []:
                            raise ValueError(
                                "An EnforceTranslation constraint seems to"
                                " clash with a DoNotTouch constraint."
                            )
                        possible_codons = reachable_possible_codons
                    if constraint.strand == -1:
                        possible_codons = [
                            reverse_complement(possible_codon)
                            for possible_codon in possible_codons
                        ]
                    unibase_mutable[cstart:cstop] = 0

                    if seq_codon in possible_codons:
                        possible_codons.remove(seq_codon)

                    if possible_codons != []:
                        self.possible_mutations[
                            (cstart, cstop)] = possible_codons
        #print unibase_mutable
        for i in unibase_mutable.nonzero()[0]:
            self.possible_mutations[i] = ["A", "T", "G", "C"]

    def mutation_space_size(self):
        return np.prod([
            len(v) + 1.0
            for v in self.possible_mutations.values()
        ])

    def iter_mutations_space(self):
        return itt.product(*[
            [None] + [(k, seq) for seq in values]
            for k, values in self.possible_mutations.items()
        ])

    def mutate_sequence(self, mutations):
        sequence_buffer = ctypes.create_string_buffer(self.sequence)
        for mutation in mutations:
            if mutation is not None:
                ind, seq = mutation
                if isinstance(ind, int):
                    sequence_buffer[ind] = seq
                else:
                    start, end = ind
                    sequence_buffer[start:end] = seq
        self.sequence = sequence_buffer.value

     # CONSTRAINTS

    def all_constraints_evaluations(self):
        return [
            constraint.evaluate(self)
            for constraint in self.constraints
        ]

    def all_constraints_pass(self):
        return all([
            evaluation.passes
            for evaluation in self.all_constraints_evaluations()
        ])

    def print_constraints_summary(self, failed_only=False):
        evaluations = self.all_constraints_evaluations()
        failed_evaluations = [e for e in evaluations if not e.passes]
        if failed_only:
            evaluations = failed_evaluations
        if failed_evaluations == []:
            print("===> SUCCESS - all constraints evaluations pass")
        else:
            print ("===> FAILURE: %d constraints evaluations failed" %
                   len(failed_evaluations))
        for evaluation in evaluations:
            print("%s %s" % (evaluation.constraint, evaluation))

    def solve_all_constraints_by_exhaustive_search(self, verbose=False):
        for mutations in self.iter_mutations_space():
            self.mutate_sequence(mutations)
            if verbose:
                self.print_constraints_summary()
            if self.all_constraints_pass():
                return
            else:
                self.sequence = self.original_sequence
        raise ValueError(
            "Exhaustive search failed to satisfy all constraints.")

    def solve_all_constraints_by_random_mutations(self, max_iter=1000,
                                                  n_mutations=3,
                                                  verbose=False):
        mutations_locs = self.possible_mutations.keys()
        evaluations = self.all_constraints_evaluations()
        score = sum([
            e.score
            for e in evaluations
            if not e.passes
        ])
        for iteration in range(max_iter):
            if score == 0:
                return

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
                self.print_constraints_summary()
            previous_sequence = self.sequence
            self.mutate_sequence(mutations)

            evaluations = self.all_constraints_evaluations()
            new_score = sum([
                e.score
                for e in evaluations
                if not e.passes
            ])
            # print "now scores with muts", map(str,evaluations), new_score, score
            if new_score > score:
                score = new_score
            else:
                self.sequence = previous_sequence
        raise ValueError(
            "Random search hit max_iterations without finding a solution.")

    def solve_constraint_by_localization(self, constraint,
                                         randomization_threshold=10000,
                                         max_random_iters=1000, verbose=False):
        """Solve the constraint using local, targeted mutations.


        """

        evaluation = constraint.evaluate(self)
        if evaluation.passes:
            return
        if evaluation.windows is not None:

            for window in evaluation.windows:
                if verbose:
                    print(window)
                do_not_modify_window = [
                    max(0, window[0] - 5),
                    min(window[1] + 5, len(self.sequence))
                ]
                localized_constraints = [
                    _constraint.localized(window)
                    for _constraint in self.constraints
                ]
                passing_localized_constraints = [
                    _constraint
                    for _constraint in localized_constraints
                    if _constraint.evaluate(self).passes
                ]
                localized_canvas = DNACanvas(
                    sequence=self.sequence,
                    constraints=[
                        cst.DoNotModifyConstraint([0, do_not_modify_window[0]]),
                        cst.DoNotModifyConstraint([do_not_modify_window[1],
                                                   len(self.sequence)]),
                    ] + [
                        constraint.localized(window)
                    ] + passing_localized_constraints
                )
                #print constraint, localized_canvas.mutation_space_size()
                #print localized_canvas.possible_mutations
                if (localized_canvas.mutation_space_size() <
                        randomization_threshold):
                    localized_canvas.solve_all_constraints_by_exhaustive_search(
                        verbose=verbose)
                    self.sequence = localized_canvas.sequence
                else:
                    localized_canvas.solve_all_constraints_by_random_mutations(
                        max_iter=max_random_iters, n_mutations=1,
                        verbose=verbose)
                    self.sequence = localized_canvas.sequence

    def solve_all_constraints_one_by_one(self, max_loops=3,
                                         randomization_threshold=10000,
                                         max_random_iters=1000, verbose=False):

        for iteration in range(max_loops):
            evaluations = self.all_constraints_evaluations()
            failed_constraints = [
                evaluation.constraint
                for evaluation in evaluations
                if not evaluation.passes
            ]
            if failed_constraints == []:
                return
            for constraint in failed_constraints:
                self.solve_constraint_by_localization(
                    constraint, randomization_threshold,
                    max_random_iters, verbose=verbose
                )

        raise ValueError(
            "Could not solve all constraints before reaching max_loops"
        )

    # OBJECTIVES
    def all_objectives_evaluations(self):
        return [
            objective.evaluate(self)
            for objective in self.objectives
        ]

    def all_objectives_score_sum(self):
        return sum([
            objective.boost * objective.evaluate(self).score
            for objective in self.objectives
        ])

    def print_objectives_summary(self, failed_only=False):
        print("===> TOTAL OBJECTIVES SCORE: %.02f" %
              self.all_objectives_score_sum())
        for evaluation in self.all_objectives_evaluations():
            print "%s: %s" % (evaluation.objective, evaluation)

    def maximize_objectives_by_exhaustive_search(self, verbose=False):
        """
        """
        if not self.all_constraints_pass():
            raise ValueError("Optimization can only be done when all"
                             " constraints are verified")
        current_score = self.all_objectives_score_sum()
        current_best_sequence = self.sequence
        for mutations in self.iter_mutations_space():
            self.mutate_sequence(mutations)
            if self.all_constraints_pass():
                score = self.all_objectives_score()
                if score > current_score:
                    current_score = score
                    current_best_sequence = self.sequence
            self.sequence = self.original_sequence
        self.sequence = current_best_sequence

    def maximize_objectives_by_random_mutations(self, max_iter=1000,
                                                n_mutations=3,
                                                verbose=False):
        """
        """
        if not self.all_constraints_pass():
            raise ValueError("Optimization can only be done when all"
                             " constraints are verified")
        mutations_locs = self.possible_mutations.keys()
        score = self.all_objectives_score_sum()
        for iteration in range(max_iter):
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
                self.print_constraints_summary()
            previous_sequence = self.sequence
            self.mutate_sequence(mutations)
            if self.all_constraints_pass():
                new_score = self.all_objectives_score_sum()
                if new_score > score:
                    score = new_score
                else:
                    self.sequence = previous_sequence
            else:
                self.sequence = previous_sequence

    def maximize_objective_by_localization(self, objective, windows=None,
        randomization_threshold=10000, max_random_iters=1000, verbose=False):
        """Maximize the objective via local, targeted mutations.
        """

        if windows is None:
            windows = objective.evaluate(self).windows
            if windows is None:
                raise ValueError(
                    "max_objective_by_localization requires either that"
                    " windows be provided or that the objective evaluation"
                    " returns windows."
                )
        for window in windows:
            if verbose:
                print(window)
            do_not_modify_window = [
                max(0, window[0] - 5),
                min(window[1] + 5, len(self.sequence))
            ]
            localized_canvas = DNACanvas(
                sequence=self.sequence,
                constraints=[
                    _constraint.localized(window)
                    for _constraint in self.constraints
                ] + [
                    cst.DoNotModifyConstraint([0, do_not_modify_window[0]]),
                    cst.DoNotModifyConstraint([do_not_modify_window[1],
                                               len(self.sequence)]),
                ],
                objectives = [
                    _objective.localized(window)
                    for _objective in self.objectives
                ]
            )

            if (localized_canvas.mutation_space_size() <
                    randomization_threshold):
                localized_canvas.maximize_objectives_by_exhaustive_search(
                    verbose=verbose)
            else:
                localized_canvas.maximize_objectives_by_random_mutations(
                    max_iter=max_random_iters, n_mutations=1,
                    verbose=verbose)
            self.sequence = localized_canvas.sequence


    def maximize_all_objectives_one_by_one(self, n_loops=1,
                                           randomization_threshold=10000,
                                           max_random_iters=1000,
                                           verbose=False):

        for iteration in range(n_loops):
            for objective in self.objectives:
                self.maximize_objective_by_localization(
                    objective,
                    randomization_threshold=randomization_threshold,
                    max_random_iters=max_random_iters,
                    verbose=verbose
                )
