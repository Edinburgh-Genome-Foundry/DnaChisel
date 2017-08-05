"""Definie MutationSpace"""

import itertools
import numpy as np
from .biotools import windows_overlap



class Mutation:

    def __init__(self, segment, subsequence):
        self.segment = segment
        self.start, self.end = segment
        self.subsequence = subsequence

    def __eq__(self, other):
        return (self.start, self.subsequence) == (other.start,
                                                  other.subsequence)

    def __hash__(self):
        return hash((self.start, self.subsequence))

    def __repr__(self):
        return "Mut(%d-%d %s)" % (self.start, self.end, self.subsequence)

    def __str__(self):
        return "Mut(%d-%d %s)" % (self.start, self.end, self.subsequence)

class MutationChoice:

    def __init__(self, segment, mutations):
        self.segment = segment
        self.start, self.end = segment
        self.mutations = mutations
        self.possible_subsequences = set(m.subsequence for m in mutations)

    def random_mutation(self, sequence):
        subsequence = sequence[self.start: self.end]
        mutations = [m for m in self.mutations if m.subsequence != subsequence]
        return mutations[np.random.randint(len(mutations))]

    def percolate_with(self, others):
        """Return a mutation restriction

        Example:
        --------

        overlap_sets = [
            ((0,3),set(['GTA', 'GCT', 'GTT'])),
            ((3,4), set("ATGC")),
            ((4,7), set(['ATG', 'ACC', 'CTG']))
        ]
        new_restriction = (2, 5), set(['ATT', 'ATA'])
        returns (0, 7): {'GTATACC', 'GTATATG'}

        """
        others = sorted(others, key=lambda o: o.start)
        others_start = others[0].start
        final_segment = others_start, others[-1].end
        final_mutations = set()
        for candidate in self.mutations:
            slots = []
            for other in others:
                istart, iend = windows_overlap(other.segment, self.segment)
                slot = []
                for mutation in other.mutations:
                    seq = mutation.subsequence
                    subseq = seq[istart - other.start: iend - other.start]
                    subcandidate = candidate.subsequence[
                        istart - self.start: iend - self.start]
                    if subseq == subcandidate:
                        slot.append(seq)
                slots.append(slot)
            for subseqs in itertools.product(*slots):
                seq = "".join(subseqs)
                matching_seq = seq[self.start - others_start:
                                   self.end - others_start]
                if matching_seq == candidate.subsequence:
                    final_mutations.add(Mutation(segment=final_segment,
                                                 subsequence=seq))

        return MutationChoice(segment=final_segment,
                              mutations=final_mutations)

    @staticmethod
    def from_tuple(segment_mutations):
        segment, mutations = segment_mutations
        if isinstance(segment, int):
            segment = segment, segment + 1
        mutations = [
            Mutation(segment=segment, subsequence=subsequence)
            for subsequence in mutations
        ]
        return MutationChoice(segment=segment, mutations=mutations)

    def __repr__(self):
        subsequences = "-".join([m.subsequence for m in self.mutations])
        return "MutChoice(%d-%d %s)" % (self.start, self.end, subsequences)

    def __str__(self):
        subsequences = "-".join([m.subsequence for m in self.mutations])
        return "MutChoice(%d-%d %s)" % (self.start, self.end, subsequences)


class MutationSpace:

    def __init__(self, choices_index):
        """

        choices_index = [MutationChoice(0-2), MutationChoice(0-2),
                         MutationChoice(3-5), MutationChoice(3-5),
                         MutationChoice(3-5), ... ]
        """
        self.choices_index = choices_index
        self.choices_list = sorted(set(choices_index), key=lambda c: c.start)
        self.unsolvable_segments = [
            choice.segment
            for choice in self.choices_list
            if len(choice.mutations) == 0
        ]
        self.determined_segments = [
            (choice.segment, list(choice.mutations)[0].subsequence)
            for choice in self.choices_list
            if len(choice.mutations) == 1
        ]
        self.choices = [
            choice
            for choice in self.choices_list
            if len(choice.mutations) > 1
        ]

    @property
    def choices_span(self):
        return self.choices_list[0].start, self.choices_list[-1].end

    def constrain_sequence(self, sequence):
        new_sequence = np.array(list(sequence))
        for choice in self.choices_list:
            possibilities = choice.possible_subsequences
            if len(possibilities) == 0:
                continue
            elif len(possibilities) == 1:
                possibilities = list(list(possibilities)[0])
                new_sequence[choice.start:choice.end] = possibilities
            elif sequence[choice.start: choice.end] not in possibilities:
                new_sequence[choice.start:choice.end] = list(possibilities)[0]
        return "".join(new_sequence)


    def localized(self, location):
        """Return a new version with only mutations overlaping the location."""
        if hasattr(location, "start"):
            start, end = location.start, location.end
        else:
            start, end = location
        return MutationSpace(self.choices_index[start:end])

    @property
    def space_size(self):
        """Return the number of possible mutations"""
        if len(self.choices) == 0:
            return 0
        return np.prod([1.0] + [
            len(choice.mutations)
            for choice in self.choices
        ])

    def pick_random_mutations(self, n_mutations, sequence):
        """Draw N random mutations"""
        n_mutations = min(len(self.choices), n_mutations)
        if n_mutations == 1:
            choice = self.choices[np.random.randint(len(self.choices))]
            return [choice.random_mutation(sequence=sequence)]

        return [
            choice_.random_mutation(sequence=sequence)
            for choice_ in [
                self.choices[i]
                for i in np.random.choice(len(self.choices), n_mutations,
                                          replace=False)
            ]
        ]


    def apply_random_mutations(self, n_mutations, sequence):
        """Return a sequence with n random mutations applied."""
        new_sequence = bytearray(sequence.encode())
        for mut in self.pick_random_mutations(n_mutations, sequence):
            new_sequence[mut.start: mut.end] = mut.subsequence.encode()
        return new_sequence.decode()

    def iterate_mutations(self, sequence):
        new_sequence = bytearray(sequence.encode())
        choice_start, choice_end = self.choices_span
        encoded_segment = sequence[choice_start: choice_end].encode()
        mutations_slots = [
             [(choice_, m.subsequence.encode()) for m in choice_.mutations]
             for choice_ in self.choices
        ]
        for mutations in itertools.product(*mutations_slots):
            new_sequence[choice_start: choice_end] = encoded_segment
            for (choice, subsequence) in mutations:
                new_sequence[choice.start: choice.end] = subsequence
            yield new_sequence.decode()


    @staticmethod
    def from_optimization_problem(optimization_problem):

        sequence = optimization_problem.sequence
        choice_index = [
            MutationChoice.from_tuple(((i, i+1), set("ATGC")))
            for i in range(len(sequence))
        ]

        mutation_choices = sorted([
            choice
            if isinstance(choice, MutationChoice)
            else MutationChoice.from_tuple(choice)
            for cst in optimization_problem.constraints
            for choice in cst.restrict_nucleotides(sequence)
        ], key=lambda choice: (choice.end - choice.start, choice.start))


        for choice in mutation_choices:
            underlying_choices = choice_index[choice.start: choice.end]
            new_choice = choice.percolate_with(set(underlying_choices))
            for i in range(new_choice.start, new_choice.end):
                choice_index[i] = new_choice
        return MutationSpace(choice_index)

    def unsolvable_nucleotides(self):
        return [
            loc for (loc, _set) in self.restrictions.items()
            if len(_set) == 0
        ]

    def adapt_to_sequence(self, sequence):

        new_sequence = list(sequence)
        for (start, end), _set in list(allowed_mutations.items()):
            if sequence[start:end] not in _set:
                new_sequence[start:end] = list(_set)[0]
            if len(_set) == 1:
                allowed_mutations.pop((start, end))

        new_sequence = "".join(new_sequence)
        return new_sequence
