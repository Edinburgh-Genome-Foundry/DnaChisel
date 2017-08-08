"""Definie MutationSpace"""

import itertools
import numpy as np
from .biotools import windows_overlap

class MutationChoice:
    __slots__ = ['segment', 'start', 'end', 'variants']

    def __init__(self, segment, variants):
        if isinstance(segment, int):
            segment = (segment, segment + 1)
        self.segment = segment
        self.start, self.end = segment
        self.variants = variants
        # self.possible_subsequences = set(m.subsequence for m in mutations)

    def random_mutation(self, sequence):
        subsequence = sequence[self.start: self.end]
        variants = [v for v in self.variants if v != subsequence]
        return variants[np.random.randint(len(variants))]

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
        final_variants = set()
        for candidate in self.variants:
            slots = []
            for other in others:
                istart, iend = windows_overlap(other.segment, self.segment)
                slot = []
                for variant in other.variants:
                    subseq = variant[istart - other.start: iend - other.start]
                    subcandidate = candidate[istart - self.start:
                                             iend - self.start]
                    if subseq == subcandidate:
                        slot.append(variant)
                slots.append(slot)
            for subseqs in itertools.product(*slots):
                seq = "".join(subseqs)
                matching_seq = seq[self.start - others_start:
                                   self.end - others_start]
                if matching_seq == candidate:
                    final_variants.add(seq)
        return MutationChoice(segment=final_segment,
                              variants=final_variants)

    def __repr__(self):
        subsequences = "-".join(self.variants)
        return "MutChoice(%d-%d %s)" % (self.start, self.end, subsequences)

    def __str__(self):
        subsequences = "-".join(self.variants)
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
            if len(choice.variants) == 0
        ]
        self.determined_segments = [
            (choice.segment, list(choice.variants)[0])
            for choice in self.choices_list
            if len(choice.variants) == 1
        ]
        self.choices = [
            choice
            for choice in self.choices_list
            if len(choice.variants) > 1
        ]

    @property
    def choices_span(self):
        if self.choices == []:
            return None
        return self.choices_list[0].start, self.choices_list[-1].end

    def constrain_sequence(self, sequence):
        new_sequence = np.array(list(sequence))
        for choice in self.choices_list:
            variants = choice.variants
            if len(choice.variants) == 0:
                continue
            elif len(variants) == 1:
                variants = list(list(variants)[0])
                new_sequence[choice.start:choice.end] = variants
            elif sequence[choice.start: choice.end] not in variants:
                new_sequence[choice.start:choice.end] = list(variants)[0]
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
            len(choice.variants)
            for choice in self.choices
        ])

    def pick_random_mutations(self, n_mutations, sequence):
        """Draw N random mutations"""
        n_mutations = min(len(self.choices), n_mutations)
        if n_mutations == 1:
            choice = self.choices[np.random.randint(len(self.choices))]
            return [
                (choice.segment,
                 choice.random_mutation(sequence=sequence))
            ]

        return [
            (choice_.segment, choice_.random_mutation(sequence=sequence))
            for choice_ in [
                self.choices[i]
                for i in np.random.choice(len(self.choices), n_mutations,
                                          replace=False)
            ]
        ]


    def apply_random_mutations(self, n_mutations, sequence):
        """Return a sequence with n random mutations applied."""
        new_sequence = bytearray(sequence.encode())
        for segment, seq in self.pick_random_mutations(n_mutations, sequence):
            start, end = segment
            new_sequence[start: end] = seq.encode()
        return new_sequence.decode()

    def all_variants(self, sequence):
        new_sequence = bytearray(sequence.encode())
        choice_start, choice_end = self.choices_span
        encoded_segment = sequence[choice_start: choice_end].encode()
        variants_slots = [
             [(choice_.segment, v.encode()) for v in choice_.variants]
             for choice_ in self.choices
        ]
        for variants in itertools.product(*variants_slots):
            new_sequence[choice_start: choice_end] = encoded_segment
            for ((start, end), variant) in variants:
                new_sequence[start: end] = variant
            yield new_sequence.decode()


    @staticmethod
    def from_optimization_problem(optimization_problem):

        sequence = optimization_problem.sequence
        choice_index = [
            MutationChoice((i, i+1), variants=set("ATGC"))
            for i in range(len(sequence))
        ]

        mutation_choices = sorted([
            choice
            if isinstance(choice, MutationChoice)
            else MutationChoice(segment=choice[0], variants=set(choice[1]))
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
