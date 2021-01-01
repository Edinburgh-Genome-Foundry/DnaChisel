"""Define MutationSpace"""

import itertools
import numpy as np

from .MutationChoice import MutationChoice


class MutationSpace:
    """Class for mutation space (set of sequence segments and their variants).

    Parameters
    ----------

    choices_index
      A list L such that L[i] gives the MutationChoice governing the mutations
      allowed at position i (and possibly around i).


    Examples
    --------

    >>> # BEWARE: below, similar mutation choices are actually the SAME OBJECT
    >>> space = MutationSpace([
            MutationChoice((0, 2), {'AT', 'TG'}),
            MutationChoice((0, 2), {'AT', 'TG'}),
            MutationChoice((2, 5), {'TTC', 'TTA', 'TTT'}), # same
            MutationChoice((2, 5), {'TTC', 'TTA', 'TTT'}), #
            MutationChoice((2, 5), {'TTC', 'TTA', 'TTT'}),
        ])
    """

    def __init__(self, choices_index, left_padding=0):
        """

        choices_index = [MutationChoice(0-2), MutationChoice(0-2),
                         MutationChoice(3-5), MutationChoice(3-5),
                         MutationChoice(3-5), ... ]
        """
        self.choices_index = left_padding * [None] + choices_index
        self.choices_list = []
        self.unsolvable_segments = []
        self.determined_segments = []
        self.multichoices = []
        for c in choices_index:
            if c is None:
                continue
            if len(self.choices_list) == 0 or (c != self.choices_list[-1]):
                self.choices_list.append(c)
                nvariants = len(c.variants)
                if nvariants == 0:
                    self.unsolvable_segments.append(c.segment)
                elif nvariants == 1:
                    self.determined_segments.append((c.segment, list(c.variants)[0]))
                else:
                    self.multichoices.append(c)

    @property
    def choices_span(self):
        """Return (start, end), segment where multiple choices are possible."""
        if self.multichoices == []:
            return None
        return self.multichoices[0].start, self.multichoices[-1].end

    def constrain_sequence(self, sequence):
        """Return a version of the sequence compatible with the mutation space.

        All nucleotides of the sequence that are incompatible with the
        mutation space are replaced by nucleotides compatible with the space.
        """
        new_sequence = bytearray(sequence.encode())
        for choice in self.choices_list:
            variants = list(choice.variants)
            if len(choice.variants) == 0:
                raise ValueError(
                    "Cannot constrain a sequence when some "
                    "positions are unsolvable, in location "
                    "(%d-%d)" % (choice.start, choice.end)
                )
            elif len(variants) == 1:
                variant = variants[0]
                new_sequence[choice.start : choice.end] = variant.encode()
            elif sequence[choice.start : choice.end] not in variants:
                variant = sorted(variants)[np.random.randint(0, len(variants))]
                new_sequence[choice.start : choice.end] = variant.encode()
        return new_sequence.decode()

    def localized(self, location):
        """Return a new version with only mutations overlapping the location."""
        if hasattr(location, "start"):
            start, end = location.start, location.end
        else:
            start, end = location
        return MutationSpace(self.choices_index[start:end], left_padding=start)

    @property
    def space_size(self):
        """Return the number of possible mutations."""
        if len(self.multichoices) == 0:
            return 0
        choices = [len(choice.variants) for choice in self.multichoices]
        # np.prod(choices) can create overflows and warnings, so instead we use
        # the mechanism below with log/exp and a min.
        return np.exp(min(100, np.log(choices).sum()))

    def pick_random_mutations(self, n_mutations, sequence):
        """Draw N random mutations."""
        n_mutations = min(len(self.multichoices), n_mutations)
        if n_mutations == 1:
            index = np.random.randint(len(self.multichoices))
            choice = self.multichoices[index]
            return [(choice.segment, choice.random_variant(sequence=sequence))]

        return [
            (choice_.segment, choice_.random_variant(sequence=sequence))
            for choice_ in [
                self.multichoices[i]
                for i in np.random.choice(
                    len(self.multichoices), n_mutations, replace=False
                )
            ]
        ]

    def apply_random_mutations(self, n_mutations, sequence):
        """Return a sequence with n random mutations applied."""
        new_sequence = bytearray(sequence.encode())
        for segment, seq in self.pick_random_mutations(n_mutations, sequence):
            start, end = segment
            new_sequence[start:end] = seq.encode()
        return new_sequence.decode()

    def all_variants(self, sequence):
        """Iterate through all sequence variants in this mutation space."""
        new_sequence = bytearray(sequence.encode())
        choice_start, choice_end = self.choices_span
        encoded_segment = sequence[choice_start:choice_end].encode()

        def sort_variants_by_distance_to_current(choice):
            """This function iterates through the variants of a given choice
            using not the alphabetical (which would bias AC over GT) but rather
            a kind of 'least-change' order, which biases towards solutions
            close to the current sequence.

            Impact on overall algorithm speed is < 0.5%."""
            current = sequence[choice.segment[0] : choice.segment[1]]
            alphasort = {v: i for i, v in enumerate(sorted(choice.variants))}

            def sort_key(v):
                return (abs(alphasort[v] - alphasort[current]), v)

            return sorted(choice.variants, key=sort_key)

        variants_slots = [
            [
                (choice_.segment, v.encode())
                for v in sort_variants_by_distance_to_current(choice_)
            ]
            for choice_ in self.multichoices
        ]
        for variants in itertools.product(*variants_slots):
            new_sequence[choice_start:choice_end] = encoded_segment
            for ((start, end), variant) in variants:
                new_sequence[start:end] = variant
            yield new_sequence.decode()

    @staticmethod
    def from_optimization_problem(problem, new_constraints=None):
        """Create a mutation space from a DNA optimization problem.

        This can be used to initialize mutation spaces for new problems.
        """

        sequence = problem.sequence

        if new_constraints is None:
            variants = {"A": "ATGC", "T": "TACG", "G": "GCAT", "C": "CGTA"}
            choices_index = [
                MutationChoice((i, i + 1), variants=variants[c], is_any_nucleotide=True)
                for i, c in enumerate(sequence)
            ]
            constraints = problem.constraints
        else:
            choices_index = [c for c in problem.mutation_space.choices_index]
            constraints = new_constraints
        mutation_choices = sorted(
            [
                MutationChoice(segment=choice[0], variants=set(choice[1]))
                for cst in constraints
                for choice in cst.restrict_nucleotides(sequence)
            ],
            key=lambda choice: (choice.end - choice.start, choice.start),
        )
        for choice in mutation_choices:
            underlying_choices = choices_index[choice.start : choice.end]
            if underlying_choices == []:
                new_choice = choice
            elif all(c.is_any_nucleotide for c in underlying_choices):
                new_choice = choice
            else:
                new_choice = choice.merge_with(set(underlying_choices))
            for choice in new_choice.extract_varying_region():
                if choice.end > len(choices_index):
                    choices_index += (choice.end - len(choices_index)) * [None]
                for i in range(choice.start, choice.end):
                    choices_index[i] = choice
            # for i in range(new_choice.start, new_choice.end):
            #     choices_index[i] = new_choice
        return MutationSpace(choices_index)

    def string_representation(self, separator="|"):
        """Generates a string of the mutation space.

        Examples
        --------
        >>> print (mutation_space.string_representation())
        >>> A|G|ATG|CA|A|CA|T|GTC|AC|T|A|C|T|T|T
        >>> T|C|   |  |G|  |C|   |  |C|T|A|C|A|C
        >>> G|A|   |  | |  | |   |  |G|G|G|G|G|G
        """
        depth = max([len(c.variants) for c in self.choices_list])
        sequences = ["" for i in range(depth)]
        for choice in self.choices_list:
            variants = list(choice.variants)
            start, end = choice.segment
            length = end - start
            if len(variants) < depth:
                variants += (depth - len(variants)) * [length * " "]
            for i, variant in enumerate(variants):
                if len(sequences[i]):
                    sequences[i] += separator
                sequences[i] += variant
        return "\n".join(sequences)

    def plot(self, ax, color="red", y_offset=-2):
        for choice in self.choices_list:
            start, end = choice.segment
            N = len(choice.variants)
            delta = 0.2
            for x, sign in [(start, +1), (end, -1)]:
                _x = x - 0.5 + sign * delta
                _y = y_offset + 0.5
                ax.plot(
                    [_x, _x], [_y, _y - N], linewidth=0.5, color=color,
                )
            for y, variant in enumerate(choice.variants):
                for x, nucleotide in enumerate(variant):
                    ax.plot(
                        [start - 0.5 + delta, end - 0.5 - delta],
                        [y_offset - y - 0.5, y_offset - y - 0.5],
                        color=color,
                        linewidth=0.5,
                    )
                    ax.text(
                        x + start,
                        -y + y_offset,
                        nucleotide,
                        horizontalalignment="center",
                        verticalalignment="center",
                        color=color,
                    )
        max_variants = max([len(c.variants) for c in self.choices_list])
        ax.set_ylim(bottom=min(-max_variants - 2, ax.get_ylim()[0]))
