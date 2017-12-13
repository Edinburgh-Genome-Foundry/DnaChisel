import numpy as np

from .CodonSpecification import CodonSpecification
from ..SpecEvaluation import SpecEvaluation
from ..biotools import (CODON_USAGE_TABLES, CODONS_TRANSLATIONS,
                        group_nearby_indices, codons_frequencies_and_positions,
                        CODON_USAGE_BY_AA)
from ..Location import Location


class CodonOptimize(CodonSpecification):
    """Specification to codon-optimize a coding sequence for a particular species.

    Several codon-optimization policies exist. At the moment this Specification
    implements a method in which codons are replaced by the most frequent
    codon in the species.

    (as long as this doesn't break any Specification or lowers the global
    optimization objective)

    Supported speciess are ``E. coli``, ``S. cerevisiae``, ``H. Sapiens``,
    ``C. elegans``, ``D. melanogaster``, ``B. subtilis``.

    Parameters
    ----------

    species
      Name of the species to codon-optimize for. Supported speciess are
      ``E. coli``, ``S. cerevisiae``, ``H. Sapiens``, ``C. elegans``,
      ``D. melanogaster``, ``B. subtilis``.
      Note that the species can be omited if a ``codon_usage_table`` is
      provided instead

    mode
      Either 'best_codon' or 'harmonized'. For 'best_codon', the optimization
      will always replace a codon with the most-frequent triplet possible.
      For 'harmonized', the optimization will bring the relative frequencies of
      the different triplets as close as possible as the frequencies in the
      reference species.

    location
      Either a DnaChisel Location or a tuple of the form (start, end, strand)
      or just (start, end), with strand defaulting to +1, indicating the
      position of the gene to codon-optimize. If not provided, the whole
      sequence is considered as the gene. The location should have a length
      that is a multiple of 3. The location strand is either 1 if the gene is
      encoded on the (+) strand, or -1 for antisense.

    codon_usage_table
      A dict of the form ``{"TAC": 0.112, "CCT": 0.68}`` giving the RSCU table
      (relative usage of each codon). Only provide if no ``species`` name
      was provided.

    Examples
    --------

    >>> objective = CodonOptimizationSpecification(
    >>>     species = "E. coli",
    >>>     location = (150, 300), # coordinates of a gene
    >>>     strand = -1
    >>> )


    """

    best_possible_score = 0
    localization_group_spread = 3

    def __init__(self, species=None, location=None, mode='best_codon',
                 codon_usage_table=None, boost=1.0):
        self.mode = mode
        self.boost = boost
        if isinstance(location, tuple):
            location = Location.from_tuple(location, default_strand=+1)
        self.location = location
        self.species = species
        if species is not None:
            codon_usage_table = CODON_USAGE_TABLES[self.species]
        if codon_usage_table is None:
            raise ValueError("Provide either an species name or a codon "
                             "usage table")
        self.codon_usage_table = codon_usage_table

    def initialize_on_problem(self, problem, role):
        """Get location from sequence if no location provided."""
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            return self.copy_with_changes(location=location)
        else:
            return self

    def evaluate(self, problem):
        """ Return the sum of all codons frequencies.

        Note: no smart localization currently, the sequence is improved via

        """
        return {
            'best_codon': self.evaluate_best_codon,
            'harmonized': self.evaluate_harmonized,

        }[self.mode](problem)

    @staticmethod
    def codon_harmonization_stats(sequence, species):
        """Return a codon harmonisation score and a suboptimal locations list.

        Parameters
        ----------

        sequence
          An ATGC string

        species
          Any species name from the DnaChisel codon tables, such as ``e_coli``.

        Returns
        -------
        score, list_of_over_represented_codons_positions
          ``score`` is a negative number equals to sum(fi - ei) where for the
          i-th codon in the sequence fi is the relative frequency of this
          triplet in the sequence and ei is the relative frequency in the
          reference species. The ``list_of_suboptimal_codons_positions`` is
          of the form [1, 4, 5, 6...] a number k in that list indicates that
          the k-th codon is over-represented, and that a synonymous mutation
          of this codon can improve the harmonization score.

        """
        length = len(sequence)
        if (length % 3):
            raise ValueError(
                "Coding sequence with size %d not multiple of 3)" % length
            )
        codons_frequencies, codons_positions = \
            codons_frequencies_and_positions(sequence)
        codon_usage = CODON_USAGE_BY_AA[species]

        score = 0
        nonoptimal_indices = []
        for aa, usage_data in codon_usage.items():
            sequence_data = codons_frequencies[aa]
            for codon, usage_freq in usage_data.items():
                sequence_freq = sequence_data.get(codon, 0)
                frequency_diff = abs(sequence_freq - usage_freq)
                score -= frequency_diff * sequence_data['total']
                if sequence_freq > usage_freq:
                    nonoptimal_indices += codons_positions[codon]
        return score, nonoptimal_indices

    def codons_indices_to_locations(self, indices):
        """Convert a list of codon positions to a list of Locations"""
        indices = np.array(indices)
        if self.location.strand == -1:
            indices = sorted(self.location.end - indices)
            return [
                Location(group[0] - 3, group[-1], strand=-1)
                for group in group_nearby_indices(
                    indices,
                    max_group_spread=self.localization_group_spread)
            ]
        else:
            indices += self.location.start
            return [
                Location(group[0], group[-1] + 3)
                for group in group_nearby_indices(
                    indices,
                    max_group_spread=self.localization_group_spread)
            ]

    def evaluate_best_codon(self, problem):
        """Return the evaluation for mode==best_codon."""
        subsequence = self.location.extract_sequence(problem.sequence)
        length = len(subsequence)
        if (length % 3):
            raise ValueError(
                "CodonOptimizationSpecification on a window/sequence"
                "with size %d not multiple of 3)" % length
            )
        codons = [
            subsequence[3 * i: 3 * (i + 1)]
            for i in range(int(length / 3))
        ]
        CT = CODONS_TRANSLATIONS
        current_usage, optimal_usage = [np.array(e) for e in zip(*[
            (self.codon_usage_table[codon],
             self.codon_usage_table['best_frequencies'][CT[codon]])
            for codon in codons
        ])]
        non_optimality = optimal_usage - current_usage
        nonoptimal_indices = 3 * np.nonzero(non_optimality)[0]
        locations = self.codons_indices_to_locations(nonoptimal_indices)
        score = -non_optimality.sum()
        return SpecEvaluation(
            self, problem, score=score, locations=locations,
            message="Codon opt. on window %s scored %.02E" %
                    (self.location, score)
        )
    def evaluate_harmonized(self, problem):
        """Return the evaluation for mode==harmonized."""
        subsequence = self.location.extract_sequence(problem.sequence)
        score, nonoptimal_indices = self.codon_harmonization_stats(subsequence,
                                                                   self.species)
        locations = self.codons_indices_to_locations(nonoptimal_indices)
        np.random.shuffle(locations)
        return SpecEvaluation(
            self, problem, score=score, locations=locations,
            message="Codon opt. on window %s scored %.02E" %
                    (self.location, score)
        )

    def localized_on_window(self, new_location, start_codon, end_codon):
        """Relocate without changing much."""
        if self.mode == 'harmonized':
            return self
        else:
            return self.__class__(species=self.species, location=new_location,
                                  boost=self.boost)

    def label_parameters(self):
        return [self.species]
