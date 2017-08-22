import numpy as np

from .CodonSpecification import CodonSpecification
from ..SpecEvaluation import SpecEvaluation
from ..biotools import (CODON_USAGE_TABLES, CODONS_TRANSLATIONS,
                        group_nearby_indices)
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

    def __init__(self, species=None, location=None,
                 codon_usage_table=None, boost=1.0):
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
        if self.location.strand == -1:
            nonoptimal_indices = sorted(self.location.end - nonoptimal_indices)
            locations = [
                Location(group[0] - 3, group[-1], strand=-1)
                for group in group_nearby_indices(
                    nonoptimal_indices,
                    max_group_spread=self.localization_group_spread)
            ]
        else:
            nonoptimal_indices += self.location.start
            locations = [
                Location(group[0], group[-1] + 3)
                for group in group_nearby_indices(
                    nonoptimal_indices,
                    max_group_spread=self.localization_group_spread)
            ]
        score = -non_optimality.sum()
        return SpecEvaluation(
            self, problem, score=score, locations=locations,
            message="Codon opt. on window %s scored %.02E" %
                    (self.location, score)
        )

    def localized_on_window(self, new_location, start_codon, end_codon):
        """Relocate without changing much."""
        return self.__class__(species=self.species, location=new_location,
                              boost=self.boost)

    def label_parameters(self):
        return [self.species]
