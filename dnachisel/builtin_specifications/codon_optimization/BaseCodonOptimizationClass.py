from ..CodonSpecification import CodonSpecification
from python_codon_tables import get_codons_table
import numpy as np
from ...Location import Location
from ...biotools import group_nearby_indices


class BaseCodonOptimizationClass(CodonSpecification):

    best_possible_score = 0  # Don't forget to change in subclasses if needed
    localization_group_spread = 3

    def __init__(
        self, species=None, location=None, codon_usage_table=None, boost=1.0
    ):
        self.boost = boost
        self.location = Location.from_data(location)
        self.species = species
        self.codon_usage_table = self.get_codons_table(
            species, codon_usage_table
        )

    def get_codons(self, problem):
        subsequence = self.location.extract_sequence(problem.sequence)
        if len(subsequence) % 3:
            raise ValueError(
                "Spec. %s is on a window/sequence with size not multiple of 3)"
                % (self.label())
            )
        return [
            subsequence[3 * i : 3 * (i + 1)]
            for i in range(int(len(subsequence) / 3))
        ]

    @staticmethod
    def get_codons_table(species, codon_usage_table):
        if codon_usage_table is None:
            if species is None:
                raise ValueError(
                    "Provide either an species name or a codon usage table"
                )
            else:
                codon_usage_table = get_codons_table(species)
        return codon_usage_table

    def initialized_on_problem(self, problem, role):
        """Get location from sequence if no location provided."""
        return self._copy_with_full_span_if_no_location(problem)

    def codons_indices_to_locations(self, indices):
        """Convert a list of codon positions to a list of Locations"""
        indices = np.array(indices)
        if self.location.strand == -1:
            indices = sorted(self.location.end - 3 * indices)
            return [
                Location(group[0] - 3, group[-1], strand=-1)
                for group in group_nearby_indices(
                    indices, max_group_spread=self.localization_group_spread
                )
            ]
        else:
            indices = self.location.start + 3 * indices
            return [
                Location(group[0], group[-1] + 3)
                for group in group_nearby_indices(
                    indices, max_group_spread=self.localization_group_spread
                )
            ]

    def get_codons_synonyms(self):
        """Return a dict {"GTG": [GTG, GTC, ...]} of synonymous codons."""
        return {
            codon: [c for c in aa_codons]
            for aa, aa_codons in self.codon_usage_table.items()
            if len(aa) == 1
            for codon in aa_codons
        }

    def get_codons_translations(self):
        """Return a dict {"ATG": "M", "TAG": "*", ...}."""
        return {
            codon: aa
            for aa, aa_codons in self.codon_usage_table.items()
            if len(aa) == 1
            for codon in aa_codons.keys()
        }

    def localized_on_window(self, new_location, start_codon, end_codon):
        """Relocate without changing much."""
        # The "new_location" already has exactly the right span and strand
        # thanks to superclass CodonSpecification
        return self.copy_with_changes(location=new_location)
