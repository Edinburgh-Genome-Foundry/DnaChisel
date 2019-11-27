from ..CodonSpecification import CodonSpecification
from python_codon_tables import get_codons_table
import numpy as np
from dnachisel.Location import Location
from dnachisel.biotools import group_nearby_indices
from copy import deepcopy


class BaseCodonOptimizationClass(CodonSpecification):
    def __init__(
        self, species=None, location=None, codon_usage_table=None, boost=1.0,
        codons_usage_threshold = 0,
    ):
        self.boost = boost
        self.location = Location.from_data(location)
        self.species = species
        self.codon_usage_table = self.get_codons_table(
            species, codon_usage_table
        )
        if codons_usage_threshold > 0:
            self.remove_codons_below_threshold(codons_usage_threshold)

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

    def remove_codons_below_threshold(self,codons_usage_threshold):
        result_table = deepcopy(self.codon_usage_table)

        for aa in 'ACDEFGHIKLMNPQRSTVWY*':
            threshold = codons_usage_threshold
            
            # ensure every aa has a codon which frequency above threshold
            if max(result_table[aa].values()) < threshold:
                threshold = max(result_table[aa].values())

            for codon in result_table[aa]:
                if result_table[aa][codon] < threshold:
                    result_table[aa][codon] = 0
            freq_sum = sum(result_table[aa].values())
            if freq_sum != 1 :
                for codon in result_table[aa]:
                    result_table[aa][codon] =  result_table[aa][codon] /freq_sum
        self.codon_usage_table =  result_table

    def initialized_on_problem(self, problem, role):
        """Get location from sequence if no location provided."""
        return self._copy_with_full_span_if_no_location(problem)

    def codons_indices_to_locations(self, indices):
        """Convert a list of codon positions to a list of Locations"""
        indices = np.array(indices)
        if self.location.strand == -1:
            indices = sorted(self.location.end - indices)
            return [
                Location(group[0] - 3, group[-1], strand=-1)
                for group in group_nearby_indices(
                    indices, max_group_spread=self.localization_group_spread
                )
            ]
        else:
            indices += self.location.start
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
