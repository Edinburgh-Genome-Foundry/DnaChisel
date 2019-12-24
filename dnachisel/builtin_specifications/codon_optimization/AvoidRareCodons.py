"Implement AvoidRareCodons."

from ...Specification import SpecEvaluation
from ...biotools import reverse_complement
from .BaseCodonOptimizationClass import BaseCodonOptimizationClass


class AvoidRareCodons(BaseCodonOptimizationClass):
    """Avoid the use of codons with low frequency.

    This can be seen as a "mild" form of codon optimization where only rare
    codons (which slow down protein synthesis) are considered.

    WARNING: Make sure to always use this specification with EnforceTranslation
    to preserve the amino-acid sequence.

    Shorthand for annotations: "no_rare_codons".

    Parameters
    -----------
    min_frequency
      Minimal frequency accepted for a given codon.

    species
      Name or TaxID of the species for which to optimize the sequence. A custom
      codon_usage_table can be provided instead (or in addition, for species
      names whose codon usage table cannot be imported).
    
    codon_usage_table
      Optional codon usage table of the species for which the sequence will be
      codon-optimized, which can be provided instead of ``species``. A dict of
      the form ``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}`` giving the RSCU
      table (relative usage of each codon). See parameter ``species`` above.

    location
      Either a DnaChisel Location or a tuple of the form (start, end, strand)
      or just (start, end), with strand defaulting to +1, indicating the
      position of the gene to codon-optimize. If not provided, the whole
      sequence is considered as the gene. The location should have a length
      that is a multiple of 3. The location strand is either 1 if the gene is
      encoded on the (+) strand, or -1 for antisense.

    boost
      Score multiplicator (=weight) for when the specification is used as an
      optimization objective alongside competing objectives.
    """

    best_possible_score = 0
    enforced_by_nucleotide_restrictions = True
    shorthand_name = "no_rare_codons"

    def __init__(
        self,
        min_frequency,
        species=None,
        codon_usage_table=None,
        location=None,
        boost=1.0,
    ):
        """Initialize."""

        BaseCodonOptimizationClass.__init__(
            self,
            species=species,
            codon_usage_table=codon_usage_table,
            location=location,
            boost=boost,
        )
        self.min_frequency = min_frequency
        self.codons_frequencies = {
            codon: freq
            for aa, aa_data in self.codon_usage_table.items()
            if len(aa) == 1
            for codon, freq in aa_data.items()
        }
        self.rare_codons = sorted(
            [
                codon
                for codon, frequency in self.codons_frequencies.items()
                if frequency < min_frequency
            ]
        )
        self.nonrare_codons = sorted(
            [
                codon
                for codon, frequency in self.codons_frequencies.items()
                if frequency >= min_frequency
            ]
        )

    def evaluate(self, problem):
        """Score is the sum of (freq - min_frequency) for all rare codons."""
        # Note: this method is actually very little used as this specification
        # class sets the enforced_by_nucleotide_restrictions attribute.
        codons = self.get_codons(problem)
        rare_codons_indices = [
            i for i, codon in enumerate(codons) if codon in self.rare_codons
        ]
        locations = self.codons_indices_to_locations(rare_codons_indices)
        score = (
            0
            if (len(locations) == 0)
            else sum(
                (self.codons_frequencies[codons[i]] - self.min_frequency)
                for i in rare_codons_indices
            )
        )
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=locations,
            message="All OK."
            if len(locations) == 0
            else "Rare codons at locations %s" % locations,
        )

    def restrict_nucleotides(self, sequence, location=None):
        nonrare_codons = list(self.nonrare_codons)
        if self.location.strand == -1:
            nonrare_codons = sorted(
                [reverse_complement(c) for c in nonrare_codons]
            )
        return [
            ((i, i + 3), nonrare_codons)
            for i in range(self.location.start, self.location.end, 3)
        ]

    def _params_string(self):
        """Parameters representation used in __repr__, __str__, etc."""
        return "%d%%, %s" % (100 * self.min_frequency, str(self.species))

    def __repr__(self):
        return "AvoidRareCodons(%s)" % self._params_string()

    def __str__(self):
        return "AvoidRareCodons(%s)" % self._params_string()

    def short_label(self):
        return "no_rare_codons(%s)" % self._params_string()
