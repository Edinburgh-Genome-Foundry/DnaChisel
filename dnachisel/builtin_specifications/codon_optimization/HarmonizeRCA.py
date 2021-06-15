import numpy as np

from ...Specification.SpecEvaluation import SpecEvaluation
from .BaseCodonOptimizationClass import BaseCodonOptimizationClass


class HarmonizeRCA(BaseCodonOptimizationClass):
    """Codon-Harmonize a native sequence for a new host (Claassens method).

    This specification will optimize a Sequence 1 from Host 1 into a Sequence
    2 for target Host 2.

    In simple, rare Host 1 codons will be replaced by rare Host 2 codons, and
    high-frequency Host 1 codons will get replaced by codons that are
    high-frequency in Host 2.

    In more specific, each codon along Sequence 1 gets replaced by the codon
    whose Relative Codon Adaptiveness (RCA) in Host 2 is the closest from the
    RCA of the original codon in Host 1. A codon's RCA in a given organism is
    defined by f/fmax where f is the codon's frequency in the organism and fmax
    is the highest frequency of all synonymous codons.

    The minimized quantity is sum_i abs(RCA(c_i, H1) - RCA(c'_i, H2))
    where c_i, c'_i represent the i-th codon before and after optimization

    This method is taken from Claassens 2017, where they simplify a previous
    algorithm (Angov 2008), which was much more complicated as it involved
    predicting "ribosome pausing" sites in the sequence.

    Warning: always use with an EnforceTranslation constraint.


    Parameters
    ----------
    species
      Name or TaxID of the species for which to optimize the sequence. A custom
      codon_usage_table can be provided instead (or in addition, for species
      names whose codon usage table cannot be imported).

    codon_usage_table
      Optional - can be provided instead of ``species``. A dict of the form
      ``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}`` giving the RSCU table
      (relative usage of each codon).

    original_species
      Name or TaxID of the species the original sequence was taken from. This
      information will be used to spot codons which are supposed to be rare
      or common. A codon_usage_table can be provided instead (or in addition,
      for species names whose codon usage table cannot be imported).

    original_codon_usage_table
      A dict of the form ``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}``
      giving the RSCU table (relative usage of each codon).

    location
      Location on which the specification applies

    boost
      Score multiplicator (=weight) for when the specification is used as an
      optimization objective alongside competing objectives.


    References
    ----------
    Claassens et. al., Improving heterologous membrane protein
    production in Escherichia coli by combining transcriptional tuning and
    codon usage algorithms. PLOS One, 2017
    """

    shorthand_name = "harmonize_rca"

    def __init__(
        self,
        species=None,
        codon_usage_table=None,
        original_species=None,
        original_codon_usage_table=None,
        location=None,
        boost=1,
    ):
        if isinstance(species, str) and "->" in species:
            original_species, species = species.split("->")
            species = species.strip()
            original_species = original_species.strip()
        BaseCodonOptimizationClass.__init__(
            self,
            species=species,
            codon_usage_table=codon_usage_table,
            location=location,
            boost=boost,
        )
        self.codons_synonyms = self.get_codons_synonyms()
        self.original_species = original_species
        self.original_codon_usage_table = self.get_codons_table(
            original_species, original_codon_usage_table
        )
        for table in [self.codon_usage_table, self.original_codon_usage_table]:
            if "RCA" not in table:
                table["RCA"] = {
                    codon: frequency / max(codons_frequencies.values())
                    for aa, codons_frequencies in table.items()
                    for codon, frequency in codons_frequencies.items()
                    if len(aa) == 1
                }

    def initialized_on_problem(self, problem, role):
        new_spec = self._copy_with_full_span_if_no_location(problem)
        new_spec.original_codons = new_spec.get_codons(problem)
        rca = new_spec.codon_usage_table["RCA"]
        rca_o = new_spec.original_codon_usage_table["RCA"]
        new_spec.smallest_possible_discrepancies = [
            min([abs(rca[c] - rca_o[c]) for c in self.codons_synonyms[codon]])
            for codon in new_spec.original_codons
        ]
        return new_spec

    def evaluate(self, problem):
        """Return the evaluation for mode==best_codon."""
        codons = self.get_codons(problem)

        if len(codons) == 1:
            # We are evaluating a single codon. Easy!
            codon = codons[0]
            original = self.original_codons[0]
            rca_codon = self.codon_usage_table["RCA"][codon]
            rca_original = self.original_codon_usage_table["RCA"][original]
            score = -abs(rca_codon - rca_original)
            return SpecEvaluation(
                self,
                problem,
                score=score,
                locations=[] if (score == 0) else [self.location],
                message="Codon harmonization on window %s scored %.02E"
                % (self.location, score),
            )
        # print (len(codons))
        rca_in_original_species = [
            self.original_codon_usage_table["RCA"][original_codon]
            for original_codon in self.original_codons
        ]
        rca_in_target_species = [
            self.codon_usage_table["RCA"][codon] for codon in codons
        ]
        discrepancies = abs(
            np.array(rca_in_original_species) - np.array(rca_in_target_species)
        )
        non_optimality = self.smallest_possible_discrepancies - discrepancies
        nonoptimal_indices = np.nonzero(non_optimality)[0]
        locations = self.codons_indices_to_locations(nonoptimal_indices)
        score = -discrepancies.sum()
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=locations,
            message="Codon harmonization on %s scored %.02E" % (self.location, score),
        )

    def label_parameters(self):
        if self.species is None:
            return ["(custom table)"]
        else:
            return [self.original_species + " -> " + self.species]

    def short_label(self):
        result = "best-codon-optimize"
        if self.species is not None:
            result += " (%s)" % self.species
        return result
