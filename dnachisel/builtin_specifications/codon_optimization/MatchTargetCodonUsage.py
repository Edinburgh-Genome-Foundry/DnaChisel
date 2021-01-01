import numpy as np
from ...Specification.SpecEvaluation import SpecEvaluation
from ...biotools import dict_to_pretty_string

from .BaseCodonOptimizationClass import BaseCodonOptimizationClass


class MatchTargetCodonUsage(BaseCodonOptimizationClass):
    """Codon-optimize a sequence so it has the same codon usage as a target.

    The objective minimized here is the sum of the discrepancies, over every
    possible triplet ATG, CCG, etc. between the codon frequency of this triplet
    in the sequence, and its frequency in the target organism.

    This method has had several names through the ages. It may have been first
    proposed by Hale and Thompson, 1998. It is called Individual Codon Usage
    Optimization in Chung 2012, Global CAI Harmonization in Mignon 2018, and
    Codon Harmonization in Jayaral 2005. We didn't call it "harmonization"
    in DNA Chisel to avoid any confusion with the now more common
    host-to-target codon harmonization. See DnaChisel's HarmonizeRCA class
    for Codon Harmonization.

    Parameters
    ----------

    species
      Species for which the sequence will be codon-optimized.
      Either a TaxID (this requires a web connection as the corresponding table
      will be downloaded from the internet) or the name of the species to
      codon-optimize for (the name must be supported by ``python_codon_tables``
      e.g. ``e_coli``, ``s_cerevisiae``, ``h_sapiens``, ``c_elegans``,
      ``b_subtilis``, ``d_melanogaster``).
      Note that a ``codon_usage_table`` can be provided instead, or even in
      addition, for species whose codon usage table cannot be auto-imported.

    location
      Either a DnaChisel Location or a tuple of the form (start, end, strand)
      or just (start, end), with strand defaulting to +1, indicating the
      position of the gene to codon-optimize. If not provided, the whole
      sequence is considered as the gene. The location should have a length
      that is a multiple of 3. The location strand is either 1 if the gene is
      encoded on the (+) strand, or -1 for antisense.

    codon_usage_table
      A dict of the form ``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}``
      giving the RSCU table (relative usage of each codon). Only provide if
      no ``species`` parameter was provided.

    boost
      Score multiplicator (=weight) for when the specification is used as an
      optimization objective alongside competing objectives.

    References
    ----------
    Hale and Thompson, Codon Optimization of the Gene Encoding a
    Domain from Human Type 1 Neurofibromin Protein... Protein Expression and
    Purification 1998.

    Jayaraj et. al. GeMS: an advanced software package for designing synthetic
    genes, Nucleic Acids Research, 2005

    Mignon et. al. Codon harmonization – going beyond the speed limit for
    protein expression. FEBS Lett, 2018

    Chung BK, Lee DY. Computational codon optimization of synthetic gene for
    protein expression. BMC Syst Biol. 2012


    """

    shorthand_name = "match_codon_usage"

    def __init__(
        self, species=None, location=None, codon_usage_table=None, boost=1.0
    ):
        BaseCodonOptimizationClass.__init__(
            self,
            species=species,
            location=location,
            codon_usage_table=codon_usage_table,
            boost=boost,
        )
        self.codons_translations = self.get_codons_translations()

    def codon_usage_matching_stats(self, problem):
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
        codons = self.get_codons(problem)
        codons_positions, aa_comparisons = self.compare_frequencies(codons)
        score = 0
        nonoptimal_aa_indices = []
        for aa, data in aa_comparisons.items():
            total = data.pop("total")
            for codon, codon_freq in data.items():
                frequency_diff = codon_freq["sequence"] - codon_freq["table"]
                score -= total * abs(frequency_diff)
                if codon_freq["sequence"] > codon_freq["table"]:
                    nonoptimal_aa_indices += codons_positions[codon]
        return score, nonoptimal_aa_indices

    def evaluate(self, problem):
        """Evaluate on a problem"""
        score, nonoptimal_indices = self.codon_usage_matching_stats(problem)
        locations = self.codons_indices_to_locations(nonoptimal_indices)
        np.random.shuffle(locations)
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=locations,
            message="Codon opt. on window %s scored %.02E"
            % (self.location, score),
        )

    def localized_on_window(self, new_location, start_codon, end_codon):
        """Relocate without changing much."""
        return self

    def label_parameters(self):
        return ["(custom table)" if self.species is None else self.species]

    def compare_frequencies(self, codons, text_mode=False):
        """Return a dict indicating differences between codons frequencies in
        the sequence and in this specifications's codons usage table.

        Examples
        --------

        >>> codons = spec.get_codons(problem)
        >>> print(spec.compare_frequencies(codons)

        Returns
        -------

        positions, comparisons
          (if text_mode = False)

        a formatted print-ready string
          (if text_mode = True)

        >>> {
        >>>   "K": {
        >>>     "total": 6,
        >>>     "AAA": {
        >>>         "sequence": 1.0,
        >>>         "table": 0.7
        >>>     },
        >>>     ...
        >>>   },
        >>>   "D": ...
        >>> }

        """
        codons_positions = {cod: [] for cod in self.codons_translations}
        for i, codon in enumerate(codons):
            codons_positions[codon].append(i)
        # aa: amino-acid
        codons_frequencies = {
            aa: {"total": 0} for aa in self.codon_usage_table
        }
        for codon, positions in codons_positions.items():
            count = len(positions)
            aa = self.codons_translations[codon]
            codons_frequencies[aa][codon] = count
            codons_frequencies[aa]["total"] += count
        for aa, data in codons_frequencies.items():
            total = max(1, data["total"])
            for codon, value in data.items():
                if codon != "total":
                    data[codon] = 1.0 * value / total
        codons_frequencies = {
            aa: data
            for aa, data in codons_frequencies.items()
            if data["total"]
        }
        comparisons = {
            aa: {
                "total": seq_data["total"],
                **{
                    codon: {"sequence": seq_data[codon], "table": table_data}
                    for codon, table_data in self.codon_usage_table[aa].items()
                },
            }
            for aa, seq_data in codons_frequencies.items()
        }
        if text_mode:
            return dict_to_pretty_string(comparisons)
        else:
            return codons_positions, comparisons
    def short_label(self):
        result = "match-codon-usage"
        if self.species is not None:
            result += " (%s)" % self.species
        return result