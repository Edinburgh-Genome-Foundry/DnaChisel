import numpy as np

from .BaseCodonOptimizationClass import BaseCodonOptimizationClass
from ...Specification.SpecEvaluation import SpecEvaluation


class MaximizeCAI(BaseCodonOptimizationClass):
    """Codon-optimize a coding sequence for a given species. Maximizes the CAI.

    To be precise, the score computed by this specification is N*log(CAI) where
    N is the number of codons. Maximizing this score also maximizes the CAI.

    Index (CAI). For a sequence with N codons, the CAI is the geometric mean
    of the Relative Codon Adaptiveness (RCA) of the different codons. The RCA
    of a codon is (f_i/fmax_i) were fi is the frequency of an oligo in the
    codon usage table, and fmax is the maximal frequency of the synonymous
    codons.

    So N*log(CAI) = sum_i ( log(f_i) - log(fmax_i) )

    This score is between -inf. and 0 (0 meaning a perfectly optimal sequence).

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

    Examples
    --------

    >>> objective = MaximizeCAI(
    >>>     species = "E. coli",
    >>>     location = (150, 300), # coordinates of a gene
    >>>     strand = -1
    >>> )


    """

    shorthand_name = "use_best_codon"

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
        if "log_best_frequencies" not in self.codon_usage_table:
            self.codon_usage_table["log_best_frequencies"] = {
                aa: np.log(max(aa_data.values()))
                for aa, aa_data in self.codon_usage_table.items()
                if len(aa) == 1
            }
        if "log_codons_frequencies" not in self.codon_usage_table:
            self.codon_usage_table["log_codons_frequencies"] = {
                codon: np.log(frequency or 0.001)
                for aa, frequencies in self.codon_usage_table.items()
                for codon, frequency in frequencies.items()
                if len(aa) == 1
            }

    def evaluate(self, problem):
        """Evaluate!"""
        codons = self.get_codons(problem)
        ct = self.codons_translations
        if len(codons) == 1:
            # We are evaluating a single codon. Easy!
            codon = codons[0]
            freq = self.codon_usage_table["log_codons_frequencies"][codon]
            optimal = self.codon_usage_table["log_best_frequencies"][ct[codon]]
            score = freq - optimal
            return SpecEvaluation(
                self,
                problem,
                score=freq - optimal,
                locations=[] if (freq == optimal) else [self.location],
                message="Codon opt. on window %s scored %.02E"
                % (self.location, score),
            )
        current_usage = [
            self.codon_usage_table["log_codons_frequencies"][codon]
            for codon in codons
        ]
        optimal_usage = [
            self.codon_usage_table["log_best_frequencies"][ct[codon]]
            for codon in codons
        ]
        non_optimality = np.array(optimal_usage) - np.array(current_usage)
        nonoptimal_indices = np.nonzero(non_optimality)[0]
        locations = self.codons_indices_to_locations(nonoptimal_indices)
        score = -non_optimality.sum()
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=locations,
            message="Codon opt. on window %s scored %.02E"
            % (self.location, score),
        )

    def label_parameters(self):
        return ["(custom table)" if self.species is None else self.species]

    def short_label(self):
        result = "best-codon-optimize"
        if self.species is not None:
            result += " (%s)" % self.species
        return result

