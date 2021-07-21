from .HarmonizeRCA import HarmonizeRCA
from .MatchTargetCodonUsage import MatchTargetCodonUsage
from .MaximizeCAI import MaximizeCAI


def CodonOptimize(
    species=None,
    method="use_best_codon",
    location=None,
    codon_usage_table=None,
    original_species=None,
    original_codon_usage_table=None,
    boost=1.0
):
    """Codon-optimize a coding sequence using a user-selected method.

    This pseudo-specification is actually a function which returns an instance
    of another specification class depending on the selected "method":

    - For method="use_best_codon", every codon will be replaced by the "best"
      (i.e. most frequent) synonymous codon in the target organism. This is
      equivalent to Codon Adaptation Index (CAI) optimization.
    - For method="match_codon_usage", the final sequence's codon usage will
      match as much as possible the codon usage profile of the target species
      (this method is used throughout the literature, see for instance Hale
      and Thomson 1998).
    - For method="harmonize_rca", Each codon will be replaced by a synonymous
      codon whose usage in the target organism matches the usage of the
      original codon in its host organism (as per Claassens 2017).

    Parameters
    ==========
    species
      Species for which the sequence will be codon-optimized.
      Either a TaxID (this requires a web connection as the corresponding table
      will be downloaded from the internet) or the name of the species to
      codon-optimize for (the name must be supported by ``python_codon_tables``
      e.g. ``e_coli``, ``s_cerevisiae``, ``h_sapiens``, ``c_elegans``,
      ``b_subtilis``, ``d_melanogaster``).
      Note that a ``codon_usage_table`` can be provided instead, or even in
      addition, for species whose codon usage table cannot be auto-imported.

    method
      Either 'use_best_codon', 'match_codon_usage', or 'harmonize_rca'
      (see above for details)

    location
      Either a DnaChisel Location or a tuple of the form (start, end, strand)
      or just (start, end), with strand defaulting to +1, indicating the
      position of the gene to codon-optimize. If not provided, the whole
      sequence is considered as the gene. The location should have a length
      that is a multiple of 3. The location strand is either 1 if the gene is
      encoded on the (+) strand, or -1 for antisense.

    codon_usage_table
      Optional codon usage table of the species for which the sequence will be
      codon-optimized, which can be provided instead of ``species``. A dict of
      the form ``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}`` giving the RSCU
      table (relative usage of each codon). See parameter ``species`` above.

    original_species
      When the method is 'harmonize_rca', this is the native species of the
      original coding sequence. Same characteristics as parameter ``species``
      above.
  
    original_codon_usage_table
      Optional codon usage table of the original sequence's native species.
      A dict of the form ``{'*': {"TGA": 0.112, "TAA": 0.68}, 'K': ...}``
      giving the codon usage table.
  
    References
    ==========
  
    Claassens et. al., Improving heterologous membrane protein
    production in Escherichia coli by combining transcriptional tuning and
    codon usage algorithms. PLOS One, 2017

    Hale and Thompson, Codon Optimization of the Gene Encoding a
    Domain from Human Type 1 Neurofibromin Protein... Protein Expression and
    Purification 1998.

    """
    if method == "use_best_codon":
        return MaximizeCAI(
            species=species,
            location=location,
            codon_usage_table=codon_usage_table,
            boost=boost,
        )

    elif method == "match_codon_usage":
        return MatchTargetCodonUsage(
            species=species,
            location=location,
            codon_usage_table=codon_usage_table,
            boost=boost,
        )
    elif method == "harmonize_rca":
        return HarmonizeRCA(
            species=species,
            location=location,
            codon_usage_table=codon_usage_table,
            original_species=original_species,
            original_codon_usage_table=original_codon_usage_table,
            boost=boost,
        )
    raise ValueError("Parameter 'mode' should be one of best_codon, "
                     "match_usage, ")
