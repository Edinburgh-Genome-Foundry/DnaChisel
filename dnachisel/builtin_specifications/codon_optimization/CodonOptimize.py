from .HarmonizeRCA import HarmonizeRCA
from .MatchTargetCodonUsage import MatchTargetCodonUsage
from .MaximizeCAI import MaximizeCAI


def CodonOptimize(
    species=None,
    mode="best_codon",
    location=None,
    codon_usage_table=None,
    original_species=None,
    original_codon_usage_table=None,
    boost=1.0,
):
    """Peudo-Specification that returns a subclassed Specification"""
    if mode == "best_codon":
        return MaximizeCAI(
            species=species,
            location=location,
            codon_usage_table=codon_usage_table,
            boost=boost,
        )

    elif mode == "match_usage":
        return MatchTargetCodonUsage(
            species=species,
            location=location,
            codon_usage_table=codon_usage_table,
            boost=boost,
        )
    elif mode == "harmonize_rca":
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
