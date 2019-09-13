"""Methods for enzymes operations."""

from .sequences_operations import all_iupac_variants
from Bio import Restriction


def list_common_enzymes(
    site_length=(6,), opt_temp=(37,), min_suppliers=1, site_unlike=()
):
    """Return a list of enzyme names with the given constraints.

    Parameters
    ----------

    site_length
      List of accepted site lengths (6, 4, ...)

    opt_temp
      List of accepted optimal temperatures for the enzyme

    min_suppliers
      Minimal number registered suppliers in the Biopython data. A minimum
      of 3 known suppliers returns the most common enzymes.

    site_unlike
      List of (ambiguous or unambiguous) DNA sequences that should NOT be
      recognized by the selected enzymes.
    """
    site_unlike = set(
        [
            variant
            for enzyme in site_unlike
            for variant in all_iupac_variants(
                Restriction.__dict__[enzyme].site
            )
        ]
    )

    def is_valid(enzyme_name):
        enzyme = Restriction.__dict__[enzyme_name]
        return (
            len(enzyme.site) in site_length
            and enzyme.opt_temp in opt_temp
            and len(enzyme.supplier_list()) >= min_suppliers
            and len(
                set(all_iupac_variants(enzyme.site)).intersection(site_unlike)
            )
            == 0
        )

    return [
        enzyme_name
        for enzyme_name in Restriction.AllEnzymes.elements()
        if is_valid(enzyme_name)
    ]
