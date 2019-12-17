import re
from Bio.Restriction.Restriction_Dictionary import rest_dict
from .DnaNotationPattern import DnaNotationPattern

class EnzymeSitePattern(DnaNotationPattern):
    """Class to represent Enzyme site patterns

    Examples
    --------

    >>> enzyme_pattern = EnzymeSitePattern("BsaI")
    >>> constraint = AvoidPattern(enzyme_pattern)

    """

    def __init__(self, enzyme_name):
        self.enzyme_site = rest_dict[enzyme_name]["site"]
        DnaNotationPattern.__init__(self, self.enzyme_site, name=enzyme_name)

    @staticmethod
    def from_string(string):
        """Convert BsmBI_site to EnzymeSitePattern(BsmBI)"""
        match = re.match(r"(\S+)_site", string)
        if match is not None:
            enzyme_name = match.groups()[0]
            if enzyme_name in rest_dict:
                return EnzymeSitePattern(enzyme_name)

    def __str__(self):
        return "%s(%s)" % (self.name, self.enzyme_site)