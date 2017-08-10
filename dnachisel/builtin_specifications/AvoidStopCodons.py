from .EnforceTranslation import EnforceTranslation
from dnachisel.biotools import CODONS_TRANSLATIONS

class AvoidStopCodons(EnforceTranslation):
    """Do not introduce any new stop codon in that frame.

    This can be used for research purposes, to avoid breaking a reading frame
    when editing it with quasi-synonymous mutations.
    """
    codons_translations = {
        codon: "*" if (translation == '*') else "_"
        for codon, translation in CODONS_TRANSLATIONS.items()
    }
    codons_sequences = None
    enforced_by_nucleotide_restrictions = False

    def __init__(self, location, translation=None, boost=1.0):
        if (len(location) % 3) != 0:
            raise ValueError("Loc. %s for AvoidStopCodon is not 3x" % location)
        self.translation = '_' * int(len(location) / 3)
        self.boost = boost
        self.location = location

    def __str__(self):
        """Represent."""
        return "AvoidStopCodons(%s)" % self.location

    def __str__(self):
        """Represent."""
        return "AvoidStopCodons(%s)" % self.location
