from .EnforceTranslation import EnforceTranslation
from dnachisel.biotools import CODONS_TRANSLATIONS
from ..Location import Location

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

    def __init__(self, location=None, boost=1.0, translation='legacy'):
        if (location is not None):
            if (len(location) % 3) != 0:
                message = "Loc. %s for AvoidStopCodon is not 3x" % location
                raise ValueError(message)
            else:
                self.translation = '_' * int(len(location) / 3)
        else:
            self.translation = None
        self.boost = boost
        self.location = location

    def initialize_on_problem(self, problem, role):
        """Get translation from the sequence if it is not already set."""
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            result = self.copy_with_changes()
            result.set_location(location)
            result.translation = '_' * int(len(location) / 3)
        else:
            result = self
        return result

    def __str__(self):
        """Represent."""
        return "AvoidStopCodons(%s)" % self.location

    def __str__(self):
        """Represent."""
        return "AvoidStopCodons(%s)" % self.location
