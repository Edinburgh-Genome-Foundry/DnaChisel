"""Built-in genetic specifications."""

from .AllowPrimer import AllowPrimer
from .AvoidBlastMatches import AvoidBlastMatches
from .AvoidChanges import AvoidChanges
from .AvoidHairpins import AvoidHairpins
from .AvoidHeterodimerization import AvoidHeterodimerization
from .AvoidMatches import AvoidMatches
from .UniquifyAllKmers import UniquifyAllKmers
from .AvoidPattern import AvoidPattern
from .AvoidStopCodons import AvoidStopCodons
from .EnforceChanges import EnforceChanges
from .EnforceChoice import EnforceChoice
from .EnforceSequence import EnforceSequence
from .EnforceGCContent import EnforceGCContent
from .EnforceMeltingTemperature import EnforceMeltingTemperature
from .EnforcePatternOccurence import EnforcePatternOccurence
from .EnforceTranslation import EnforceTranslation
from .EnforceRegionsCompatibility import EnforceRegionsCompatibility
from .EnforceTerminalGCContent import EnforceTerminalGCContent
from .SequenceLengthBounds import SequenceLengthBounds

from .codon_optimization import (
    AvoidRareCodons,
    CodonOptimize,
    HarmonizeRCA,
    MaximizeCAI,
    MatchTargetCodonUsage,
)


DEFAULT_SPECIFICATIONS_DICT = {
    "AvoidRareCodons": AvoidRareCodons,
    "AllowPrimer": AllowPrimer,
    "AvoidBlastMatches": AvoidBlastMatches,
    "AvoidChanges": AvoidChanges,
    "AvoidHairpins": AvoidHairpins,
    "UniquifyAllKmers": UniquifyAllKmers,
    "AvoidPattern": AvoidPattern,
    "CodonOptimize": CodonOptimize,
    "EnforceGCContent": EnforceGCContent,
    "EnforcePatternOccurence": EnforcePatternOccurence,
    "EnforceTranslation": EnforceTranslation,
    "EnforceRegionsCompatibility": EnforceRegionsCompatibility,
    "EnforceSequence": EnforceSequence,
    "EnforceChoice": EnforceChoice,
    "EnforceChanges": EnforceChanges,
    "HarmonizeRCA": HarmonizeRCA,
    "MaximizeCAI": MaximizeCAI,
    "MatchTargetCodonUsage": MatchTargetCodonUsage
}

# Add the shorthands to the specifications dict.
for spec in list(DEFAULT_SPECIFICATIONS_DICT.values()):
    if spec.__dict__.get("shorthand_name", None) is not None:
        DEFAULT_SPECIFICATIONS_DICT[spec.shorthand_name] = spec

__all__ = [
    "AllowPrimer",
    "AvoidBlastMatches",
    "AvoidChanges",
    "AvoidHairpins",
    "AvoidHeterodimerization",
    "AvoidMatches",
    "UniquifyAllKmers",
    "AvoidPattern",
    "AvoidStopCodons",
    "CodonOptimize",
    "EnforceChanges",
    "EnforceChoice",
    "EnforceSequence",
    "EnforceGCContent",
    "EnforceMeltingTemperature",
    "EnforcePatternOccurence",
    "EnforceTranslation",
    "EnforceRegionsCompatibility",
    "EnforceTerminalGCContent",
    "EnforceSequence",
    "SequenceLengthBounds",
    "DEFAULT_SPECIFICATIONS_DICT",
]
