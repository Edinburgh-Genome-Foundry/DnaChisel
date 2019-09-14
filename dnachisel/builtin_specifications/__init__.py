"""Built-in genetic specifications."""

from .AllowPrimer import AllowPrimer
from .AvoidBlastMatches import AvoidBlastMatches
from .AvoidChanges import AvoidChanges
from .AvoidHairpins import AvoidHairpins
from .AvoidNonUniqueSegments import AvoidNonUniqueSegments
from .AvoidPattern import AvoidPattern
from .AvoidStopCodons import AvoidStopCodons
from .CodonOptimize import CodonOptimize
from .EnforceChoice import EnforceChoice
from .EnforceSequence import EnforceSequence
from .EnforceGCContent import EnforceGCContent
from .EnforcePatternOccurence import EnforcePatternOccurence
from .EnforceTranslation import EnforceTranslation
from .EnforceRegionsCompatibility import EnforceRegionsCompatibility
from .EnforceTerminalGCContent import EnforceTerminalGCContent
from .SequenceLengthBounds import SequenceLengthBounds

__all__ = [
    AllowPrimer,
    AvoidBlastMatches,
    AvoidChanges,
    AvoidHairpins,
    AvoidNonUniqueSegments,
    AvoidPattern,
    AvoidStopCodons,
    CodonOptimize,
    EnforceChoice,
    EnforceSequence,
    EnforceGCContent,
    EnforcePatternOccurence,
    EnforceTranslation,
    EnforceRegionsCompatibility,
    EnforceTerminalGCContent,
    EnforceSequence,
    SequenceLengthBounds,
]
