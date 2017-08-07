"""Built-in genetic specifications."""


from .AvoidBlastMatches import AvoidBlastMatches
from .AvoidChanges import AvoidChanges
from .AvoidHairpins import AvoidHairpins
from .AvoidNonUniqueSegments import AvoidNonuniqueSegments
from .AvoidPattern import AvoidPattern
from .CodonOptimize import CodonOptimize
from .EnforceGCContent import EnforceGCContent
from .EnforcePattern import EnforcePattern
from .EnforceTranslation import EnforceTranslation
from .EnforceRegionsCompatibility import EnforceRegionsCompatibility
from .EnforceSequence import EnforceSequence
from .SequenceLengthBounds import SequenceLengthBounds

DEFAULT_SPECIFICATIONS_DICT = {
   'AvoidBlastMatches': AvoidBlastMatches,
   'AvoidChanges': AvoidChanges,
   'AvoidHairpins': AvoidHairpins,
   'AvoidNonUniqueSegments': AvoidNonuniqueSegments,
   'AvoidPattern': AvoidPattern,
   'CodonOptimize': CodonOptimize,
   'EnforceGCContent': EnforceGCContent,
   'EnforcePattern': EnforcePattern,
   'EnforceTranslation': EnforceTranslation,
   'EnforceRegionsCompatibility': EnforceRegionsCompatibility,
   'EnforceSequence': EnforceSequence,
   'SequenceLengthBounds': SequenceLengthBounds
}
