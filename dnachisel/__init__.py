from .DnaOptimizationProblem import (
    DnaOptimizationProblem,
    NoSolutionError,
    DEFAULT_SPECIFICATIONS_DICT,
)
from .CircularDnaOptimizationProblem import CircularDnaOptimizationProblem
from .Location import Location

from .builtin_specifications import (
    AllowPrimer,
    AvoidBlastMatches,
    AvoidChanges,
    AvoidHairpins,
    AvoidNonUniqueSegments,
    AvoidPattern,
    AvoidStopCodons,
    CodonOptimize,
    EnforceChoice,
    EnforceGCContent,
    EnforceMeltingTemperature,
    EnforcePatternOccurence,
    EnforceRegionsCompatibility,
    EnforceSequence,
    EnforceTerminalGCContent,
    EnforceTranslation,
    SequenceLengthBounds,
)

from .Specification import Specification, SpecificationsSet
from .SpecEvaluation import SpecEvaluation

from .SequencePattern import (
    DnaNotationPattern,
    HomopolymerPattern,
    RepeatedKmerPattern,
    EnzymeSitePattern,
)

from .biotools import (
    blast_sequence,
    complement,
    list_common_enzymes,
    load_record,
    random_dna_sequence,
    random_protein_sequence,
    reverse_complement,
    reverse_translate,
    sequences_differences,
    translate,
    change_biopython_record_sequence,
    annotate_differences,
    annotate_pattern_occurrences,
    annotate_record,
    sequence_to_biopython_record,
    sequences_differences_segments,
)

from .utils import random_compatible_dna_sequence


from .version import __version__

DEFAULT_SPECIFICATIONS_DICT.update(
    {
        "AvoidBlastMatches": AvoidBlastMatches,
        "AvoidChanges": AvoidChanges,
        "AvoidHairpins": AvoidHairpins,
        "AvoidNonUniqueSegments": AvoidNonUniqueSegments,
        "AvoidPattern": AvoidPattern,
        "CodonOptimize": CodonOptimize,
        "EnforceGCContent": EnforceGCContent,
        "EnforcePatternOccurence": EnforcePatternOccurence,
        "EnforceTranslation": EnforceTranslation,
        "EnforceRegionsCompatibility": EnforceRegionsCompatibility,
        "EnforceSequence": EnforceSequence,
        "EnforceChoice": EnforceChoice,
        # SHORTHAND NOTATIONS
        "cds": EnforceTranslation,
        "choice": EnforceChoice,
        "gc": EnforceGCContent,
        "insert": EnforcePatternOccurence,
        "keep": AvoidChanges,
        "no": AvoidPattern,
        "tm": EnforceMeltingTemperature,
        "primer": AllowPrimer,
        "sequence": EnforceSequence,
    }
)


__all__ = [
    "DnaOptimizationProblem",
    "NoSolutionError",
    "CircularDnaOptimizationProblem",
    "Location",
    "AllowPrimer",
    "AvoidBlastMatches",
    "AvoidChanges",
    "AvoidHairpins",
    "AvoidNonUniqueSegments",
    "AvoidPattern",
    "AvoidStopCodons",
    "CodonOptimize",
    "EnforceChoice",
    "EnforceGCContent",
    "EnforceMeltingTemperature",
    "EnforcePatternOccurence",
    "EnforceRegionsCompatibility",
    "EnforceSequence",
    "EnforceTerminalGCContent",
    "EnforceTranslation",
    "SequenceLengthBounds",
    "Specification",
    "SpecificationsSet",
    "SpecEvaluation",
    "DnaNotationPattern",
    "HomopolymerPattern",
    "RepeatedKmerPattern",
    "EnzymeSitePattern",
    "blast_sequence",
    "complement",
    "list_common_enzymes",
    "load_record",
    "random_dna_sequence",
    "random_protein_sequence",
    "reverse_complement",
    "reverse_translate",
    "sequences_differences",
    "translate",
    "change_biopython_record_sequence",
    "annotate_differences",
    "annotate_pattern_occurrences",
    "annotate_record",
    "sequence_to_biopython_record",
    "sequences_differences_segments",
    "random_compatible_dna_sequence",
    "__version__",
    "DEFAULT_SPECIFICATIONS_DICT",
]
