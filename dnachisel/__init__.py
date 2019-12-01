from .DnaOptimizationProblem import (
    DnaOptimizationProblem,
    CircularDnaOptimizationProblem,
    NoSolutionError,
)

from .Location import Location

from .Specification import (
    Specification,
    SpecificationSet,
    SpecEvaluation
)

from .builtin_specifications import (
    AllowPrimer,
    AvoidBlastMatches,
    AvoidChanges,
    AvoidHairpins,
    AvoidMatches,
    AvoidRareCodons,
    UniquifyAllKmers,
    AvoidPattern,
    AvoidStopCodons,
    CodonOptimize,
    MaximizeCAI,
    HarmonizeRCA,
    MatchTargetCodonUsage,
    EnforceChanges,
    EnforceChoice,
    EnforceGCContent,
    EnforceMeltingTemperature,
    EnforcePatternOccurence,
    EnforceRegionsCompatibility,
    EnforceSequence,
    EnforceTerminalGCContent,
    EnforceTranslation,
    SequenceLengthBounds,
    DEFAULT_SPECIFICATIONS_DICT,
)

from .SequencePattern import (
    SequencePattern,
    DnaNotationPattern,
    HomopolymerPattern,
    RepeatedKmerPattern,
    EnzymeSitePattern,
    MotifPssmPattern
)

from .biotools import (
    blast_sequence,
    complement,
    list_common_enzymes,
    load_record,
    write_record,
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

Specification.DEFAULT_SPECIFICATIONS_DICT = DEFAULT_SPECIFICATIONS_DICT

__all__ = [
    "DnaOptimizationProblem",
    "NoSolutionError",
    "CircularDnaOptimizationProblem",
    "Location",
    "AllowPrimer",
    "AvoidBlastMatches",
    "AvoidChanges",
    "AvoidHairpins",
    "AvoidMatches",
    "AvoidRareCodons",
    "UniquifyAllKmers",
    "AvoidPattern",
    "AvoidStopCodons",
    "CodonOptimize",
    "MaximizeCAI",
    "HarmonizeRCA",
    "MatchTargetCodonUsage",
    "EnforceChoice",
    "EnforceChanges",
    "EnforceGCContent",
    "EnforceMeltingTemperature",
    "EnforcePatternOccurence",
    "EnforceRegionsCompatibility",
    "EnforceSequence",
    "EnforceTerminalGCContent",
    "EnforceTranslation",
    "SequenceLengthBounds",
    "Specification",
    "SpecificationSet",
    "SpecEvaluation",
    "SequencePattern",
    "DnaNotationPattern",
    "HomopolymerPattern",
    "RepeatedKmerPattern",
    "EnzymeSitePattern",
    "MotifPssmPattern",
    "blast_sequence",
    "complement",
    "list_common_enzymes",
    "load_record",
    "write_record",
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
