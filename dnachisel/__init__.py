from .DnaOptimizationProblem import (DnaOptimizationProblem, NoSolutionError,
                                     DEFAULT_SPECIFICATIONS_DICT)
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
    EnforcePatternOccurence,
    EnforceRegionsCompatibility,
    EnforceSequence,
    EnforceTerminalGCContent,
    EnforceTranslation,
    SequenceLengthBounds
)

from .Specification import Specification, SpecificationsSet
from .SpecEvaluation import SpecEvaluation

from .SequencePattern import (
    DnaNotationPattern,
    HomopolymerPattern,
    RepeatedKmerPattern,
    EnzymeSitePattern
)

from .biotools import (
    blast_sequence,
    complement,
    is_palyndromic,
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
    sequences_differences_segments
)

from .utils import (random_compatible_dna_sequence,
                    plot_constraint_breaches,
                    plot_gc_content_breaches,
                    plot_sequence_manufacturability_difficulties)
                    
                    

from .version import __version__

DEFAULT_SPECIFICATIONS_DICT.update({
    'AvoidBlastMatches': AvoidBlastMatches,
    'AvoidChanges': AvoidChanges,
    'AvoidHairpins': AvoidHairpins,
    'AvoidNonUniqueSegments': AvoidNonUniqueSegments,
    'AvoidPattern': AvoidPattern,
    'CodonOptimize': CodonOptimize,
    'EnforceGCContent': EnforceGCContent,
    'EnforcePatternOccurence': EnforcePatternOccurence,
    'EnforceTranslation': EnforceTranslation,
    'EnforceRegionsCompatibility': EnforceRegionsCompatibility,
    'EnforceSequence': EnforceSequence,
    'EnforceChoice': EnforceChoice,

    # SHORTHAND NOTATIONS
    'cds': EnforceTranslation,
    'choice': EnforceChoice,
    'gc': EnforceGCContent,
    'insert': EnforcePatternOccurence,
    'keep': AvoidChanges,
    'no': AvoidPattern
})
