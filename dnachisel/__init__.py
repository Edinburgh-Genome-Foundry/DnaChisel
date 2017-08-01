from .DnaOptimizationProblem import (DnaOptimizationProblem,
                                     NoSolutionFoundError,
                                     NoSolutionError)
from .Location import Location

from .specifications import (
    AvoidBlastMatches,
    AvoidHairpins,
    AvoidNonuniqueSegments,
    AvoidPattern,
    AvoidStopCodon,
    CodonOptimize,
    AvoidChanges,
    EnforceGCContent,
    EnforcePattern,
    EnforceRegionsCompatibility,
    EnforceTranslation,
    Specification,
    SpecEvaluation,
    PatternSpecification,
    PROVIDERS_CONSTRAINTS,
    SequenceLengthBounds,
    VoidSpecification,
    # MinimizeDifferences,
    DEFAULT_SPECIFICATIONS_DICT
)

from .biotools import (
    blast_sequence,
    complement,
    DnaNotationPattern,
    crop_record,
    enzyme_pattern,
    homopolymer_pattern,
    random_dna_sequence,
    random_protein_sequence,
    repeated_kmers,
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

from .plotting_tools import (SpecAnnotationsTranslator,
                             plot_local_gc_content,
                             plot_local_gc_with_features)

from .utils import random_compatible_dna_sequence

from .version import __version__
