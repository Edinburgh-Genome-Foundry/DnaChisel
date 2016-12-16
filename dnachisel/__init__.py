from .DnaOptimizationProblem import (DnaOptimizationProblem,
                                     NoSolutionFoundError)
from .Location import Location
import inspect

from .objectives import (
    AvoidBlastMatches,
    AvoidNonuniqueSegments,
    AvoidPattern,
    CodonOptimize,
    DoNotModify,
    EnforceGCContent,
    EnforcePattern,
    EnforceRegionsCompatibility,
    EnforceTranslation,
    Objective,
    ObjectiveEvaluation,
    PatternObjective,
    PROVIDERS_CONSTRAINTS,
    SequenceLengthBounds,
    VoidObjective,
    MinimizeDifferences
)

from .biotools import (
    blast_sequence,
    DnaNotationPattern,
    enzyme_pattern,
    homopolymer_pattern,
    random_dna_sequence,
    random_protein_sequence,
    reverse_complement,
    reverse_translate,
    sequences_differences,
    translate,
    change_biopython_record_sequence
)

objectives_dict = {
   k: v
   for (k, v) in locals().items()
   if inspect.isclass(v) and issubclass(v, Objective)
}

from .utils import random_compatible_dna_sequence

from .version import __version__
