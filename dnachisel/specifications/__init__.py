import inspect
from .Specification import (Specification, VoidSpecification,
                            PatternSpecification)
from .SpecEvaluation import (SpecEvaluation,
                             ProblemObjectivesEvaluations,
                             ProblemConstraintsEvaluations)
# from .specifications_sets import PROVIDERS_CONSTRAINTS

from .builtin_specifications import (
    AvoidBlastMatches,
    AvoidHairpins,
    AvoidNonuniqueSegments,
    AvoidPattern,
    AvoidStopCodon,
    CodonOptimize,
    AvoidChanges,
    EnforceGCContent,
    EnforcePattern,
    EnforceTranslation,
    EnforceRegionsCompatibility,
    #MinimizeDifferences,
    SequenceLengthBounds
)

DEFAULT_SPECIFICATIONS_DICT = {
   k: v
   for (k, v) in locals().items()
   if inspect.isclass(v) and issubclass(v, Specification)
}
