import inspect
from .Objective import (Objective, VoidObjective, PatternObjective)
from .ObjectiveEvaluation import (ObjectiveEvaluation,
                                  ProblemObjectivesEvaluations,
                                  ProblemConstraintsEvaluations)
from .objectives_sets import PROVIDERS_CONSTRAINTS
from .objectives import (
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
    MinimizeDifferences,
    SequenceLengthBounds
)

DEFAULT_OBJECTIVES_DICT = {
   k: v
   for (k, v) in locals().items()
   if inspect.isclass(v) and issubclass(v, Objective)
}
