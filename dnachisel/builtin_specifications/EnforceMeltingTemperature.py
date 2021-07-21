try:
    import primer3

    PRIMER3_AVAILABLE = True
except ImportError:
    primer3 = None
    PRIMER3_AVAILABLE = False

from Bio.SeqUtils import MeltingTemp as bio_mt
from ..Location import Location
from ..Specification import Specification, SpecEvaluation



class EnforceMeltingTemperature(Specification):
    """Ensure that the subsequence's Tm is in a certain segment/target.

    Shorthand for annotations: "tm".

    Parameters
    ----------
    mini, maxi
      Minimum and maximum acceptable melting temperatures in Celcius, for
      instance 55 and 70. A "target" can be provided instead when using this
      specification as an optimization objective.
    target
      Target melting temperature. Will be overridden by (mini+maxi)/2 if these
      are provided. The "target" parameter is only practical when the spec is
      used as an optimization objective.

    location
      Location of the subsequence whose melting temperature is considered.
      Can be None if the whole sequence is to be considered.

    boost
      Multiplicator for this specification's score when used in a
      multi-objective optimization.
    """

    shorthand_name = "tm"

    def __init__(
        self, mini=None, maxi=None, target=None, location=None, boost=1.0
    ):
        """Initialize."""
        if isinstance(mini, str) and mini.endswith('C'):
            # PROCESS CASES "45-55%" and "45%"
            split = mini[:-1].split('-')
            if len(split) == 2:
                mini, maxi = float(split[0]), float(split[1])
            else:
                target = float(split[0])
        if target is not None:
            mini = maxi = target
        else:
            target = 0.5 * (mini + maxi)
        self.mini = mini
        self.maxi = maxi
        self.target = target
        self.location = Location.from_data(location)
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Return the sum of breaches extent for all windowed breaches."""
        sequence = self.location.extract_sequence(problem.sequence)
        predictor = primer3.calcTm if PRIMER3_AVAILABLE else bio_mt.Tm_NN
        tm = predictor(sequence)
        score = 0.5 * (self.maxi - self.mini) - abs(tm - self.target)
        return SpecEvaluation(
            specification=self,
            problem=problem,
            score=score,
            locations=[self.location],
            message="Tm = %.1f " % tm,
        )
