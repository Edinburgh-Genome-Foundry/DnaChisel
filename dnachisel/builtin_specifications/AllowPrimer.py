try:
    import primer3
    PRIMER3_AVAILABLE = True
except:
    primer3 = None
    PRIMER3_AVAILABLE = False

from Bio.SeqUtils import MeltingTemp as bio_mt
from ..Specification import Specification, SpecificationsSet
from ..SpecEvaluation import SpecEvaluation
from ..Location import Location
from ..SequencePattern import RepeatedKmerPattern
# from .VoidSpecification import VoidSpecification
from .AvoidNonUniqueSegments import AvoidNonUniqueSegments
from .AvoidPattern import AvoidPattern

class AvoidHeterodimerization(Specification):
    """Avoid that the (sub)sequence anneals with other primers.

    Parameters
    ----------

    other_primers_sequences
      List of ATGC strings representing the sequences of the other primers
      with which annealing will be avoided.

    tmax
      Maximal melting temperature of the interaction between the (sub)sequence
      and the other_primers_sequences.
    
    location
      Location of the subsequence to be optimized.
    
    boost
      Multiplicator for this specification's score when used in a
      multi-objective optimization.
    """
    def __init__(self, other_primers_sequences, tmax=5, location=None,
                 boost=1.0):
        self.other_primers_sequences = other_primers_sequences
        self.tmax = tmax
        self.location = location
    
    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)
    
    def localize(self, location, problem=None):
        if self.location.overlap_region(location) is None:
            return None 
 # VoidSpecification(parent_specification=self)
        else:
            return self
    
    def evaluate(self, problem):
        if not PRIMER3_AVAILABLE:
            raise ImportError("Using avoid_heterodimerization requires primer3"
                              " installed (pip install primer3-py)")
        sequence = self.location.extract_sequence(problem.sequence)
        melting_temps = [
            primer3.calcHeterodimer(sequence, other_seq).tm
            for other_seq in self.other_primers_sequences
        ]
        largest_tm = max(melting_temps)
        # hackish penalty to guide optimization:
        penalty = 0.001 * sum(melting_temps) / len(melting_temps)
        score = self.tmax - largest_tm - penalty
        return SpecEvaluation(specification=self, problem=problem, score=score,
                              locations=[self.location],
                              message="Largest Tm = %.1f " % largest_tm)

    def label_parameters(self):
        return [('primers', len(self.other_primers_sequences)),
                ('tmax', self.tmax)]

class EnforceMeltingTemperature(Specification):
    """Ensure that the subsequence's Tm is in a certain segment/target.

    Parameters
    ----------

    mini, maxi
      Minimum and maximum acceptable melting temperatures in Celcius, for
      instance 55 and 70. A "target" can be provided instead when using this
      specification as an optimization objective.

    target
      Target melting temperature. Will be overriden by (mini+maxi)/2 if these
      are provided. The "target" parametr is only practical when the spec is
      used as an optimization objective.
    
    location
      Location of the subsequence whose melting temperature is considered.
      Can be None if the whole sequence is to be considered.
    
    boost
      Multiplicator for this specification's score when used in a
      multi-objective optimization.
    """
    def __init__(self, mini=None, maxi=None, target=None, location=None,
                 boost=1.0):
        """Initialize."""
        if target is not None:
            mini = maxi = target
        else:
            target = 0.5 * (mini + maxi)
        self.mini = mini
        self.maxi = maxi
        self.target = target
        if isinstance(location, tuple):
            location = Location.from_tuple(location)
        self.location = location
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Return the sum of breaches extent for all windowed breaches."""
        sequence = self.location.extract_sequence(problem.sequence)
        predictor = primer3.calcTm if PRIMER3_AVAILABLE else bio_mt.Tm_NN
        tm = predictor(sequence)
        score = 0.5 * (self.maxi - self.mini) - abs(tm - self.target)
        return SpecEvaluation(specification=self, problem=problem, score=score,
                              locations=[self.location],
                              message="Tm = %.1f " % tm)

    
class AllowPrimer(SpecificationsSet):
    """Enforce various specifications for enabling primers at the location.

    This is useful for making sure that you will be able to conduct a PCR or
    sanger sequencing with a primer annealing at a particular location of
    the sequence.

    Parameters
    ----------

    location
      The exact location where the primer will anneal, i.e. the subsequence
      under this location will be the sequence

    tmin, tmax
      Minimum and maximum acceptable melting temperatures in Celcius, for
      instance 55 and 70. When the used as an optimization objective the
      "target" will be (tmin + tmax)/2.
    
    max_homology_length
      Maximal length of any homology between the subsequence at that location
      and anywhere else in the sequence  

    avoid_heterodim_with
      List of ATGC strings representing external (primer) sequences with which
      the optimised location should have no annealing.
    
    max_heterodim_tm
      Max melting temperature of the heterodimerization between this subsegment
      and the other sequences in "avoid_heterodim_with"
    
    avoided_repeats
      Properties of the repeated patterns avoided. List of pairs (K, N)
      meaning "avoid K-mers repeated N times in a row".
    """
    
    def __init__(self, location=None, tmin=50, tmax=70,
                 max_homology_length=6, avoid_heterodim_with=None,
                 max_heterodim_tm=5,
                 avoided_repeats=((2, 5), (3, 4), (4, 3))):
        if isinstance(location, tuple):
            location = Location.from_tuple(location)
        specs = {
            'unique_sequence': AvoidNonUniqueSegments(
                min_length=max_homology_length, location=location),
            'melting_temperature': EnforceMeltingTemperature(
                mini=tmin, maxi=tmax, location=location),
           **{
               'repeats_%d_%d' % (k, n): AvoidPattern(
                   RepeatedKmerPattern(k, n), location=location)
                for (k, n) in avoided_repeats
           }
        }
        if avoid_heterodim_with is not None:
            specs['avoid_heterodimerization'] = AvoidHeterodimerization(
                other_primers_sequences=avoid_heterodim_with,
                tmax=max_heterodim_tm,
                location=location
            )
        self.register_specifications(specs)