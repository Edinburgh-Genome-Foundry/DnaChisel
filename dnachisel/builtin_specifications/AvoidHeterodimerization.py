try:
    import primer3

    PRIMER3_AVAILABLE = True
except (ImportError, OSError):
    primer3 = None

from ..Specification import Specification, SpecEvaluation


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

    def __init__(
        self, other_primers_sequences, tmax=5, location=None, boost=1.0
    ):
        self.other_primers_sequences = other_primers_sequences
        self.tmax = tmax
        self.location = location

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def localize(self, location, problem=None):
        if self.location.overlap_region(location) is None:
            return None
        else:
            return self

    def evaluate(self, problem):
        if not PRIMER3_AVAILABLE:
            raise ImportError(
                "Using avoid_heterodimerization requires primer3"
                " installed (pip install primer3-py)"
            )
        if len(self.other_primers_sequences) == 0:
            return SpecEvaluation(
                specification=self,
                problem=problem,
                score=0,
                locations=[self.location],
                message="No existing primer"
            )
        sequence = self.location.extract_sequence(problem.sequence)
        melting_temps = [
            primer3.calcHeterodimer(sequence, other_seq).tm
            for other_seq in self.other_primers_sequences
        ]
        largest_tm = max(melting_temps)
        # hackish penalty to guide optimization:
        penalty = 0.001 * sum(melting_temps) / len(melting_temps)
        score = self.tmax - largest_tm - penalty
        return SpecEvaluation(
            specification=self,
            problem=problem,
            score=score,
            locations=[self.location],
            message="Largest Tm = %.1f " % largest_tm,
        )

    def label_parameters(self):
        return [
            ("primers", len(self.other_primers_sequences)),
            ("tmax", self.tmax),
        ]
