from ..Specification import Specification, SpecEvaluation
from ..Location import Location

from .tools import analyze_cai

class tailor_Specification(Specification):

    pass

class CAI(tailor_Specification):


    def __init__(self, location=None, boost=1.0):
        """Initialize."""
        self.location = Location.from_data(location)
        self.boost = boost

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Find matches and count them negatively"""
        location = self.location
        if location is None:
            location = Location(0, len(problem.sequence))
        sequence = location.extract_sequence(problem.sequence)
        

        score = analyze_cai(sequence)

        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=[location],
            message="CAI: %s at %s" % (score, location),
        )

class GCContent(tailor_Specification):
    def __init__(self, location=None, boost=1.0):
        """Initialize."""
        self.location = Location.from_data(location)
        self.boost = boost

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Find matches and count them negatively"""
        location = self.location
        if location is None:
            location = Location(0, len(problem.sequence))
        sequence = location.extract_sequence(problem.sequence)
        
        
        num_G = sequence.count('G')
        num_C = sequence.count('C')
        
        score = (num_C + num_G) / len(sequence)

        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=[location],
            message="GC: %s at %s" % (score, location),
        )
    