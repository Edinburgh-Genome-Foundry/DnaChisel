"""Implement AvoidPattern"""

from ..Specification import Specification, VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from .EnforceSequence import EnforceSequence
from ..MutationSpace import MutationSpace
from ..SequencePattern import DnaNotationPattern, enzyme_pattern
from dnachisel.Location import Location
from ..DnaOptimizationProblem import DnaOptimizationProblem, NoSolutionError


class EnforcePattern(Specification):
    """Enforce a number of occurences of the given pattern in the sequence.

    Parameters
    ----------
    pattern
      A SequencePattern or DnaNotationPattern

    enzyme
      Enzyme name, can be provided instead of a pattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``

    """

    best_possible_score = 0
    priority = -1

    def __init__(self, pattern=None, enzyme=None,
                 location=None, occurences=1, boost=1.0):
        """Initialize."""
        if enzyme is not None:
            pattern = enzyme_pattern(enzyme)
        if isinstance(pattern, str):
            pattern = DnaNotationPattern(pattern)
        self.pattern = pattern
        self.location = location
        self.enzyme = enzyme
        self.boost = boost
        self.occurences = occurences
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        if self.location is None:
            location = Location(0, len(problem.sequence))
            return self.copy_with_changes(location=location)
        else:
            return self

    def evaluate(self, problem):
        """Score the difference between expected and observed n_occurences."""
        matches = self.pattern.find_matches(problem.sequence, self.location)
        score = -abs(len(matches) - self.occurences)

        if score == 0:
            message = "Passed. Pattern found at positions %s" % matches
        else:
            if self.occurences == 0:
                message = "Failed. Pattern not found."
            else:
                message = ("Failed. Pattern found %d times instead of %d"
                           " wanted, at locations %s") % (len(matches),
                                                          self.occurences,
                                                          matches)
        return SpecEvaluation(
            self, problem, score, message=message,
            locations=[self.location],
            data=dict(matches=matches)
        )

        def localized(self, location):
            """Localize the evaluation."""
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return VoidSpecification(parent_specification=self)
            else:
                return self


    def insert_pattern_in_problem(self, problem):
        """Insert the pattern in the problem's sequence by successive tries.

        This heuristic is attempted to get the number of occurences in the
        pattern from 0 to some number
        """
        L = self.pattern.size
        localized_constraints = [
            c#.localized(self.location)
            for c in problem.constraints
        ]
        for start in range(self.location.start, self.location.end - L):
            #print (start)
            new_location = Location(start, start + L, self.location.strand)
            new_constraint = EnforceSequence(sequence=self.pattern.sequence,
                                             location=new_location)
            new_space = MutationSpace.from_optimization_problem(
                problem, new_constraints=[new_constraint])
            #print ('new_space', new_space.choices)
            if len(new_space.unsolvable_segments) > 0:
                #print ("bummer unsolvable")
                continue
            new_sequence = new_space.constrain_sequence(problem.sequence)
            new_constraints = localized_constraints + [new_constraint]
            new_problem = DnaOptimizationProblem(
                sequence=new_sequence,
                constraints=new_constraints,
                mutation_space=new_space,
                progress_logger=None
            )
            assert self.evaluate(new_problem).passes
            try:
                new_problem.resolve_constraints()
                #print (new_problem.constraints_text_summary())
                #print (new_problem.mutation_space.choices)
                problem.sequence = new_problem.sequence
                return
            except NoSolutionError:
                # print ("bummer")
                pass
        raise NoSolutionError(problem=problem, location=self.location,
            message='Insertion of pattern %s in %s failed' % (
                self.pattern.sequence, self.location
            ))

    def resolution_heuristic(self, problem):
        evaluation = self.evaluate(problem)
        if evaluation.passes:
            return
        if len(evaluation.data["matches"]) == self.occurences - 1:
            self.insert_pattern_in_problem(problem)
        else:
            problem.resolve_constraints_locally()  # default resolution method


    def __repr__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.location)

    def __str__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.location)
