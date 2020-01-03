from Bio.SeqRecord import SeqRecord
from ..Specification.SpecificationSet import SpecificationSet
from ..DnaOptimizationProblem.mixins import ConstraintsSolverMixin
from ..MutationSpace import MutationSpace


class DnaDesignProblem(
    ConstraintsSolverMixin,
    
    ):

    def __init__(
        self,
        sequence,
        constraints=None,
        mutation_space=None,
        # tailor args
        design_space = None,
        solution_id = None,
        parent=None):

        if isinstance(sequence, SeqRecord):
            self.record = sequence
            self.sequence = str(sequence.seq).upper()
        else:
            self.record = None
            self.sequence = sequence.upper()
        self.constraints = [] if constraints is None else list(constraints)
        self.mutation_space = mutation_space

        self.design_space = design_space
        self.solution_id = solution_id
        self.parent = parent

        self.scores = {}
        self.levels = {}

        self.initialize()
    
    def _set_scores(self):
        for i in range(len(self.design_space.feature_label)):
            self.scores[self.design_space.feature_label[i]] = self.design_space.feature_specification[i].evaluate(self).score
    
    def _set_level(self):
        for i in range(len(self.design_space.feature_label)):
            feature_label = self.design_space.feature_label[i]
            self.levels[feature_label + "_Level"] = self.design_space.range_set_list[i].get_level(self.scores[feature_label])

            
            
    def initialize(self):

        # for specs in self.constraints:
        specsets = [
            spec for spec in self.constraints if isinstance(spec, SpecificationSet)
        ]
        specs_in_sets = [
            spec
            for specset in specsets
            for spec in specset.specifications.values()
        ]
        for specset in specsets:
            self.constraints.remove(specset)
        self.constraints.extend(specs_in_sets)

        # INITIALIZE THE CONSTRAINTS AND OBJECTIVES

        self.constraints = [
            constraint.initialized_on_problem(self, role="constraint")
            for constraint in self.constraints
        ]

        self.sequence_before = self.sequence
        self._constraints_before = None
        self._objectives_before = None

        # INITIALIZE THE MUTATION SPACE

        if self.mutation_space is None:
            self.mutation_space = MutationSpace.from_optimization_problem(self)
            # If the original sequence is outside of the allowed mutations
            # space, replace the sequence by a sequence which complies with
            # the mutation space.
            self.sequence = self.mutation_space.constrain_sequence(
                self.sequence
            )
        
        self._set_scores()
        self._set_level()


    def is_match_design(self,desired_design):
        levels = [
            str(self.levels[feature+'_Level']) 
            for feature in self.design_space.feature_label
        ]

        return '.'.join(levels) == desired_design