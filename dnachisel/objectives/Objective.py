import copy
import re

from ..biotools import (DnaNotationPattern, enzyme_pattern,
                        find_objective_in_feature)
from ..Location import Location
from .ObjectiveEvaluation import ObjectiveEvaluation
from Bio.SeqFeature import SeqFeature


class Objective:
    """General class to define objectives to optimize.

    Note that all objective have a ``boost`` attribute that is a multiplicator
    that will be used when computing the global objective score of a problem
    with ``problem.all_objectives_score()``.

    New types of objectives are defined by subclassing ``Objective`` and
    providing a custom ``evaluate`` and ``localized`` methods.

    """

    best_possible_score = None
    can_be_solved_locally = False

    def __init__(self, evaluate=None, boost=1.0):
        self.boost = boost
        if evaluate is not None:
            self.evaluate = evaluate

    def localized(self, location):
        """Return a modified version of the objective for the case where
        sequence modifications are only performed inside the provided location.

        For instance if an objective concerns local GC content, and we are
        only making local mutations to destroy a restriction site, then we only
        need to check the local GC content around the restriction site after
        each mutation (and not compute it for the whole sequence), so
        ``EnforceGCContent.localized(location)`` will return an objective
        that only looks for GC content around the provided location.

        If an objective concerns a DNA segment that is completely disjoint from
        the provided location, this must return a ``VoidConstraint``.

        Must return an object of class ``Constraint``.
        """
        return self

    def copy_with_changes(self, **kwargs):
        new_objective = copy.deepcopy(self)
        new_objective.__dict__.update(kwargs)
        return new_objective

    def initialize_problem(self, problem, role="constraint"):
        return self

    @staticmethod
    def from_biopython_feature(feature, objectives_dict):
        label = find_objective_in_feature(feature)
        if isinstance(label, list):
            label = label[0]
        if not label.endswith(")"):
            label += "()"
        pattern = "([@~])(\S+)(\(.*\))"
        match = re.match(pattern, label)
        role, objective, parameters = match.groups()
        role = {"@": "constraint", "~": "objective"}[role]
        kwargs = dict(e.split('=') for e in parameters[1:-1].split(', ')
                    if ("=" in e))
        for k, v in kwargs.items():
            match = re.match(r"'(.*)'", v)
            if match is not None:
                kwargs[k] = match.groups()[0]
            else:
                try:
                    kwargs[k] = int(v)
                except ValueError:
                    try:
                        kwargs[k] = float(v)
                    except:
                        pass
        kwargs["location"] = Location.from_biopython_location(feature.location)
        return role, objectives_dict[objective](**kwargs)

    def to_biopython_feature(self, feature_type="misc_feature",
                             role="constraint", colors_dict=None,
                             **qualifiers):
        if colors_dict is None:
            colors_dict = {"constraint": "#355c87", "objective": "#f9cd60"}
        qualifiers["role"] = role
        if "label" not in qualifiers:
            qualifiers['label'] = self.__repr__()

        if "color" not in qualifiers:
            qualifiers['color'] = colors_dict[role]
        return SeqFeature(self.location.to_biopython_location(),
                          type=feature_type,
                          qualifiers=qualifiers)



class VoidObjective(Objective):
    """Void Objectives are a special case of Objectives that always pass.

    Void Objectives are generally obtained when a Objective is "made void"
    by a localization. For instance, if we are optimizing the segment (10,50)
    of a DNA segment, the Objective EnforceTranslation([300,500]) does not
    apply as it concerns a Gene that is in a completely different segment.
    Therefore the localized version of EnforceTranslation will be void.

    Note: the initializer accepts starred arguments/keyword arguments to make
    it easy to void any other objective by replacing the class to Void.
    Particularly useful when importing an optimization problem from genbank.
    """
    best_possible_score = 0



    def __init__(self, parent_objective=None, boost=1.0, *a, **kw):
        self.parent_objective = parent_objective
        self.message = ("Pass (not relevant in this context)")
        self.boost = boost

    def evaluate(self, problem):
        """The evaluation of VoidObjectives always passes with score=1.0
        It returns a message indicating that the parent Objective was voided
        """
        return ObjectiveEvaluation(self, problem, score=1.0,
                                   message=self.message,
                                   locations=None)

    def __repr__(self):
        return "Voided %s" % repr(self.parent_objective)

class PatternObjective(Objective):
    """Class for Objectives such as presence or absence of a pattern.

    The particularity of the PatternObjectives is that they will either infer
    or ask for the length of the associated pattern and use this to localize
    the objective efficiently when performing local optimization or solving.

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    dna_pattern
      A string of ATGC that will be converted automatically to a DNA pattern

    enzyme
      Enzyme name, can be provided instead of pattern or dna_pattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``


    """
    shrink_when_localized = True

    def __init__(self, pattern=None, location=None, boost=1.0, enzyme=None,
                 dna_pattern=None):
        if enzyme is not None:
            pattern = enzyme_pattern(enzyme)
        if dna_pattern is not None:
            pattern = DnaNotationPattern(dna_pattern)
        self.pattern = pattern
        self.location = location
        self.dna_pattern = dna_pattern
        self.enzyme = enzyme
        self.boost = boost

    def localized(self, location):
        """Localize the pattern to the given location. Taking into account the
        objective's own location, and the size of the pattern."""
        pattern_size = self.pattern.size
        if self.location is None:
            new_location = location.extended(pattern_size - 1)
        else:
            if self.location.overlap_region(location) is None:
                return VoidObjective(parent_objective=self)
            else:
                if not self.shrink_when_localized:
                    return self
                extended_location = location.extended(pattern_size - 1)
                new_location = self.location.overlap_region(extended_location)

        return self.copy_with_changes(location=new_location)

class TerminalObjective(Objective):
    """Objectives that apply in the same way to both ends of the sequence.

    These are particularly useful for modeling constraints from providers
    who have terminal-ends constraints.

    Subclasses of these objectives should have a `location_size` and a
    `evaluate_end` method"""

    def evaluate(self, problem):
        sequence = problem.sequence
        L = len(sequence)
        wsize = self.window_size
        locations = [
            location
            for location in [Location(0, wsize), Location(L - wsize, L)]
            if not self.evaluate_end(location.extract_sequence(sequence))
        ]

        if locations == []:
            message = "Passed (no breach at the ends)"
        else:
            message = "Failed: breaches at ends %s" % str(locations)

        return ObjectiveEvaluation(self, problem, score=len(locations),
                                   locations=locations, message=message)
