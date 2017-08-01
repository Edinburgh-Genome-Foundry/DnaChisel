import copy
import re

from ..biotools import (DnaNotationPattern, enzyme_pattern,
                        find_specification_in_feature)
from ..Location import Location
from .SpecEvaluation import SpecEvaluation
from Bio.SeqFeature import SeqFeature


class Specification:
    """General class to define specifications to optimize.

    Note that all specifications have a ``boost`` attribute that is a
    multiplicator that will be used when computing the global specification
    score of a problem with ``problem.all_objectives_score()``.

    New types of specifications are defined by subclassing ``Specification`` and
    providing a custom ``evaluate`` and ``localized`` methods.

    Parameters
    -----------
    evaluate
      function (sequence) => SpecEvaluation
    boost
      Relative importance of the Specification's score in a multi-specification
      problem.
    """

    best_possible_score = None
    can_be_solved_locally = False

    def __init__(self, evaluate=None, boost=1.0):
        """Initialize."""
        self.boost = boost
        if evaluate is not None:
            self.evaluate = evaluate

    def localized(self, location):
        """Return a modified version of the specification for the case where
        sequence modifications are only performed inside the provided location.

        For instance if an specification concerns local GC content, and we are
        only making local mutations to destroy a restriction site, then we only
        need to check the local GC content around the restriction site after
        each mutation (and not compute it for the whole sequence), so
        ``EnforceGCContent.localized(location)`` will return an specification
        that only looks for GC content around the provided location.

        If an specification concerns a DNA segment that is completely disjoint
        from the provided location, this must return a ``VoidSpecification``.

        Must return an object of class ``Constraint``.
        """
        return self

    def copy_with_changes(self, **kwargs):
        """Return a copy of the Specification with modified properties.

        For instance ``new_spec = spec.copy_with_changes(boost=10)``.
        """
        new_specification = copy.deepcopy(self)
        new_specification.__dict__.update(kwargs)
        return new_specification

    def initialize_on_problem(self, problem, role="constraint"):
        """Complete specification initialization when the sequence gets known.

        Some specifications like to know what their role is and on which
        sequence they are employed before they complete some values.
        """
        return self

    @staticmethod
    def from_biopython_feature(feature, specifications_dict):
        """Parse a Biopython feature create an annotation.

        The specifications_dict enables to map specification names to the
        actual implemented class.

        """

        label = find_specification_in_feature(feature)
        if isinstance(label, list):
            label = label[0]
        if not label.endswith(")"):
            label += "()"
        pattern = "([@~])(\S+)(\(.*\))"
        match = re.match(pattern, label)
        role, specification, parameters = match.groups()
        role = {"@": "constraint", "~": "specification"}[role]
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
        return role, specifications_dict[specification](**kwargs)

    def to_biopython_feature(self, feature_type="misc_feature",
                             role="constraint", colors_dict=None,
                             **qualifiers):
        """Return a Biopython feature representing the specification.

        The feature label is a string representation of the specification,
        and its location indicates the specification's scope.

        """
        if colors_dict is None:
            colors_dict = {"constraint": "#355c87", "specification": "#f9cd60"}
        qualifiers["role"] = role
        if "label" not in qualifiers:
            qualifiers['label'] = self.__repr__()

        if "color" not in qualifiers:
            qualifiers['color'] = colors_dict[role]
        return SeqFeature(self.location.to_biopython_location(),
                          type=feature_type,
                          qualifiers=qualifiers)

    def restrict_nucleotides(self, sequence, location=None):
        """Restrict the possible nucleotides mutations in a sequence to speed
        up optimization.

        This method has no effect unless a special heuristic is implemented
        for it.
        """
        return []



class VoidSpecification(Specification):
    """Void Specifications are a special case of Specifications that always pass.

    Void Specifications are generally obtained when a Specification is "made void"
    by a localization. For instance, if we are optimizing the segment (10,50)
    of a DNA segment, the Specification EnforceTranslation([300,500]) does not
    apply as it concerns a Gene that is in a completely different segment.
    Therefore the localized version of EnforceTranslation will be void.

    Note: the initializer accepts starred arguments/keyword arguments to make
    it easy to void any other specification by replacing the class to Void.
    Particularly useful when importing an optimization problem from genbank.
    """
    best_possible_score = 0



    def __init__(self, parent_specification=None, boost=1.0, *a, **kw):
        """Initialize."""
        self.parent_specification = parent_specification
        self.message = ("Pass (not relevant in this context)")
        self.boost = boost

    def evaluate(self, problem):
        """The evaluation of VoidSpecifications always passes with score=1.0
        It returns a message indicating that the parent Specification was voided
        """
        return SpecEvaluation(self, problem, score=1.0,
                                   message=self.message,
                                   locations=None)

    def __repr__(self):
        return "Voided %s" % repr(self.parent_specification)

class PatternSpecification(Specification):
    """Class for Specifications such as presence or absence of a pattern.

    The particularity of the PatternSpecifications is that they will either
    infer or ask for the length of the associated pattern and use this to
    localize the specification efficiently when performing local optimization or
    solving.

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
        """Initialize."""
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
        specification's own location, and the size of the pattern."""
        pattern_size = self.pattern.size
        if self.location is None:
            new_location = location.extended(pattern_size - 1)
        else:
            if self.location.overlap_region(location) is None:
                return VoidSpecification(parent_specification=self)
            else:
                if not self.shrink_when_localized:
                    return self
                extended_location = location.extended(pattern_size - 1)
                new_location = self.location.overlap_region(extended_location)

        return self.copy_with_changes(location=new_location)

class TerminalSpecification(Specification):
    """Specifications that apply in the same way to both ends of the sequence.

    These are particularly useful for modeling constraints from providers
    who have terminal-ends constraints.

    Subclasses of these specifications should have a `location_size` and a
    `evaluate_end` method"""

    def evaluate(self, problem):
        """Apply method ``evaluate_end`` to both sides and compile results."""
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

        return SpecEvaluation(self, problem, score=-len(locations),
                                   locations=locations, message=message)
