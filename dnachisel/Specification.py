"""Base class for specification.

Notable features implemented here:

- Many empty methods that features subclass will overwrite
- Feature import/export from/to Genbank features.
"""
import copy
import re

from .biotools import find_specification_in_feature
from .Location import Location
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
    optimize_passively = False
    enforced_by_nucleotide_restrictions = False
    priority = 0
    is_void = False

    def __init__(self, evaluate=None, boost=1.0):
        """Initialize."""
        self.boost = boost
        if evaluate is not None:
            self.evaluate = evaluate

    def localized(self, location, problem=None):
        """Return a modified version of the specification for the case where
        sequence modifications are only performed inside the provided location.

        For instance if an specification concerns local GC content, and we are
        only making local mutations to destroy a restriction site, then we only
        need to check the local GC content around the restriction site after
        each mutation (and not compute it for the whole sequence), so
        ``EnforceGCContent.localized(location)`` will return an specification
        that only looks for GC content around the provided location.

        If an specification concerns a DNA segment that is completely disjoint
        from the provided location, this must return None.

        Must return an object of class ``Constraint``.
        """
        return self

    def copy_with_changes(self, **kwargs):
        """Return a copy of the Specification with modified properties.

        For instance ``new_spec = spec.copy_with_changes(boost=10)``.
        """
        new_specification = copy.copy(self)
        new_specification.__dict__.update(kwargs)
        return new_specification

    def shifted(self, shift):
        """Shift the location of the specification.

        Some specification classes may have a special method to do side effects
        when shifting the location.

        Location shifting is used in particular when solving circular DNA
        optimization problems.
        """
        new_location = None if self.location is None else self.location + shift
        return self.copy_with_changes(location=new_location, derived_from=self)

    def initialized_on_problem(self, problem, role="constraint"):
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

        # PARSE THE SPECIFICATION, IDENTIFY THE TYPE AND ARGUMENTS

        label = find_specification_in_feature(feature)
        if isinstance(label, list):
            label = label[0]
        if not label.endswith(")"):
            # Standardizes the expression: @cds => @cds()
            label += "()"
        # The regular expression below detects spec definitions:
        # ~Avoidpattern(ARGS) => ~, AvoidPattern, ARGS
        pattern = r"([@~])(\S+)(\(.*\))"
        match = re.match(pattern, label)
        role, specification, parameters = match.groups()
        if specification not in specifications_dict:
            raise TypeError("Unknown specification %s" % specification)
        specification_class = specifications_dict[specification]
        role = {"@": "constraint", "~": "objective"}[role]

        # PARSE THE ARGUMENTS AND KEYWORD ARGUMENTS

        def format_value(value):
            """Converts stringed integers and floats back to numerical.
            Also converts "'bla'" => "bla"
            If the value is a list, apply to all elements."""
            if isinstance(value, (list, tuple)):
                return [format_value(v) for v in value]
            match = re.match(r"'(.*)'", value)
            if match is not None:
                return match.groups()[0]
            else:
                try:
                    return int(value)
                except ValueError:
                    try:
                        return float(value)
                    except Exception:
                        return value

        args, kwargs = [], {}
        for arg in parameters[1:-1].split(", "):
            if arg == "":
                continue
            if ":" in arg:
                key, value = arg.split(":")
                if "|" in value:
                    value = value.split("|")
                kwargs[key] = format_value(value)
            elif "=" in arg:
                key, value = arg.split("=")
                if "|" in value:
                    value = value.split("|")
                kwargs[key] = format_value(value)
            else:
                args.append(format_value(arg))
        kwargs["location"] = Location.from_biopython_location(feature.location)
        # ATTEMPT TO CREATE A SPECIFICATION WITH THE GIVEN TYPE AND ARGS

        try:
            specification_instance = specification_class(*args, **kwargs)
        except TypeError as err:
            message = err.args[0]
            faulty_parameter = message.split("'")[1]
            raise TypeError(
                "Unknown parameter %s for specification %s "
                "at location %s"
                % (faulty_parameter, specification, kwargs["location"])
            )

        return role, specification_instance

    def label(self, role=None, with_location=True, assignment=":"):
        prefix = {"constraint": "@", "objective": "~", None: ""}[role]
        if with_location and hasattr(self, "location") and self.location:
            location = "[%s]" % self.location
        else:
            location = ""
        params = self.label_parameters()
        if params == []:
            params = ""
        else:
            params = "(%s)" % ", ".join(
                [
                    assignment.join(map(str, p)) if isinstance(p, tuple) else p
                    for p in params
                ]
            )
        return "".join([prefix, self.__class__.__name__, location, params])

    def label_parameters(self):
        """In subclasses, returns a list of the creation parameters.

        For instance [('pattern', 'ATT'), ('occurences', 2)]
        """
        return []

    def __str__(self):
        """By default, represent the Specification using its label()"""
        return self.label()

    def __repr__(self):
        """By default, represent the Specification using its label()"""
        return self.label()

    def short_label(self):
        """Shorter, less precise label to be used in tables, reports, etc.

        This is meant for specifications such as EnforceGCContent(0.4, 0.6)
        to be represented as '40-60% GC' in reports tables etc..
        """
        return self.label()

    def to_biopython_feature(
        self,
        feature_type="misc_feature",
        role="constraint",
        colors_dict=None,
        **qualifiers
    ):
        """Return a Biopython feature representing the specification.

        The feature label is a string representation of the specification,
        and its location indicates the specification's scope.

        This method is primarily meant to display specifications in Genbank.
        They may result in "viable", DnaChisel-compatible annotations that can
        be imported back into DNA Chisel, but this is not the intended goal.

        """
        if colors_dict is None:
            colors_dict = {"constraint": "#355c87", "objective": "#f9cd60"}
        qualifiers["role"] = role
        if "label" not in qualifiers:
            qualifiers["label"] = self.label(
                role=role, with_location=False, assignment=":"
            )
        if "color" not in qualifiers:
            qualifiers["color"] = colors_dict[role]
        return SeqFeature(
            self.location.to_biopython_location(),
            type=feature_type,
            qualifiers=qualifiers,
        )

    def restrict_nucleotides(self, sequence, location=None):
        """Restrict the mutation space to speed up optimization.

        This method only kicks in when this specification is used as a
        constraint. By default it does nothing, but subclasses such as
        EnforceTranslation, AvoidChanges, EnforceSequence, etc. have custom
        methods.

        In the code, this method is run during the initialize() step of
        DNAOptimizationProblem, when the MutationSpace is created for each constraint
        """
        return []

    def as_passive_objective(self):
        """Return a copy with optimize_passively set to true.

        "Optimize passively" means that when the specification is used as an
        objective, the solver will not do a specific pass to optimize this
        specification, however this specification's score will be taken into
        account in the global score when optimizing other objectives, and
        may therefore influence the final sequence.
        """
        return self.copy_with_changes(optimize_passively=True)

    def _copy_with_full_span_if_no_location(self, problem):
        """Return either self, or a copy with location "everywhere".

        And by "everywhere" we mean Location(0, L) where L is the problem's
        sequence length.

        Most Specifications use this method in their "initialized_on_problem()"
        custom method.
        """
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            return self.copy_with_changes(location=location)
        else:
            return self


class SpecificationsSet:
    """Generic class for writing Specs which are actually made of more specs.

    Behaves as a Specification when it comes to instanciation, reading it
    from annotated records, etc. but the initialization actually creates a
    dictionnary of standard Specifications in the DNAOptimizationProblem
    """

    def register_specifications(self, specifications):
        for name, spec in specifications.items():
            spec.parent_specification = self
            spec.name_in_parent = name
        self.specifications = specifications
