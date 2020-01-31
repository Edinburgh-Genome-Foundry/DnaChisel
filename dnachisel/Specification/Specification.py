"""Base class for specification.

Notable features implemented here:

- Many empty methods that features subclass will overwrite
- Feature import/export from/to Genbank features.
"""
import copy
from ..Location import Location
from .FeatureRepresentationMixin import FeatureRepresentationMixin


class Specification(FeatureRepresentationMixin):
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

    Attributes
    ----------

    best_possible_score
      Best score that the specification can achieve. Used by the optimization
      algorithm to understand when no more optimization is required.

    optimize_passively (boolean)
      Indicates that there should not be a pass of the optimization algorithm
      to optimize this objective. Instead, this objective is simply taken into
      account when optimizing other objectives.

    enforced_by_nucleotide_restrictions (boolean)
      When the spec is used as a constraints, this indicates that the
      constraint will initially restrict the mutation space in a way that
      ensures that the constraint will always be valid. The constraint does
      not need to be evaluated again, which speeds up the resolution algorithm.

    priority
      Value used to sort the specifications and solve/optimize them in order,
      with highest priority first.

    shorthand_name
      Shorter name for the specification class that will be recognized when
      parsing annotations from genbanks.

    """

    best_possible_score = None
    optimize_passively = False
    enforced_by_nucleotide_restrictions = False
    priority = 0
    shorthand_name = None
    is_focus = False

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

    def label(
        self,
        role=None,
        with_location=True,
        assignment_symbol=":",
        use_short_form=False,
        use_breach_form=False
    ):
        """Return a string label for this specification.

        Parameters
        ----------

        role
          Either 'constraint' or 'objective' (for prefixing the label with @
          or ~), or None.

        with_location
          If true, the location will appear in the label.

        assignment_symbol
          Indicates whether to use ":" or "=" or anything else when indicating
          parameters values.

        use_short_form
          If True, the label will use self.short_label(), so for instance
          AvoidPattern(BsmBI_site) will become "no BsmBI". How this is handled
          is dependent on the specification.
        """
        prefix = {"constraint": "@", "objective": "~", None: ""}[role]
        if use_short_form:
            label = self.short_label()
            if with_location:
                label += ", %s" % self.location
            return prefix + label
        if use_breach_form:
            label = self.breach_label()
            if with_location:
                label += ", %s" % self.location
            return label
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
                    assignment_symbol.join(map(str, p))
                    if isinstance(p, tuple)
                    else str(p)
                    for p in params
                ]
            )

        return "".join([prefix, self.__class__.__name__, location, params])

    def short_label(self):
        """Shorter, less precise label to be used in tables, reports, etc.

        This is meant for specifications such as EnforceGCContent(0.4, 0.6)
        to be represented as '40-60% GC' in reports tables etc..
        """
        return self.__class__.__name__
    
    def breach_label(self):
        """Shorter, less precise label to be used in tables, reports, etc.

        This is meant for specifications such as EnforceGCContent(0.4, 0.6)
        to be represented as '40-60% GC' in reports tables etc..
        """
        return "'%s' breach" % (self.short_label())

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

    def restrict_nucleotides(self, sequence, location=None):
        """Restrict the mutation space to speed up optimization.

        This method only kicks in when this specification is used as a
        constraint. By default it does nothing, but subclasses such as
        EnforceTranslation, AvoidChanges, EnforceSequence, etc. have custom
        methods.

        In the code, this method is run during the initialize() step of
        DNAOptimizationProblem, when the MutationSpace is created for each
        constraint
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
            location = Location(0, len(problem.sequence))
            return self.copy_with_changes(location=location)
        else:
            return self
