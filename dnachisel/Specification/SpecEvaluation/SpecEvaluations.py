from ...biotools import score_to_formatted_string
from ...Location import Location
from ...reports import colors_cycle


class SpecEvaluations:
    """Base class for handling lists of SpecEvaluations.

    See ProblemObjectivesEvaluations and ProblemConstraintsEvaluations for
    the useful subclasses.

    Parameters
    ----------
    evaluations
      list of SpecEvaluations

    problem
      (optional) problem on which the evaluations were carried.

    """

    color_lightness = 1.0
    color_shift = 0

    def __init__(self, evaluations, problem=None):
        """Initialize."""
        self.evaluations = evaluations
        self.problem = problem

    def __iter__(self):
        """Iterate over evaluations."""
        return self.evaluations.__iter__()

    def __len__(self):
        """Return the number of evaluations."""
        return len(self.evaluations)

    def all_evaluations_pass(self):
        """Return whether all evaluations pass."""
        return all([ev.passes for ev in self.evaluations])

    def scores_sum(self, as_text=False):
        """Return the sum of all evaluations scores, factoring their boost.

        Scores are multiplied by their respective boost factor.
        """
        result = sum(
            [ev.specification.boost * ev.score for ev in self.evaluations]
        )
        if as_text:
            result = score_to_formatted_string(result)
        return result

    def filter(self, eval_filter):
        """Create a new instance with a subset of the evaluations.

        ``eval_filter`` is either a function ``(evaluation) => True/False``
        or one of "passing", "failing", "optimal", "suboptimal", to obtain the
        corresponding constraints.
        """
        if isinstance(eval_filter, str):
            eval_filter = {
                "passing": lambda e: e.passes,
                "failing": lambda e: not e.passes,
                "optimal": lambda e: e.is_optimal,
                "suboptimal": lambda e: not e.is_optimal,
            }[eval_filter]
        return self.__class__(
            evaluations=[ev for ev in self.evaluations if eval_filter(ev)],
            problem=self.problem,
        )

    def all_locations(self):
        return [
            location
            for evaluation in self.evaluations_with_locations()
            for location in evaluation.locations
        ]

    def to_text(self):
        """Return a long representation of the evaluations list."""
        return (
            "\n".join(
                ["===> %s" % self.text_summary_message()]
                + [
                    e.to_text(role=self.specifications_role)
                    for e in self.evaluations
                ]
            )
            + "\n\n"
        )

    def evaluations_with_locations(self):
        """Return the list of all evaluations whose location is not None."""
        return [ev for ev in self.evaluations if ev.locations is not None]

    def success_and_failures_as_features(self, feature_type="misc_feature"):
        """Return all evaluations as Biopython features.

        With color property depending on whether the evaluation is passing,
        failing, optimal, or suboptimal, the color being determined by
        ``self.success_failure_color``.

        """
        return [
            evaluation.specification.to_biopython_feature(
                feature_type=feature_type,
                color=self.success_failure_color(evaluation),
                passes="true" if evaluation.passes else "false",
                is_optimal="true" if evaluation.is_optimal else "false",
            )
            for evaluation in self.evaluations
            if evaluation.specification.__dict__.get("location", False)
        ]

    def locations_as_features(
        self,
        features_type="misc_feature",
        with_specifications=True,
        label_prefix="From",
        colors="cycle",
        merge_overlapping=False,
        locations_filter=None,
    ):
        """Return all locations from all evaluations as biopython features.

        Parameters
        ----------
        features_type
          The Genbank feature type

        with_specifications
          If True, features are added to the list to indicate the scope of the
          different specifications. If false, only the specification breaches
          are returned.

        label_prefix
          Each breach may be labeled "prefix NameOfSpec(props)", for instance,
          "From AvoidPattern(100-200)", to indicate where the breach belongs.

        colors
          Either a list of colors (one for each specification), e.g.
          ['red', '#7aab71', ...] or "cycle" for cycling through predefined
          colors. The colors are applied to all breaches.

        locations_filter
          A function (location => True/False) deciding whether

        """
        if colors == "cycle":
            cycle = colors_cycle(
                lightness_factor=self.color_lightness,
                color_shift=self.color_shift,
            )
            colors = [next(cycle) for ev in self.evaluations]

        features = [
            location.to_biopython_feature(
                feature_type="misc_feature",
                label=" ".join(
                    [
                        label_prefix,
                        ev.specification.label(
                            use_short_form=True, with_location=False
                        ),
                    ]
                ),
                color=color,
                ApEinfo_fwdcolor=color,
                ApEinfo_revcolor=color,
            )
            for (ev, color) in zip(self.evaluations_with_locations(), colors)
            for location in (
                Location.merge_overlapping_locations(ev.locations)
                if merge_overlapping
                else ev.locations
            )
            if (locations_filter is None) or locations_filter(location)
        ]
        if with_specifications:
            features += [
                ev.specification.to_biopython_feature(
                    feature_type="misc_feature",
                    label=ev.specification.label(
                        use_short_form=True, with_location=False
                    ),
                    role=self.specifications_role,
                    color=color,
                    ApEinfo_fwdcolor=color,
                    ApEinfo_revcolor=color,
                )
                for ev, color in zip(self.evaluations, colors)
                if ev.specification.__dict__.get("location", False)
            ]
        return features
