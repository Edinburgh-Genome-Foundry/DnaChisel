from ...biotools import (
    load_record,
    write_record,
    sequence_to_biopython_record,
    find_specification_label_in_feature,
    sequences_differences_segments,
)
from ...Specification.Specification import Specification
from ...Location import Location


class RecordRepresentationMixin:
    """Mixin for DnaOptimizationProblem gathering methods for converting
    DnaOptimizationProblems from/to Biopython records and Genbank records."""

    @classmethod
    def from_record(
        cls,
        record,
        specifications_dict="default",
        logger="bar",
        extra_constraints=(),
        extra_objectives=(),
    ):
        """Create a DnaOptimizationProblem by parsing a record's annotations.

        Examples
        --------

        >>> problem = DnaOptimizationProblem.from_record("my_record.gb")
        >>> problem = DnaOptimizationProblem.from_record(some_biopython_record)

        Parameters
        ----------

        record
          Either a biopython record or path to a genbank/snapgene file.

        specifications_dict
          Provide a custom dict with user-defined specifications instead of the
          default dict, which contains the built-in specifications.

        logger
          Logger of the DnaOptimizationProblem

        extra_constraints, extra_objectives
          List of Specifications to be added to the problem, in addition to
          the specifications parsed from the genbank.

        """
        # unfortunately the local import below is the most elegant found so
        # far. builtin_specifications cannot be imported at the top of this
        # file as some built-in specifications use DnaOptimizationProblem
        # internally to resolve constructs (see EnforcePatternOccurences)
        if isinstance(record, str):
            record = load_record(record)
        parameters = dict(
            sequence=record,
            constraints=[] + list(extra_constraints),  # shallow copy
            objectives=[] + list(extra_objectives),  # shallow copy
            logger=logger,
        )
        for feature in record.features:
            if feature.type != "misc_feature":
                continue
            label = find_specification_label_in_feature(feature)
            if label is None:
                continue
            specs = Specification.list_from_biopython_feature(
                feature, specifications_dict=specifications_dict
            )
            for role, specification in specs:
                parameters[role + "s"].append(specification)
        return cls(**parameters)

    def to_record(
        self,
        filepath=None,
        features_type="misc_feature",
        with_original_features=True,
        with_original_spec_features=False,
        with_constraints=True,
        with_objectives=True,
        with_sequence_edits=False,
        colors_dict=None,
        use_short_labels=True,
        record_id = None
    ):
        """Return/write record representing the final sequence and problem.

        the many options enable to also annotate specifications, sequence
        edits, etc.

        Parameters
        ----------
        filepath
          Path to a target genbank file where the record will be written. If
          none is provided, a Biopython record is returned instead

        features_type
          Genbank standard type to give to all extra genbank annotation created
          by this method (to indicate constraints, objectives, edits, etc.)


        with_original_features
          If True, the features from the original record provided at problem
          creation will be included in the record (if a simple sequence, and
          not a record, was originally provided, then there is no such
          features)


        with_original_spec_features
          If False, any feature from the original record provided at problem
          creation that defines a DNAChisel Specification be stripped off the
          record returned by this method (to make space for the annotations)
          created by this method

        with_constraints
          If True, annotations representing the constraints will be added to
          the record


        with_objectives
          If True, annotations representing the objectives will be added to the
          record


        with_sequence_edits
          If True, annotations representing each nucleotide change will be
          added to the record.


        colors_dict
          A dict indicating the feature color for constraints and objectives.
          The default is {"constraint": "#355c87", "objective": "#f9cd60"}.



        use_short_labels
          If True, the annotations representing constraints and objectives will
          use shorter labels to indicate the type of specification.


        Notes
        -----

        If the original problem was created from a Genbank, it is a good idea
        to set with_original_spec_features=True and with_constraints=False,
        with_objectives=False.


        """
        record = sequence_to_biopython_record(self.sequence)
        if record_id is not None:
            record.id = record_id

        record.features = []
        if with_constraints:
            record.features += [
                cst.to_biopython_feature(
                    role="constraint",
                    feature_type=features_type,
                    colors_dict=colors_dict,
                    use_short_label=use_short_labels,
                )
                for cst in self.constraints
                if cst.__dict__.get("location", False)
            ]
        if with_objectives:
            record.features += [
                obj.to_biopython_feature(
                    role="objective",
                    feature_type=features_type,
                    colors_dict=colors_dict,
                    use_short_label=use_short_labels,
                )
                for obj in self.objectives
            ]
        if with_original_features and (self.record is not None):
            record.features += [
                f
                for f in self.record.features
                if with_original_spec_features
                or not find_specification_label_in_feature(f)
            ]
        if with_sequence_edits:
            record.features += self.sequence_edits_as_features()

        if filepath is not None:
            write_record(record=record, target=filepath, file_format="genbank")
        else:
            return record

    def sequence_edits_as_features(self, feature_type="misc_feature"):
        """Return a list of Biopython Record Features indicating each of the
        edits."""
        segments = sequences_differences_segments(
            self.sequence, self.sequence_before
        )
        return [
            Location(start, end).to_biopython_feature(
                label="%s=>%s"
                % (self.sequence_before[start:end], self.sequence[start:end]),
                is_edit="true",
                ApEinfo_fwdcolor="#ff0000",
                color="#ff0000",
            )
            for start, end in segments
        ]
