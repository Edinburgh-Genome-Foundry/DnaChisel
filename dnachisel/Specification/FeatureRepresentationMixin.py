
import re
from Bio.SeqFeature import SeqFeature
from ..Location import Location
from ..biotools import find_specification_in_feature


class FeatureRepresentationMixin:
    """Mixin for class Specification. Defines methods for converting
    Specifications from and to Biopython records."""

    def to_biopython_feature(
        self,
        feature_type="misc_feature",
        role="constraint",
        colors_dict=None,
        use_short_label=True,
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
                role=role,
                with_location=False,
                assignment_symbol=":",
                use_short_form=use_short_label,
            )

        if "color" not in qualifiers:
            qualifiers['color'] = colors_dict[role]
        qualifiers.update(
            dict(
                ApEinfo_fwdcolor=qualifiers['color'],
                ApEinfo_revcolor=qualifiers['color'],
            )
        )
        return SeqFeature(
            self.location.to_biopython_location(),
            type=feature_type,
            qualifiers=qualifiers,
        )

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
