from ..biotools import find_specification_label_in_feature
from .tools import install_extras_message

DFV_AVAILABLE = False
try:
    from dna_features_viewer import BiopythonTranslator

    DFV_AVAILABLE = True
except ImportError:

    class BiopythonTranslator:
        "Class unavailable. Install DNA Features Viewer."

        def __init__(self):
            raise ImportError(install_extras_message("DNA Features Viewer"))


class SpecAnnotationsTranslator(BiopythonTranslator):
    """Translator of DnaChisel feature-constraints for DNA Features Viewer"""

    feature_prefixes_colors = {
        "@": "#ce5454",
        "~": "#e5be54",
        "#": "#8edfff",
        "!": "#fcff75",
    }

    def compute_filtered_features(self, features):
        """Do not display edits."""
        return [
            feature
            for feature in features
            if "".join(feature.qualifiers.get("is_edit", "false")) != "true"
        ]

    def compute_feature_color(self, f):
        color = f.qualifiers.get("color", None)
        if color is not None:
            if isinstance(color, list):
                color = color[0]
            return color

        if f.type == "misc_feature":
            specification = find_specification_label_in_feature(f)
            if specification is not None:
                return self.feature_prefixes_colors.get(
                    specification[0], "#f4df42"
                )
        return "#eeeafa"

    def compute_feature_label(self, f):
        is_edit = f.qualifiers.get("is_edit", "false")
        if "true" in [is_edit, is_edit[0]]:
            return None
        default = BiopythonTranslator.compute_feature_label(self, f)
        label = None if (f.type != "misc_feature") else default
        if label == "misc_feature":
            label = None
        return label
