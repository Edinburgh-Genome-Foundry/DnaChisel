try:

    from dna_features_viewer import BiopythonTranslator

    DFV_AVAILABLE = True
except ImportError:
    BiopythonTranslator = object
    DFV_AVAILABLE = False


class GraphicTranslator(BiopythonTranslator):
    """A Biopython record translator for DNA Features Viewer.

        This translator produces label-free plots.
        """

    @staticmethod
    def compute_feature_box_linewidth(f):
        return 1 if f.qualifiers.get("is_a_breach", False) else 0
    
    @staticmethod
    def compute_feature_fontdict(f):
        return {
            "fontsize": 12 if f.qualifiers.get("is_a_breach", False) else 9
        }

    def compute_feature_label(self, f):
        label = BiopythonTranslator.compute_feature_label(self, f)
        if not f.qualifiers.get("is_a_breach", False):
            if len(label) > 20:
                label = label[:19] + "…"
        return label

    def compute_feature_color(self, f):
        if f.qualifiers.get("is_a_breach", False):
            return f.qualifiers.get("color", "red")
        else:
            return "#ffffff"
