try:
    from dna_features_viewer import BiopythonTranslator
    DSV_AVAILABLE = True
except:
    BiopythonTranslator = object
    DSV_AVAILABLE = False

class GraphicTranslator(BiopythonTranslator):
    @staticmethod
    def compute_feature_color(f):
        default = BiopythonTranslator.compute_feature_color(f)
        if f.type == "misc_feature":
            label = BiopythonTranslator.compute_feature_label(f)
            return {"@": "#ce5454", "~": "#e5be54"}.get(label[0], "#f4df42")
        else:
            return "#d9e0fc"
    @staticmethod
    def compute_feature_label(f):
        default = BiopythonTranslator.compute_feature_label(f)
        return None if (f.type != "misc_feature") else default
