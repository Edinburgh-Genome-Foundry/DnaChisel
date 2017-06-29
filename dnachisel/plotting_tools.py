import itertools
from .biotools import find_objective_in_feature

MATPLOTLIB_AVAILABLE = False
DFV_AVAILABLE = False

import numpy as np
try:
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
    from dna_features_viewer import BiopythonTranslator
    DFV_AVAILABLE = True
except:
    class BiopythonTranslator:
        "Class unavailable. Install DNA Features Viewer."
        def __init__(self):
            raise ImportError("BiopythonTranslator unavailable. Install "
                              "DNA Features Viewer.")


def colors_cycle():
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("colors_cycle requires Matplotlib installed")
    cycle = itertools.cycle([cm.Paired(0.21 * i % 1.0) for i in range(30)])
    return (
        '#%02x%02x%02x' % tuple([int(255 * c) for c in rgb_tuple[:3]])
        for rgb_tuple in cycle
    )


class ObjectivesAnnotationsTranslator(BiopythonTranslator):
    """Translator of DnaChisel feature-constraints for DNA Features Viewer"""
    @staticmethod
    def compute_feature_color(f):
        color = f.qualifiers.get('color', None)
        if color is not None:
            if isinstance(color, list):
                color = color[0]
            return color

        if f.type == "misc_feature":
            objective = find_objective_in_feature(f)
            if objective is None:
                return "#f4df42"
            else:
                return {
                    "@": "#ce5454",
                    "~": "#e5be54",
                    "#": "#8edfff",
                    "!": "#fcff75",
                }.get(objective[0], "#f4df42")
        else:
            return "#eeeafa"

    @staticmethod
    def compute_feature_label(f):
        is_edit = f.qualifiers.get("is_edit", "false")
        if "true" in [is_edit, is_edit[0]]:
            return None
        default = BiopythonTranslator.compute_feature_label(f)
        label = None if (f.type != "misc_feature") else default
        if label == "misc_feature":
            label = None
        return label


def plot_local_gc_content(sequence, window_size, ax=None):
    """Plot the local (windowed) GC content of the sequence

    Parameters
    ----------

    sequence
      An ATGC string of a sequence

    window_size
      Size of the GC window

    ax
      Matplotlib ax on which to plot. If none is provided, one is created.
    """
    if ax is None:
        fig, ax = plt.subplots(1)

    def gc_content(sequence):
        return 100.0 * len([c for c in sequence if c in "GC"]) / len(sequence)
    yy = [gc_content(sequence[i:i + window_size])
          for i in range(len(sequence) - window_size)]
    xx = np.arange(len(sequence) - window_size) + 25
    ax.fill_between(xx, yy, alpha=0.3)
    for x in range(10, 100, 10):
        ax.axhline(x, lw=0.5, alpha=0.3)
    ax.set_ylabel("GC(%)")
    ax.set_xlim(0, len(sequence))
    ax.set_ylim(ymin=0)


def plot_local_gc_with_features(record, window_size, axes=None):
    """Plot the local GC content curve below a plot of the record's features.

    Parameters
    ----------

    record
      A biopython record

    window_size
      Size of the GC window

    axes
      A tuple of two matplotlib axes. If none provided, they will be created.

    """
    if axes is None:
        fig, axes = plt.subplots(2, figsize=(20, 5.5), sharex=True)
    graphic_record = BiopythonTranslator().translate_record(record)
    graphic_record.plot(ax=axes[0], with_ruler=False)
    ax = plot_local_gc_content(record.seq, 40, ax=axes[1])
    return ax
