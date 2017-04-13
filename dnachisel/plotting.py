"""Misc. plotting and reporting methods"""
from copy import deepcopy

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from dna_features_viewer import BiopythonTranslator
    DFV_AVAILABLE = True
except:
    BiopythonTranslator = object
    DFV_AVAILABLE = False
import numpy as np

from .DnaOptimizationProblem import DnaOptimizationProblem


class GraphicTranslator(BiopythonTranslator):
    """Translator of DnaChisel feature-constraints for DNA Features Viewer"""
    @staticmethod
    def compute_feature_color(f):
        if f.type == "misc_feature":
            label = BiopythonTranslator.compute_feature_label(f)
            return {
                "@": "#ce5454",
                "~": "#e5be54",
                "#": "#8edfff",
                "!": "#fcff75",

            }.get(label[0], "#f4df42")
        else:
            return "#d9e0fc"

    @staticmethod
    def compute_feature_label(f):
        default = BiopythonTranslator.compute_feature_label(f)
        return None if (f.type != "misc_feature") else default


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
        return 100.0*len([c for c in sequence if c in "GC"]) / len(sequence)
    yy = [gc_content(sequence[i:i+window_size])
          for i in range(len(sequence)-window_size)]
    xx = np.arange(len(sequence)-window_size)+25
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

def plot_constraints_breaches(record, constraints, color="red",
                              ax=None, figure_width=20):
    """Plot the record and highlight locations of constraints breaches.

    Requires "DNA Features Viewer" installed.

    Parameters
    ----------

    record
      A Biopython record

    constraints
      A list or tuple of DnaChisel objectives to use as constraints

    color
      Color for highlighting the constraint breaches in the plot. All
      other features will appear in white.

    ax
       Matplotlib ax on which to plot. If none is provided, one is created.

    figure_width
      Width of the final figure in inches, if no ax is provided


    """
    if not DFV_AVAILABLE:
        raise ImportError("Plotting constraint breaches requires "
                          "Matplotlib and dna_features_viewer installed.")
    class MyTranslator(BiopythonTranslator):
        def compute_feature_color(self, f):
            return color if (f.type == "breach") else "white"

    new_record = deepcopy(record)
    sequence = str(record.seq).upper()
    pb = DnaOptimizationProblem(sequence, constraints=constraints)
    rec = pb.constraints_breaches_as_biopython_record("breach")
    for f in rec.features:
        new_record.features.append(f)
    gr_record = MyTranslator().translate_record(new_record)
    ax, _ = gr_record.plot(ax=ax, figure_width=20)
    return ax

def make_constraints_breaches_pdf(constraints_sets, record, pdf_path):
    """Plot the record and highlight locations of constraints breaches.

    Requires "DNA Features Viewer" installed.

    Parameters
    ----------

    record
      A Biopython record

    constraints
      A list or tuple of DnaChisel objectives to use as constraints

    color
      Color for highlighting the constraint breaches in the plot. All
      other features will appear in white.

    ax
       Matplotlib ax on which to plot. If none is provided, one is created.

    figure_width
      Width of the final figure in inches, if no ax is provided


    """
    if not DFV_AVAILABLE:
        raise ImportError("Package ``dna_features_viewer`` must be installed"
                          "to use method ``make_constraints_breaches_pdf``")
    with PdfPages(pdf_path) as pdf:
        for title, constraints in constraints_sets.items():
            ax = plot_constraints_breaches(record, constraints,
                                           figure_width=20)
            ax.set_title(title, fontsize=16, fontweight="bold")
            pdf.savefig(ax.figure, bbox_inches="tight")
            plt.close(ax.figure)
