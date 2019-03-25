"""Misc. plotting and reporting methods, some of which are really arbitrary.
"""

from copy import deepcopy

from Bio.SeqRecord import SeqRecord
import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from dna_features_viewer import BiopythonTranslator
    DFV_AVAILABLE = True
except:
    BiopythonTranslator = object
    DFV_AVAILABLE = False

def install_extras_message(libname):
    return (
        "Could not load %s (is it installed ?). You can install it separately "
        " with:  pip install %s\n\n"
        "Install all dependencies for generating DNA Chisel reports with:"
        "\n\npip install dnachisel[reports]" % (
            libname, libname.lower().replace(" ", "_")))

from ..DnaOptimizationProblem import DnaOptimizationProblem
from ..biotools import gc_content
from ..SequencePattern import RepeatedKmerPattern, HomopolymerPattern
from ..builtin_specifications import (EnforceGCContent, AvoidPattern,
                                      AvoidHairpins)


def plot_constraint_breaches(constraint, sequence, title=None, ax=None,
                             show_locations=False, show_feature_labels=False):
    """Plot all breaches in a sequence for a single constraint.

    Parameters
    -----------
    constraint
      A DnaChisel Specification, to be used as a constraint.

    sequence
      An ATGC string, or a record. If it is a record, its annotations will also
      be shown in the plot.

    title
      Title to be added to the plot

    ax
      Matplotlib ax on which to plot the result. If no ax is provided, one
      will be created, and it will be returned at the end

    show_locations
      Whether or not to label the breach features with their coordinates

    show_feature_labels
      Whether or not to show the labels of all features.


    Returns
    -------
    ax
      The matplotlib ax on which the plot was drawn

    """
    if not DFV_AVAILABLE:
      raise ImportError(install_extras_message("DNA Features Viewer"))
    class Translator(BiopythonTranslator):
        """A Biopython record translator for DNA Features Viewer.

        This translator produces label-free plots.
        """

        default_feature_color = "#ffaaaa"

        def compute_feature_label(self, f):
            if f.type != "original" and show_locations:
                return str(int(f.location.start))
            elif show_feature_labels and f.type == "original":
                return BiopythonTranslator.compute_feature_label(f)
            else:
                return None

        def compute_feature_color(self, f):
            return "#ffffff" if f.type == "original" else "#ff9999"

    if isinstance(sequence, SeqRecord):
        record = sequence
        sequence = str(sequence.seq).upper()
    else:
        record = None

    problem = DnaOptimizationProblem(sequence, constraints=[constraint])
    evals = problem.constraints_evaluations().filter("failing")
    breaches_record = problem.to_record(
        with_original_features=False, with_objectives=False,
        with_constraints=False
    )
    breaches_record.features = evals.locations_as_features(
        with_specifications=False, merge_overlapping=True)

    if record is not None:
        for feature in record.features:
            feature = deepcopy(feature)
            feature.type = "original"
            breaches_record.features.append(feature)

    translator = Translator()
    graphic_record = translator.translate_record(breaches_record)
    ax, _ = graphic_record.plot(ax=ax)
    ax.set_ylim(bottom=-1, top=min(5, ax.get_ylim()[1]))
    if title is not None:
        ax.set_title(title, fontweight="bold", loc="left")
    return ax

def plot_gc_content_breaches(sequence, window=70, mini=0.35, maxi=0.65,
                             ax=None, title=None):
    """Plot a profile of GC content along the sequence.
    The regions with out-of-bound GC are highlighted.

    Returns the Matplotlib ax on which the figure was plotted.

    Parameters
    ----------
    sequence
      An ATGC string

    window
      Number of nucleotides to use for local averaging of GC content

    mini
      minimal allowed proportion of gc (between 0 and 1)

    maxi
      maximal allowed proportion of gc (between 0 and 1)

    ax
      Matplotlib ax on which to plot the figure. If none is provided, one
      is created and returned in the end.

    title
      Title to be added to the figure.

    """
    gc = gc_content(sequence, window_size=window)
    xx = np.arange(len(gc)) + window / 2

    if ax is None:
        fig, ax = plt.subplots(1, figsize=(12, 3))
    gc_inbound = +gc
    gc_inbound[(gc < mini) | (gc > maxi)] = np.nan
    ax.plot(xx, gc_inbound, alpha=0.2, c='b')

    for limit, whr in (mini, gc < mini), (maxi, gc > maxi):
        ax.axhline(limit, c="k", lw=0.5, ls=":")
        gc_bad = +gc
        gc_bad[1 - whr] = np.nan
        ax.plot(xx, gc_bad, alpha=0.1, c='r')
        ax.fill_between(xx, gc, y2=limit, where=whr, alpha=0.6,
                        facecolor="r")

    ax.set_xlim(0, len(gc))
    ax.set_ylim(0, 1)
    if title is not None:
        ax.set_title(title, fontweight="bold", loc="left")
    return ax


def plot_sequence_manufacturability_difficulties(sequence):
    """Plot a series of manufacturability factor checks, on a single figure.

    The factors considered are mostly arbitrary. Most users won't be interested
    by this function. The function looks for presence of Type IIs restriction
    sites, homopolymers, hairpins, extreme GC content, repeats.
    The breaches for each constraint type are plotted on separate axis.

    Parameters
    ----------
    sequence
      An ATGC string or a Biopythn sequence record (at which case the features
      will be drawn along with the constraints breaches).

    Returns
    -------

    axes
      List of the different axes (in a same figure) on which the plots were
      drawn
    """
    if isinstance(sequence, SeqRecord):
        record = sequence
        sequence = str(sequence.seq).upper()
    else:
        record = sequence

    nplots = 8
    fig, axes = plt.subplots(nplots, 1, figsize=(10, 1.4 * nplots),
                             sharex=True, facecolor="white")

    mini, maxi, window = 0.25, 0.80, 50
    plot_gc_content_breaches(
        sequence, window=window, mini=mini,
        maxi=maxi, ax=axes[0],
        title="GC content (window= %d)" % window
    )

    plot_constraint_breaches(
        EnforceGCContent(mini=mini, maxi=maxi, window=window),
        record, ax=axes[1],
        title="Zones of extreme GC content (Gen9-type short window)"
    )

    plot_constraint_breaches(
        AvoidPattern("BsmBI_site"), sequence,
        title="BsmBI sites", ax=axes[2]
    )

    plot_constraint_breaches(
        AvoidPattern("BsaI_site"),
        record, title="BsaI sites", ax=axes[3]
    )

    plot_constraint_breaches(
        AvoidPattern("BbsI_site"),
        record, title="BbsI sites", ax=axes[4]
    )

    for l, n in [("A", 9), ("T", 9), ("G", 6), ("C", 6)]:
        constraint = AvoidPattern(HomopolymerPattern(l, n))
        plot_constraint_breaches(
            constraint, record,
            title="Homopolymers (6+ G or C | 9+ A or T)",
            ax=axes[5]
        )

    for kmer_size, n_repeats in (3, 5), (2, 9):
        pattern = RepeatedKmerPattern(n_repeats=n_repeats, kmer_size=kmer_size)
        constraint = AvoidPattern(pattern)
        plot_constraint_breaches(
            constraint, record,
            title="Repeats (5 x 3bp, 9 x 2bp)", ax=axes[6]
        )

    plot_constraint_breaches(
        AvoidHairpins(stem_size=20, hairpin_window=200),
        record, title="Hairpins", ax=axes[7]
    )

    fig.tight_layout()
    return axes
