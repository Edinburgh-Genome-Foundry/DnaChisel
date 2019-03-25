"""Methods to generate optimization reports."""

import os
import textwrap

from Bio import SeqIO

import flametree

from ..biotools import (sequence_to_biopython_record,
                        find_specification_in_feature)
from ..version import __version__

def install_extras_message(libname):
    return (
        "Could not load %s (is it installed ?). You can install it separately "
        " with:  pip install %s\n\n"
        "Install all dependencies for generating DNA Chisel reports with:"
        "\n\npip install dnachisel[reports]" % (
            libname, libname.lower().replace(" ", "_")))

try:
    from sequenticon import sequenticon
    SEQUENTICON_AVAILABLE = True
except:
    SEQUENTICON_AVAILABLE = False

MATPLOTLIB_AVAILABLE = DFV_AVAILABLE = False
try:
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    MATPLOTLIB_AVAILABLE = True
    from dna_features_viewer import BiopythonTranslator
    DFV_AVAILABLE = True
except ImportError:
    class BiopythonTranslator:
        "Class unavailable. Install DNA Features Viewer."
        def __init__(self):
            raise ImportError(install_extras_message("DNA Features Viewer"))

try:
    from geneblocks import DiffBlocks
    GENEBLOCKS_AVAILABLE = True
except:
    GENEBLOCKS_AVAILABLE = False

try:
    from pdf_reports import ReportWriter
    import pdf_reports.tools as pdf_tools
    PDF_REPORTS_AVAILABLE = True
except:
    def ReportWriter(*a, **kw):
        return None
    PDF_REPORTS_AVAILABLE = False

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
ASSETS_DIR = os.path.join(THIS_DIR, "assets")
TITLE_FONTDICT = fontdict = dict(size=14, weight="bold")

report_writer = ReportWriter(
    dnachisel_logo_url=os.path.join(ASSETS_DIR, 'logo.png'),
    version=__version__,
    default_stylesheets=(os.path.join(ASSETS_DIR, "style.css"),)
)

install_reports_extra_message =(
    "Could not load %s (is it installed ?). You can install all "
    "dependencies for generating reports in DNA Chisel with this command:\n\n "
    "pip install dnachisel[reports]")

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
            feature for feature in features
            if "".join(feature.qualifiers.get("is_edit", "false")) != 'true'
        ]


    def compute_feature_color(self, f):
        color = f.qualifiers.get('color', None)
        if color is not None:
            if isinstance(color, list):
                color = color[0]
            return color

        if f.type == "misc_feature":
            specification = find_specification_in_feature(f)
            if specification is None:
                return "#f4df42"
            else:
                return self.feature_prefixes_colors.get(specification[0],
                                                        "#f4df42")
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

def write_no_solution_report(target, problem, error):
    """Write a report on incompatibility found in the problem's constraints.

    The report comprises a PDF of plots of the sequence (global constraints,
    local constraints around the problem) and an annotated genbank.

    Parameters
    ----------
    target
      Either a path to a folder, or a path to a zip archive, or "@memory" to
      return raw data of a zip archive containing the report.

    problem
      A DnaOptimizationProblem

    error
      A NoSolutionError (carries a message and a location)
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError(install_extras_message("Matplotlib"))
    if isinstance(target, str):
        root = flametree.file_tree(target, replace=True)
    else:
        root = target
    translator = SpecAnnotationsTranslator()
    with PdfPages(root._file("plots.pdf").open("wb")) as pdf_io:

        # PLOT GLOBAL LOCATION OF ERROR

        record = problem.to_record()
        translator = SpecAnnotationsTranslator()
        graphical_record = translator.translate_record(record)
        ax, _ = graphical_record.plot(figure_width=min(20, 0.3*len(record)))
        if len(record) < 60:
            graphical_record.plot_sequence(ax)
        if error.location is None:
            raise error
        start, end, strand = error.location.to_tuple()
        ax.fill_between([start, end], -10, 10, zorder=-1000,
                        facecolor='#ffeeee')
        title = "\n".join(textwrap.wrap(
            "No solution found in zone [%d, %d]: %s" %
            (start, end, str(error)), width=120)
        )
        ax.set_title(title, fontdict=TITLE_FONTDICT)
        pdf_io.savefig(ax.figure, bbox_inches="tight", alpha=0.5)
        plt.close(ax.figure)

        # PLOT LOCAL CONSTRAINTS BREACHES

        evals = error.problem.constraints_evaluations()
        record = error.problem.to_record(
            with_original_spec_features=False,
            with_constraints=False, with_objectives=False)
        record.features += evals.filter('passing') \
                                .success_and_failures_as_features()
        record.features += evals.filter('failing') \
                                .locations_as_features(label_prefix="BREACH")
        start = max(0, error.location.start - 5)
        end = min(len(record), error.location.end + 4)
        graphical_record = translator.translate_record(record)
        graphical_record = graphical_record.crop((start, end))
        ax, _ = graphical_record.plot(figure_width=min(20, 0.3*(end - start)))
        graphical_record.plot_sequence(ax)
        ax.set_title("Local constraints breaches in [%d, %d]" % (start, end) +
                     "     (green = passing constraints)",
                     fontdict=TITLE_FONTDICT)
        pdf_io.savefig(ax.figure, bbox_inches="tight", alpha=0.5)
        plt.close(ax.figure)

        # WRITE GENBANK

        record = problem.to_record(with_original_spec_features=False,
                                   with_constraints=True,
                                   with_objectives=True)
        evals = problem.constraints_evaluations()
        record.features += evals.filter('passing') \
                                .success_and_failures_as_features()
        record.features += evals.filter('failing') \
                                .locations_as_features(label_prefix="BREACH")
        SeqIO.write(record, root._file("constraints breaches.gb").open("w"),
                    "genbank")
    root._file('logs.txt').write(problem.logger.dump_logs())

    # returns zip data if target == '@memory'
    if isinstance(target, str):
        return root._close()


def write_optimization_report(target, problem, project_name="unnammed",
                              constraints_evaluations=None,
                              objectives_evaluations=None,
                              figure_width=20, max_features_in_plots=300):
    """Write an optimization report with a PDF summary, plots, and genbanks.

    Parameters
    ----------
    target
      Path to a directory or zip file, or "@memory" for returning raw data of
      a zip file created in-memory.
    
    problem
      A DnaOptimizationProblem to be solved and optimized
    
    project_name
      Name of the project that will appear on the PDF report
    
    constraints_evaluations
      Precomputed constraints evaluations. If None provided, they will be
      computed again from the problem.
    
    objectives_evaluations
      Precomputed objectives evaluations. If None provided, they will be
      computed again from the problem.
      
    
    figure_width
      Width of the report's figure, in inches. The more annotations there will
      be in the figure, the wider it should be. The default should work for
      most cases.
    
    max_features_in_plots
      Limit to the number of features to plot (plots with thousands of features
      may take ages to plot)

    """
    if not PDF_REPORTS_AVAILABLE:
        raise ImportError(install_extras_message("PDF Reports"))
    if not SEQUENTICON_AVAILABLE:
        raise ImportError(install_extras_message("Sequenticon"))
    if constraints_evaluations is None:
        constraints_evaluations = problem.constraints_evaluations()
    if objectives_evaluations is None:
        objectives_evaluations = problem.objectives_evaluations()
    if isinstance(target, str):
        root = flametree.file_tree(target, replace=True)
    else:
        root = target
    translator = SpecAnnotationsTranslator()
    # CREATE FIGURES AND GENBANKS
    diffs_figure_data = None
    sequence_before = sequence_to_biopython_record(problem.sequence_before)
    if GENEBLOCKS_AVAILABLE:
        sequence_after = problem.to_record()
        contract_under = max(3, int(len(sequence_after) / 10))
        diffs = DiffBlocks.from_sequences(sequence_before, sequence_after,
                                          use_junk_over=50,
                                          contract_under=contract_under)
        _, diffs_ax = diffs.plot()
        diffs_figure_data = pdf_tools.figure_data(diffs_ax.figure, fmt='svg')
        plt.close(diffs_ax.figure)

    with PdfPages(root._file("before_after.pdf").open("wb")) as pdf_io:

        figures_data = [
            (
                "Before",
                sequence_before,
                problem.constraints_before,
                problem.objectives_before,
                []
            ),
            (
                "After",
                sequence_to_biopython_record(problem.sequence),
                constraints_evaluations,
                objectives_evaluations,
                problem.sequence_edits_as_features()
            )
        ]

        plot_height = None
        for (title, record, constraints, objectives, edits) in figures_data:

            full_title = (
                "{title}:        {nfailing} constraints failing (in red)"
                "        Total Score: {score:.01E} {bars}").format(
                title=title, score=objectives.scores_sum(),
                nfailing=len(constraints.filter("failing").evaluations),
                bars="" if (title == "Before") else
                "       (bars indicate edits)"
            )
            ax = None
            if title == "After":
                record.features += edits
                graphical_record = translator.translate_record(record)
                fig, ax = plt.subplots(1, figsize=(figure_width, plot_height))
                graphical_record.plot(ax=ax, level_offset=-0.3)
                record.features = []

            record.features += constraints.success_and_failures_as_features()
            record.features += objectives.success_and_failures_as_features()

            graphical_record = translator.translate_record(record)
            ax, _ = graphical_record.plot(ax=ax, figure_width=figure_width)
            ax.set_title(full_title, loc="left", fontdict=TITLE_FONTDICT)
            plot_height = ax.figure.get_size_inches()[1]
            pdf_io.savefig(ax.figure, bbox_inches="tight")
            plt.close(ax.figure)

            record.features += edits
            breaches_locations = \
                constraints.filter("failing") \
                           .locations_as_features(label_prefix="Breach from",
                                                  merge_overlapping=True)
            record.features += breaches_locations

            SeqIO.write(record, root._file(title.lower() + ".gb").open("w"),
                        "genbank")

            if breaches_locations != []:
                record.features = breaches_locations
                graphical_record = translator.translate_record(record)
                if len(graphical_record.features) > max_features_in_plots:
                    features = sorted(graphical_record.features,
                                      key=lambda f: f.start - f.end)
                    new_ft = features[:max_features_in_plots]
                    graphical_record.features = new_ft
                    message = "(only %d features shown)" % \
                              max_features_in_plots
                else:
                    message = ""
                ax, _ = graphical_record.plot(figure_width=figure_width)
                ax.set_title(title + ": Constraints breaches locations"
                             + message, loc="left", fontdict=TITLE_FONTDICT)
                pdf_io.savefig(ax.figure, bbox_inches="tight")
                plt.close(ax.figure)

    # CREATE PDF REPORT
    html = report_writer.pug_to_html(
        path=os.path.join(ASSETS_DIR, "optimization_report.pug"),
        project_name=project_name,
        problem=problem,
        constraints_evaluations=constraints_evaluations,
        objectives_evaluations=objectives_evaluations,
        edits=sum(len(f) for f in edits),
        diffs_figure_data=diffs_figure_data,
        sequenticons={
            label: sequenticon(seq, output_format="html_image", size=24)
            for label, seq in [("before", problem.sequence_before),
                               ("after", problem.sequence)]
        }
    )
    problem.to_record(root._file("final_sequence.gb").open("w"),
                      with_constraints=False,
                      with_objectives=False)

    report_writer.write_report(html, root._file("Report.pdf"))
    if isinstance(target, str):
        return root._close()