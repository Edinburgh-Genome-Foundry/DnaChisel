"""Methods to generate optimization reports."""

import os
import textwrap
from collections import OrderedDict
import hashlib

from Bio import SeqIO
import flametree
import numpy as np

from ..biotools import (
    sequence_to_biopython_record,
    find_specification_label_in_feature,
)
from ..version import __version__
from .SpecAnnotationsTranslator import SpecAnnotationsTranslator
from .tools import install_extras_message
from ..Location import Location

try:
    import pandas
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

try:
    from sequenticon import sequenticon

    SEQUENTICON_AVAILABLE = True
except ImportError:
    SEQUENTICON_AVAILABLE = False

MATPLOTLIB_AVAILABLE = False
try:
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    pass


try:
    from geneblocks import DiffBlocks

    GENEBLOCKS_AVAILABLE = True
except ImportError:
    GENEBLOCKS_AVAILABLE = False

try:
    from pdf_reports import ReportWriter
    import pdf_reports.tools as pdf_tools

    PDF_REPORTS_AVAILABLE = True
except ImportError:

    def ReportWriter(*a, **kw):
        return None

    PDF_REPORTS_AVAILABLE = False

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
ASSETS_DIR = os.path.join(THIS_DIR, "assets")
TITLE_FONTDICT = fontdict = dict(size=14, weight="bold")

report_writer = ReportWriter(
    dnachisel_logo_url=os.path.join(ASSETS_DIR, "logo.png"),
    version=__version__,
    default_stylesheets=(os.path.join(ASSETS_DIR, "style.css"),),
)

install_reports_extra_message = (
    "Could not load %s (is it installed ?). You can install all "
    "dependencies for generating reports in DNA Chisel with this command:\n\n "
    "pip install dnachisel[reports]"
)


def write_no_solution_report(
    target, problem, error, file_content=None, file_path=None
):
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

    # TRANSFER THE ORIGINAL FILE
    file_hash = None
    if file_path is not None:
        if file_content is None:
            with open(file_path, "rb") as f:
                file_content = f.read()
        basename = os.path.basename(file_path)
        file_hash = hashlib.md5(file_content).hexdigest()[:8]
        root._file("_".join([file_hash, basename])).write(file_content)

    translator = SpecAnnotationsTranslator()
    with PdfPages(root._file("plots.pdf").open("wb")) as pdf_io:

        # PLOT GLOBAL LOCATION OF ERROR

        record = problem.to_record()
        translator = SpecAnnotationsTranslator()
        graphical_record = translator.translate_record(record)
        ax, _ = graphical_record.plot(figure_width=min(20, 0.3 * len(record)))
        if len(record) < 60:
            graphical_record.plot_sequence(ax)
        if error.location is None:
            raise error
        start, end, strand = error.location.to_tuple()
        ax.fill_between(
            [start, end], -10, 10, zorder=-1000, facecolor="#ffcccc"
        )
        title = "\n".join(
            textwrap.wrap(
                "No solution found in zone [%d, %d]:%s"
                % (start, end, str(error)),
                width=120,
            )
        )
        ax.set_title(title, fontdict=TITLE_FONTDICT)
        pdf_io.savefig(ax.figure, bbox_inches="tight", alpha=0.5)
        plt.close(ax.figure)

        # CREATE AND SAVE THE LOCAL CONSTRAINTS BREACHES RECORD

        record = error.problem.to_record(
            with_original_spec_features=False,
            with_constraints=False,
            with_objectives=False,
        )

        start = max(0, error.location.start - 5)
        end = min(len(record), error.location.end + 4)
        focus_location = Location(start, end)

        def is_in_focus(location):
            return location.overlap_region(focus_location) is not None

        evals = error.problem.constraints_evaluations()
        passing = evals.filter("passing")
        record.features += passing.success_and_failures_as_features()
        failing = evals.filter("failing")
        record.features += failing.locations_as_features(
            label_prefix="BREACH", locations_filter=is_in_focus
        )
        SeqIO.write(
            record,
            root._file("local_constraints_breaches.gb").open("w"),
            "genbank",
        )

        # CREATE A FIGURE OF THE LOCAL CONSTRAINTS BREACHES AS A NEW PDF PAGE

        graphical_record = translator.translate_record(record)
        graphical_record = graphical_record.crop((start, end))
        figure_width = min(20, 0.3 * (end - start))
        ax, _ = graphical_record.plot(figure_width=figure_width)
        graphical_record.plot_sequence(ax)
        ax.set_title(
            "Local constraints breaches in [%d, %d]" % (start, end)
            + "     (green = passing constraints)",
            fontdict=TITLE_FONTDICT,
        )
        ax.set_ylim(top=ax.get_ylim()[1] + 1)
        pdf_io.savefig(ax.figure, bbox_inches="tight", alpha=0.5)
        plt.close(ax.figure)

    root._file("logs.txt").write(problem.logger.dump_logs())

    # returns zip data if target == '@memory'
    if isinstance(target, str):
        return root._close()


def constraints_before_after_dataframe(problem, constraints_evaluations=None):
    if not PANDAS_AVAILABLE:
        raise ImportError("Install pandas to use this method.")
    if constraints_evaluations is None:
        constraints_evaluations = problem.constraints_evaluations()
    edits = problem.sequence_edits_as_array()

    def constraint_record(evaluation_before, evaluation_after):
        constraint = evaluation_before.specification
        start, end, _ = constraint.location.to_tuple()
        edits_sum = edits[start:end].sum()
        edits_percent = 100 * edits_sum / (end - start)
        label = constraint.label(use_short_form=True, with_location=False)
        return OrderedDict(
            [
                ("constraint", label),
                ("start", start),
                ("end", end),
                ("before", "PASS" if evaluation_before.passes else "FAIL"),
                ("after", "PASS" if evaluation_after.passes else "FAIL"),
                ("edits", edits_sum),
                ("% edited", np.round(edits_percent, 2)),
            ]
        )

    dataframe = pandas.DataFrame.from_records(
        [
            constraint_record(before, after)
            for (before, after) in zip(
                problem.constraints_before, constraints_evaluations
            )
        ]
    )
    if len(dataframe):
        dataframe = dataframe.sort_values(by="start")
    return dataframe


def objectives_before_after_dataframe(problem, objectives_evaluations=None):
    if objectives_evaluations is None:
        objectives_evaluations = problem.objectives_evaluations()
    edits = problem.sequence_edits_as_array()

    def objective_record(evaluation_before, evaluation_after):
        objective = evaluation_before.specification
        start, end, _ = objective.location.to_tuple()
        edits_sum = edits[start:end].sum()
        edits_percent = 100 * edits_sum / (end - start)
        label = objective.label(use_short_form=True, with_location=False)
        return OrderedDict(
            [
                ("objective", label),
                ("boost", objective.boost),
                ("start", start),
                ("end", end),
                ("before", evaluation_before.score_to_formatted_string),
                ("after", evaluation_after.score_to_formatted_string),
                ("edits", edits_sum),
                ("% edited", np.round(edits_percent, 2)),
            ]
        )

    dataframe = pandas.DataFrame.from_records(
        [
            objective_record(before, after)
            for (before, after) in zip(
                problem.objectives_before, objectives_evaluations
            )
        ]
    )
    if len(dataframe):
        dataframe = dataframe.sort_values(by="start")
    return dataframe


def plot_optimization_changes(problem):
    if not GENEBLOCKS_AVAILABLE:
        raise ImportError("Install Geneblocks to use plot_differences()")
    sequence_before = sequence_to_biopython_record(problem.sequence_before)
    sequence_after = problem.to_record()
    diffs = DiffBlocks.from_sequences(sequence_before, sequence_after)
    span = max(2, len(sequence_after) / 20)
    diffs = diffs.merged(
        blocks_per_span=(3, span), replace_gap=span / 2, change_gap=span / 2
    )
    _, diffs_ax = diffs.plot(
        translator_class=SpecAnnotationsTranslator,
        annotate_inline=True,
        figure_width=15,
    )
    return diffs_ax


def write_optimization_report(
    target,
    problem,
    project_name="unnamed",
    plot_figure=True,
    constraints_evaluations=None,
    objectives_evaluations=None,
    figure_width=20,
    max_features_in_plots=300,
    file_path=None,
    file_content=None,
):
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
    
    file_path
      Path to the file from which the problem was created
    
    

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

    # TRANSFER THE ORIGINAL FILE
    file_hash = None
    if file_path is not None:
        if file_content is None:
            with open(file_path, "rb") as f:
                file_content = f.read()
        basename = os.path.basename(file_path)
        file_hash = hashlib.md5(file_content).hexdigest()[:8]
        root._file("_".join([file_hash, basename])).write(file_content)

    # CREATE FIGURES AND GENBANKS
    diffs_figure_data = None
    if GENEBLOCKS_AVAILABLE and plot_figure:
        diffs_ax = plot_optimization_changes(problem)
        diffs_figure_data = pdf_tools.figure_data(diffs_ax.figure, fmt="svg")
        plt.close(diffs_ax.figure)

    # GENERATE AND SAVE THE CONSTRAINTS SUMMARY

    constraints_before_after = constraints_before_after_dataframe(
        problem=problem, constraints_evaluations=constraints_evaluations
    )
    filename = "constraints_before_and_after.csv"
    constraints_before_after.to_csv(
        root._file(filename).open("w"), index=False
    )

    # GENERATE AND SAVE THE OBJECTIVES SUMMARY

    objectives_before_after = objectives_before_after_dataframe(
        problem=problem, objectives_evaluations=objectives_evaluations
    )
    filename = "objectives_before_and_after.csv"
    objectives_before_after.to_csv(root._file(filename).open("w"), index=False)

    # CREATE PDF REPORT
    html = report_writer.pug_to_html(
        path=os.path.join(ASSETS_DIR, "optimization_report.pug"),
        project_name=project_name,
        problem=problem,
        constraints_evaluations=constraints_evaluations,
        objectives_evaluations=objectives_evaluations,
        constraints_before_after=constraints_before_after,
        objectives_before_after=objectives_before_after,
        edits=problem.sequence_edits_as_array().sum(),
        diffs_figure_data=diffs_figure_data,
        file_hash=file_hash,
        sequenticons={
            label: sequenticon(seq, output_format="html_image", size=24)
            for label, seq in [
                ("before", problem.sequence_before),
                ("after", problem.sequence),
            ]
        },
    )
    report_writer.write_report(html, root._file("Report.pdf"))

    # CREATE THE "SEQUENCE EDITS" REPORT

    record = problem.to_record(with_sequence_edits=True)
    breaches = problem.constraints_before.filter("failing")
    breaches_locations = breaches.locations_as_features(
        label_prefix="Breach from", merge_overlapping=True
    )
    record.features += breaches_locations
    SeqIO.write(
        record, root._file("final_sequence_with_edits.gb").open("w"), "genbank"
    )

    # CREATE THE "FINAL SEQUENCE" REPORT

    problem.to_record(
        root._file("final_sequence.gb").open("w"),
        with_constraints=False,
        with_objectives=False,
    )

    if isinstance(target, str):
        return root._close()

