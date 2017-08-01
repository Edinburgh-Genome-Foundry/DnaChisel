"""Methods to generate optimization reports."""

import os
import textwrap

from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import jinja2
import weasyprint
import flametree

from ..biotools import sequence_to_biopython_record, crop_record
from ..plotting_tools import SpecAnnotationsTranslator
from ..version import __version__
from ..DnaOptimizationProblem import (DnaOptimizationProblem, NoSolutionError,
                                      NoSolutionFoundError)

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATES_DIR = os.path.join(THIS_DIR, "templates")
TITLE_FONTDICT = fontdict = dict(size=14, weight="bold")


def optimization_with_report(target, problem=None, record=None,
                             project_name="unnamed",
                             specifications_dict='default',
                             **solver_parameters):
    """Optimize a sequence and write a multi-file report.

    The report's content may vary depending on the success of the optimization

    Parameters
    ----------
    target
      Either a path to a folder that will containt the report, or a path to
      a zip archive, or "@memory" to return raw data of a zip archive
      containing the report.

    problem
      A DnaOptimizationProblem. A ``record`` can be provided instead

    record
      Either a Biopython record or the path to a Genbank file, provided to
      serve as a problem definition. The record should have specifications
      starting with @ or ~ as explained in the documentation.

    project_name
      Project name to write on PDF reports

    specifications_dict
      a dictionnary {'@SpecName': SpecClass} allowing the parser to map
      specifications found in annotations to implemented specifications
      classes. By default, the dictionnary contains the built-in objectives
      of DNA chisel, but an extended dictionary can be provided to support
      custom specifications.

    **solver_parameters
      Extra keyword arguments passed to ``problem.solve_constraints()``

    Returns
    -------
    (success, message, zip_data)
      Triplet where success is True/False, message is a one-line string
      summary indication whether some clash was found, or some solution, or
      maybe no solution was found because the random searches were too short

    """
    if problem is None:
        problem = DnaOptimizationProblem.from_record(
            record, specifications_dict=specifications_dict)

    try:
        problem.progress_logger(message="Solving constraints")
        problem.solve_constraints(**solver_parameters)
    except NoSolutionError as error:
        problem.progress_logger(message="No solution found: making report")
        data = write_no_solution_report(target, problem, error)
        start, end, s = error.location.as_tuple()
        message = ("No solution found in zone [%d, %d]: %s." %
                   (start, end, str(error)))
        return False, message, data
    except NoSolutionFoundError as error:
        problem.progress_logger(
            message="Failed to solve constraints. Generating report.")
        data = write_optimization_report(
            target, problem, project_name=project_name)
        return False, "No solution found before the end of search.", data

    problem.progress_logger(message="Now optimizing the sequence")
    problem.maximize_objectives()
    problem.progress_logger(message="Success ! Generating report.")
    data = write_optimization_report(
        target, problem, project_name=project_name)
    return True, "Optimization successful.", data


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
    root = flametree.file_tree(target, replace=True)
    translator = SpecAnnotationsTranslator()
    with PdfPages(root._file("plots.pdf").open("wb")) as pdf_io:

        # PLOT GLOBAL LOCATION OF ERROR

        record = problem.to_record()
        translator = SpecAnnotationsTranslator()
        graphical_record = translator.translate_record(record)
        ax, _ = graphical_record.plot(figure_width=20)
        start, end, strand = error.location.to_tuple()
        ax.fill_between([start, end], -10, 10)
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
            with_original_objective_features=False,
            with_constraints=False, with_objectives=False)
        record.features += evals.filter('passing') \
                                .success_and_failures_as_features()
        record.features += evals.filter('failing') \
                                .locations_as_features(label_prefix="BREACH")
        start, end = error.location.start-5, error.location.end+4
        local_record = crop_record(record, start, end)
        graphical_record = translator.translate_record(local_record)
        ax, _ = graphical_record.plot(figure_width=20)
        ax.set_title("Local constraints breaches in [%d, %d]" % (start, end) +
                     "     (green = passing constraints)",
                     fontdict=TITLE_FONTDICT)
        pdf_io.savefig(ax.figure, bbox_inches="tight", alpha=0.5)
        plt.close(ax.figure)

        # WRITE GENBANK

        record = problem.to_record(with_original_objective_features=False,
                                   with_constraints=True,
                                   with_objectives=True)
        evals = problem.constraints_evaluations()
        record.features += evals.filter('passing') \
                                .success_and_failures_as_features()
        record.features += evals.filter('failing') \
                                .locations_as_features(label_prefix="BREACH")
        SeqIO.write(record, root._file("constraints breaches.gb").open("w"),
                    "genbank")

    # returns zip data if target == '@memory'
    return root._close()


def write_optimization_report(target, problem, project_name="unnammed",
                              constraints_evaluations=None,
                              objectives_evaluations=None,
                              figure_width=20, max_features_in_plots=300):
    """Write an optimization report with a PDF summary, plots, and genbanks.

    Parameters
    ----------




    """
    if constraints_evaluations is None:
        constraints_evaluations = problem.constraints_evaluations()
    if objectives_evaluations is None:
        objectives_evaluations = problem.objectives_evaluations()
    root = flametree.file_tree(target, replace=True)
    translator = SpecAnnotationsTranslator()
    # CREATE FIGURES AND GENBANKS

    with PdfPages(root._file("before_after.pdf").open("wb")) as pdf_io:

        figures_data = [
            (
                "Before",
                sequence_to_biopython_record(problem.sequence_before),
                problem.constraints_before,
                problem.objectives_before,
                []
            ),
            (
                "After",
                sequence_to_biopython_record(problem.sequence_before),
                constraints_evaluations,
                objectives_evaluations if objectives_evaluations is not None
                else problem.objectives_evaluations(),
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
                           .locations_as_features(label_prefix="Breach from")
            record.features += breaches_locations

            SeqIO.write(record, root._file(title.lower() + ".gb").open("w"),
                        "genbank")

            if breaches_locations != []:
                record.features = breaches_locations
                graphical_record = translator.translate_record(record)
                if len(graphical_record.features) > max_features_in_plots:
                    features = sorted(graphical_record.features,
                                      key= lambda f: f.start - f.end)
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

    jinja_env = jinja2.Environment()
    jinja_env.globals.update(zip=zip, len=len)
    template_path = os.path.join(TEMPLATES_DIR, "optimization_report.html")
    with open(template_path, "r") as f:
        REPORT_TEMPLATE = jinja_env.from_string(f.read())

    html = REPORT_TEMPLATE.render(
        dnachisel_version=__version__,
        project_name="bla",
        problem=problem,
        outcome="SUCCESS" if constraints_evaluations.all_evaluations_pass()
                else "FAIL",
        constraints_after=constraints_evaluations,
        objectives_after=objectives_evaluations,
        edits=sum(len(f) for f in edits)
    )
    weasy_html = weasyprint.HTML(string=html,  base_url=TEMPLATES_DIR)
    weasy_html.write_pdf(root._file("Report.pdf"))


    problem.to_record(root._file("final_sequence.gb").open("w"),
                      with_constraints=False,
                      with_objectives=False)

    # returns zip data if target == '@memory'
    return root._close()
