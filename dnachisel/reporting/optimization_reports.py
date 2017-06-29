"""Methods to generate reports from optimizations"""

import os
import textwrap

from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import jinja2
import weasyprint
import flametree

from ..biotools import sequence_to_biopython_record, crop_record
from ..plotting_tools import ObjectivesAnnotationsTranslator
from ..version import __version__
from ..DnaOptimizationProblem import NoSolutionError

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEMPLATES_DIR = os.path.join(THIS_DIR, "templates")
TITLE_FONTDICT = fontdict = dict(size=14, weight="bold")

def optimization_with_report(problem, target, project_name="unnamed",
                             **solver_parameters):
    """

    Returns
    -------

    success, message, zip_dat
    """
    try:
        problem.solve_all_constraints_one_by_one(**solver_parameters)
    except NoSolutionError as error:
        data = write_no_solution_report(target, problem, error)
        start, end, s = error.location.as_tuple()
        message = ("No solution found in zone [%d, %d]: %s" %
                   (start, end, str(error)))
        return False, message, data
    constraints_evaluations = problem.constraints_evaluations()
    if constraints_evaluations.all_evaluations_pass():
        problem.maximize_all_objectives_one_by_one()
    data = write_optimization_report(
        target, problem, project_name=project_name,
        constraints_evaluations=constraints_evaluations)
    return True, "Optimization successful", data





def write_no_solution_report(target, problem, error):
    root = flametree.file_tree(target, replace=True)
    translator = ObjectivesAnnotationsTranslator()
    with PdfPages(root._file("plots.pdf").open("wb")) as pdf_io:

        # PLOT GLOBAL LOCATION OF ERROR

        record = problem.to_record()
        translator = ObjectivesAnnotationsTranslator()
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


def write_optimization_report(target, problem, project_name="unnammed",
                              constraints_evaluations=None,
                              objectives_evaluations=None,
                              figure_width=20):
    """Write a report with a PDF summary, plots, and genbanks."""
    if constraints_evaluations is None:
        constraints_evaluations = problem.constraints_evaluations()
    if objectives_evaluations is None:
        objectives_evaluations = problem.objectives_evaluations()
    root = flametree.file_tree(target, replace=True)
    translator = ObjectivesAnnotationsTranslator()


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
        objectives_after=objectives_evaluations
    )
    weasy_html = weasyprint.HTML(string=html,  base_url=TEMPLATES_DIR)
    weasy_html.write_pdf(root._file("Report.pdf"))

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
                "{title}:\t\t{nfailing} constraints failing"
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
                ax, _ = graphical_record.plot(figure_width=figure_width)
                ax.set_title(title + ": Constraints breaches locations",
                             loc="left", fontdict=TITLE_FONTDICT)
                height = max(ax.figure.get_size_inches()[1], plot_height)
                ax.figure.set_size_inches((figure_width, height))
                pdf_io.savefig(ax.figure, bbox_inches="tight")
                plt.close(ax.figure)

    problem.to_record(root._file("final_sequence.gb").open("w"),
                      with_constraints=False,
                      with_objectives=False)

    return root._close()
