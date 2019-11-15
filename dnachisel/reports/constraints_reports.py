"""Misc. plotting and reporting methods, some of which are really arbitrary.


Here is a typical example of use:

>>> import dnachisel.reports.constraint_reports as cr
>>> dataframe = cr.constraints_breaches_dataframe(constraints, sequences)
>>> dataframe.to_excel("output_breaches.xlsx")
>>> records = cr.records_from_breaches_dataframe(dataframe, sequences)
>>> cr.breaches_records_to_pdf(records, 'output_breaches_plots.pdf')
"""

from copy import deepcopy
import re
from io import BytesIO

import proglog

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from dna_features_viewer import BiopythonTranslator

    DFV_AVAILABLE = True
except ImportError:
    BiopythonTranslator = object
    DFV_AVAILABLE = False

try:
    import pandas

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

from ..DnaOptimizationProblem import DnaOptimizationProblem
from ..biotools import sequence_to_biopython_record, annotate_record
from ..builtin_specifications import (
    EnforceGCContent,
    AvoidPattern,
    AvoidHairpins,
)
from ..Location import Location
from .colors_cycle import colors_cycle


__all__ = [
    "records_from_breaches_dataframe",
    "plot_breaches_record",
    "breaches_records_to_pdf",
    "constraints_breaches_dataframe",
]


def _sequences_to_new_records(sequences):
    """Turn acceptable sequences input into a records list.

    Acceptable formats are
    - ('name', 'sequence')
    - {'name': 'sequence'}
    - [records] (will be deepcopied)
    """
    if isinstance(sequences, dict):
        sequences = list(sequences.items())
    records = []
    for seq in sequences:
        if hasattr(seq, "id"):
            records.append(deepcopy(seq))
        else:
            name, seq = seq
            records.append(
                sequence_to_biopython_record(seq, id=name, name=name)
            )
    return records


def _parse_location(location_string):
    """Parses locations like 235-350(+)"""
    location_regex = r"(\d+)-(\d+)(\(+\)|\(-\)|)"
    match = re.match(location_regex, location_string.strip())
    start, end, strand = match.groups()
    return int(start), int(end), -1 if strand == "(-)" else 1


def _install_extras_message(libname):
    return (
        "Could not load %s (is it installed ?). You can install it separately "
        " with:  pip install %s\n\n"
        "Install all dependencies for generating DNA Chisel reports with:"
        "\n\npip install dnachisel[reports]"
        % (libname, libname.lower().replace(" ", "_"))
    )


def _breaches(constraint, sequence):
    problem = DnaOptimizationProblem(sequence, mutation_space={})
    new_constraint = constraint.initialized_on_problem(problem, role=None)
    evaluation = new_constraint.evaluate(problem)
    locations = Location.merge_overlapping_locations(evaluation.locations)
    return ", ".join(map(str, locations))


def constraints_breaches_dataframe(
    constraints, sequences, display_constraints_locations=False,
):
    """Return a dataframe summarizing constraints breaches in the sequences.

    Output dataframe schema (cst = constraint):

    =====  ========  ======================
      /    Cst1       Cst2
    =====  ========  ======================
    Seq1   10-50(+)  100-200(+), 300-350(+)
    seq2   
    Seq3   2-10(+)
    Seq4             500-1000(-)
    =====  ========  ======================


    Parameters
    ----------

    constraints
      A list of DNA Chisel Specifications.

    sequences
      Either a list [("name", "sequence")...] or a dict {"name": "sequence"}
      or a list of biopython records whole id is the sequence name.

    Examples
    --------

    >>> import dnachisel as dc
    >>> from dnachisel.utils import constraints_breaches_dataframe
    >>> sequences = [
    >>>     ("SEQ 1", "ATTGTGCATGTGACAC"),
    >>>     ("SEQ 2", "ACATGTGTTGTGACAC"),
    >>>     ("SEQ 3", "TTGTGCACACATGTGA"),
    >>> ]
    >>> constraints = [
    >>>     dc.AvoidPattern('ATTG'),
    >>>     dc.EnforceGCContent(0.4, 0.6),
    >>>     dc.UniquifyAllKmers(5)
    >>> ]
    >>> dataframe = constraints_breaches_dataframe(constraints, sequences)
    >>> dataframe.to_excel('summary_spreadsheet.xlsx')
    """
    if not PANDAS_AVAILABLE:
        raise ImportError(_install_extras_message("Pandas"))
    if isinstance(sequences, dict):
        sequences = list(sequences.items())
    if hasattr(sequences[0], "id"):
        sequences = [(s.id, s) for s in sequences]

    dataframe_records = [
        dict(
            [("sequence", name)]
            + [
                (
                    constraint.label(
                        use_short_form=True,
                        with_location=display_constraints_locations,
                    ),
                    _breaches(constraint, sequence),
                )
                for constraint in constraints
            ]
        )
        for (name, sequence) in sequences
    ]

    return pandas.DataFrame.from_records(dataframe_records, index="sequence")


def records_from_breaches_dataframe(dataframe, sequences):
    """Generate records with annotations indicating constraints breaches.
    
    Parameters
    ----------

    dataframe
      A breaches dataframe returned by ``constraints_breaches_dataframe``
    
    sequences
      Either a list [("name", "sequence")...] or a dict {"name": "sequence"}
      or a list of biopython records whole id is the sequence name.

    """
    records = _sequences_to_new_records(sequences)
    for record in records:
        record.features = [
            f
            for f in record.features
            if not f.qualifiers.get("is_a_breach", False)
        ]
    colors_cycle_iterator = colors_cycle()
    columns_colors = {
        c: next(colors_cycle_iterator) for c in dataframe.columns
    }
    for rec, (i, row) in zip(records, dataframe.iterrows()):
        for column in dataframe.columns:
            locations = row[column]
            if not locations:
                continue
            for location in locations.split(","):
                annotate_record(
                    rec,
                    location=_parse_location(location),
                    label=column,
                    color=columns_colors[column],
                    ApEinfo_fwdcolor=columns_colors[column],
                    ApEinfo_revcolor=columns_colors[column],
                    is_a_breach=True,
                )
    return records


class Translator(BiopythonTranslator):
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


def plot_breaches_record(record, ax=None, figure_width=10):
    translator = Translator()
    graphic_record = translator.translate_record(record)
    ax, _ = graphic_record.plot(ax=ax, figure_width=figure_width)
    ax.set_title(record.id, loc="left", fontweight="bold")
    ax.set_ylim(top=ax.get_ylim()[1] + 1)
    return ax


def breaches_records_to_pdf(
    breaches_records, pdf_path=None, figure_width=10, logger="bar"
):
    """Plots figures of the breaches annotated in the records into a PDF file.
    
    Parameters
    ----------

    breaches_records
      A least of records annotated with breaches, as returned by the
    
    pdf_path
      Either the path to a PDF, or a file handle (open in wb mode) or None
      for this method to return binary PDF data.
    
    logger
      Either "bar" for a progress bar, None for no logging, or any Proglog
      logger. The bar name is "sequence".
    """
    pdf_io = BytesIO() if pdf_path is None else pdf_path
    logger = proglog.default_bar_logger(logger, min_time_interval=0.2)

    with PdfPages(pdf_io) as pdf:
        for record in logger.iter_bar(sequence=breaches_records):
            ax = plot_breaches_record(record, figure_width=figure_width)
            pdf.savefig(ax.figure, bbox_inches="tight")
            plt.close(ax.figure)
    if pdf_path is None:
        return pdf_io.getvalue()


EXAMPLE_MANUFACTURING_CONSTRAINTS = [
    AvoidPattern("BsaI_site"),
    AvoidPattern("BsmBI_site"),
    AvoidPattern("BbsI_site"),
    AvoidPattern("SapI_site"),
    AvoidPattern("9xA"),
    AvoidPattern("9xT"),
    AvoidPattern("6xG"),
    AvoidPattern("6xC"),
    AvoidPattern("5x3mer"),
    AvoidPattern("9x2mer"),
    AvoidHairpins(stem_size=20, hairpin_window=200),
    EnforceGCContent(mini=0.3, maxi=0.7, window=100),
]
