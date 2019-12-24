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

from ...biotools import sequence_to_biopython_record, annotate_record
from ...builtin_specifications import (
    EnforceGCContent,
    AvoidPattern,
    AvoidHairpins,
)
from ..colors_cycle import colors_cycle

from .GraphicTranslator import GraphicTranslator

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    MPL_AVAILABLE = True
except ImportError:
    MPL_AVAILABLE = False

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


def plot_breaches_record(record, ax=None, figure_width=10):
    translator = GraphicTranslator()
    graphic_record = translator.translate_record(record)
    ax, _ = graphic_record.plot(
        ax=ax, figure_width=figure_width, strand_in_label_threshold=7
    )
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
