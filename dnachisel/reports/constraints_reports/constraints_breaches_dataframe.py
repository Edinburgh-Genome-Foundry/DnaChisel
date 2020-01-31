try:
    import pandas

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

from ...DnaOptimizationProblem import DnaOptimizationProblem
from ...Location import Location

def _breaches(constraint, sequence):
    problem = DnaOptimizationProblem(sequence, mutation_space={})
    new_constraint = constraint.initialized_on_problem(problem, role=None)
    evaluation = new_constraint.evaluate(problem)
    locations = Location.merge_overlapping_locations(evaluation.locations)
    return ", ".join(map(str, locations))

def _install_extras_message(libname):
    return (
        "Could not load %s (is it installed ?). You can install it separately "
        " with:  pip install %s\n\n"
        "Install all dependencies for generating DNA Chisel reports with:"
        "\n\npip install dnachisel[reports]"
        % (libname, libname.lower().replace(" ", "_"))
    )

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
                        use_breach_form=True,
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