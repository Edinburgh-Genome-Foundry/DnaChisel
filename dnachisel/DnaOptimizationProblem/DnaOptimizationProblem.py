"""Define the DnaOptimizationProblem class.

DnaOptimizationProblem is where the whole problem is defined: sequence,
constraints, objectives.
"""

from Bio.SeqRecord import SeqRecord
from proglog import default_bar_logger
from ..Specification.SpecificationSet import SpecificationSet
from ..biotools import sequences_differences_array
from ..MutationSpace import MutationSpace
from ..reports.optimization_reports import (
    write_optimization_report,
    write_no_solution_report,
)
from .NoSolutionError import NoSolutionError
from . import mixins


class DnaOptimizationProblem(
    mixins.ConstraintsSolverMixin,
    mixins.ObjectivesMaximizerMixin,
    mixins.RecordRepresentationMixin,
):
    """Problem specifications: sequence, constraints, optimization objectives.

    The original constraints, objectives, and original sequence of the problem
    are stored in the DNA problem. This class also has methods to display
    reports on the constraints and objectives, as well as solving the
    constraints and objectives.

    Examples
    --------

    >>> from dnachisel import *
    >>> problem = DnaOptimizationProblem(
    >>>     sequence = "ATGCGTGTGTGC...",
    >>>     constraints = [constraint1, constraint2, ...],
    >>>     objectives = [objective1, objective2, ...]
    >>> )
    >>> problem.resolve_constraints()
    >>> problem.optimize()
    >>> print(problem.constraints_text_summary())
    >>> print(problem.objectives_text_summary())


    Parameters
    ----------

    sequence
      A string of ATGC characters (they must be upper case!), e.g. "ATTGTGTA"

    constraints
      A list of objects of type ``Specification``.

    objectives
      A list of objects of type ``Specification`` specifying what must be
      optimized in the problem. Note that each objective has a float ``boost``
      parameter. The larger the boost, the more the objective is taken into
      account during the optimization.

    logger
      Either None for no logger, 'bar' for a tqdm progress bar logger, or
      any ProgLog progress bar logger.

    mutations_space
      A MutationSpace indicating the possible mutations. In most case the
      mutation space will be left to None and computed at problem
      initialization (which can be slightly compute-intensive), however some
      core DNA Chisel methods will create optimization problems with a provided
      mutation_space to save computing time.

    Attributes
    ----------

    randomization_threshold
      The algorithm will use an exhaustive search when the size of the mutation
      space (=the number of possible variants) is above this threshold, and
      a (guided) random search when it is above.

    max_random_iters
      When using a random search, stop after this many iterations

    mutations_per_iteration
      When using a random search, produce this many sequence mutations each
      iteration.

    optimization_stagnation_tolerance
      When using a random search, stop if the score hasn't improved in
      the last "this many" iterations

    local_extensions
      Try local resolution several times if it fails, increasing the mutable
      zone by [N1, N2...] nucleotides on each side, until resolution works.
      (by default, an extension of 0bp is tried, then 5bp.

    Notes
    -----

    The dictionary ``self.possible_mutations`` is of the form
    ``{location1 : list1, location2: list2...}``
    where ``location`` is either a single index (e.g. 10) indicating the
    position of a nucleotide to be muted, or a couple ``(start, end)``
    indicating a whole segment whose sub-sequence should be replaced.
    The ``list`` s are lists of possible sequences to replace each location,
    e.g. for the mutation of a whole codon ``(3,6): ["ATT", "ACT", "AGT"]``.
    """

    randomization_threshold = 10000
    max_random_iters = 1000
    mutations_per_iteration = 2
    optimization_stagnation_tolerance = 100
    local_extensions = (0, 5)

    def __init__(
        self,
        sequence,
        constraints=None,
        objectives=None,
        logger="bar",
        mutation_space=None,
    ):
        """Initialize"""
        if isinstance(sequence, SeqRecord):
            self.record = sequence
            self.sequence = str(sequence.seq).upper()
        else:
            self.record = None
            self.sequence = sequence.upper()
        self.constraints = [] if constraints is None else list(constraints)
        self.objectives = [] if objectives is None else list(objectives)
        self.logger = default_bar_logger(
            logger,
            bars=("objective", "constraint", "location"),
            ignored_bars=("mutation",),
            min_time_interval=0.2,
        )
        self.mutation_space = mutation_space
        self.initialize()

    def initialize(self):
        """Precompute specification sets, evaluations, and mutation space."""

        # Find the specifications (objectives, constraints) which are actually
        #  SpecificationSets, and unpack these to complete the lists of
        # objectives and constraints.
        for specs in (self.constraints, self.objectives):
            specsets = [
                spec for spec in specs if isinstance(spec, SpecificationSet)
            ]
            specs_in_sets = [
                spec
                for specset in specsets
                for spec in specset.specifications.values()
            ]
            for specset in specsets:
                specs.remove(specset)
            specs.extend(specs_in_sets)

        # INITIALIZE THE CONSTRAINTS AND OBJECTIVES

        self.constraints = [
            constraint.initialized_on_problem(self, role="constraint")
            for constraint in self.constraints
        ]
        self.objectives = [
            objective.initialized_on_problem(self, role="objective")
            for objective in self.objectives
        ]

        # INITIALIZE THE "BEFORE" CLASS ATTRIBUTES, USED IN REPORTS

        self.sequence_before = self.sequence
        self._constraints_before = None
        self._objectives_before = None

        # INITIALIZE THE MUTATION SPACE

        if self.mutation_space is None:
            self.mutation_space = MutationSpace.from_optimization_problem(self)
            # If the original sequence is outside of the allowed mutations
            # space, replace the sequence by a sequence which complies with
            # the mutation space.
            self.sequence = self.mutation_space.constrain_sequence(
                self.sequence
            )

    def _replace_sequence(self, new_sequence):
        """Replace the current sequence of the problem.

        This method is subclassed in CircularDnaOptimization problem where
        is is more complex (changing the sequence in one location changes
        it in more locations).
        """
        self.sequence = new_sequence

    def sequence_edits_as_array(self):
        """Return an array [False, False, True...] where True indicates an edit
        (i.e. a change at this position between the original problem sequence
        and the current one)."""
        return sequences_differences_array(self.sequence, self.sequence_before)

    def number_of_edits(self):
        """Return the number of nucleotide differences between the original
        and current sequence."""
        return self.sequence_edits_as_array().sum()

    def optimize_with_report(
        self,
        target,
        project_name="My project",
        file_path=None,
        file_content=None,
    ):
        """Resolve constraints, optimize objectives, write a multi-file report.

        The report's content may vary depending on the optimization's success.

        Parameters
        ----------

        target
          Either a path to a folder that will contain the report, or a path to
          a zip archive, or "@memory" to return raw data of a zip archive
          containing the report.

        project_name
          Project name to write on PDF reports

        Returns
        -------

        (success, message, zip_data)
          Triplet where success is True/False, message is a one-line string
          summary indication whether some clash was found, or some solution, or
          maybe no solution was found because the random searches were too short
        """
        self.logger(message="Solving constraints")
        try:
            self.resolve_constraints()
        except NoSolutionError as error:
            self.logger(message="No solution found: making report")
            data = write_no_solution_report(
                target,
                self,
                error,
                file_path=file_path,
                file_content=file_content,
            )
            start, end, s = error.location.to_tuple()
            message = "No solution found in zone [%d, %d]: %s." % (
                start,
                end,
                str(error),
            )
            return False, message, data
        self.logger(message="Now optimizing the sequence")
        self.optimize()
        self.logger(message="Success! Generating report.")
        data = write_optimization_report(
            target,
            self,
            project_name=project_name,
            file_path=file_path,
            file_content=file_content,
        )
        return True, "Optimization successful.", data
