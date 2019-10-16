"""Implementation of EnforceChanges."""

import numpy as np

from ..Specification import Specification
from ..SpecEvaluation import SpecEvaluation

# from .VoidSpecification import VoidSpecification
from dnachisel.biotools import (
    sequences_differences_array,
    group_nearby_indices,
    OTHER_BASES,
)

other_bases_sets = {
    base: other_bases for (base, other_bases) in OTHER_BASES.items()
}


from dnachisel.Location import Location


class EnforceChanges(Specification):
    """Specify that some locations of the sequence should be changed/different.

    Note: for now this class is simply a derivative of AvoidChanges where
    the scoring function penalizes equalities instead of differences.

    Parameters
    ----------
    location
      Location object indicating the position of the segment that must be
      different from the original sequence. Alternatively,
      indices can be provided. If neither is provided, the assumed location
      is the whole sequence.

    indices
      List of indices that must be different from the original sequence.

    reference
      At the moment, this is rather an internal variable. Do not use unless
      you're not afraid of side effects.

    """

    localization_interval_length = 7  # used when optimizing
    best_possible_score = 0
    enforced_by_nucleotide_restrictions = True

    def __init__(self, location=None, indices=None, reference=None, boost=1.0):
        """Initialize."""
        # raise NotImplementedError("This class is not yet implemented")
        if isinstance(location, tuple):
            location = Location.from_tuple(location)
        self.location = location
        if (self.location is not None) and self.location.strand == -1:
            self.location.strand = 1
        self.indices = np.array(indices) if (indices is not None) else None
        self.reference = reference
        # self.passive_objective = passive_objective
        self.boost = boost

    def extract_subsequence(self, sequence):
        """Extract a subsequence from the location or indices.

        Used to initialize the function when the sequence is provided.

        """
        if (self.location is None) and (self.indices is None):
            return sequence
        elif self.indices is not None:
            return "".join(np.array(sequence)[self.indices])
        else:  # self.location is not None:
            return self.location.extract_sequence(sequence)

    def initialized_on_problem(self, problem, role=None):
        """Find out what sequence it is that we are supposed to conserve."""
        result = self._copy_with_full_span_if_no_location(problem)

        # Initialize the "reference" in two cases:
        # - Always at the very beginning
        # - When the new sequence is bigger than the previous one
        #   (used in CircularDnaOptimizationProblem)
        if self.reference is None or (
            len(self.reference) < len(self.location)
        ):
            result = result.copy_with_changes()
            result.reference = self.extract_subsequence(problem.sequence)
        return result

    def evaluate(self, problem):
        """Return a score equal to -number_of_equalities.

        Locations are "binned" equality regions. Each bin has a length
        in nucleotides equal to ``localization_interval_length`.`
        """
        target = self.reference
        sequence = self.extract_subsequence(problem.sequence)
        equalities = np.nonzero(
            1 - sequences_differences_array(sequence, target)
        )[0]

        if self.indices is not None:
            equalities = self.indices[equalities]
        elif self.location is not None:
            if self.location.strand == -1:
                equalities = self.location.end - equalities
            else:
                equalities = equalities + self.location.start

        intervals = [
            (r[0], r[-1] + 1)
            for r in group_nearby_indices(
                equalities,
                max_group_spread=self.localization_interval_length,
            )
        ]
        locations = [Location(start, end, 1) for start, end in intervals]

        return SpecEvaluation(
            self, problem, score=-len(equalities), locations=locations
        )

    def localized(self, location, problem=None):
        """Localize the spec to the overlap of its location and the new.
        """
        start, end = location.start, location.end
        if self.location is not None:
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return None
            # VoidSpecification(parent_specification=self)
            else:
                # return self
                # TODO: refine using the code hereunder, which sometimes
                # creates exceptions like "different sequences"

                new_constraint = self.copy_with_changes(location=new_location)
                relative_location = new_location + (-self.location.start)
                new_constraint.reference = relative_location.extract_sequence(
                    self.reference
                )
                return new_constraint

        elif self.indices is not None:
            inds = self.indices
            new_indices = inds[(start <= inds) & (inds <= end)]
            return self.copy_with_changes(indices=new_indices)
        else:
            return self

    def restrict_nucleotides(self, sequence, location=None):
        """When localizing, forbid any nucleotide but the one already there."""
        if location is not None:
            start = max(location.start, self.location.start)
            end = min(location.end, self.location.end)
        else:
            start, end = self.location.start, self.location.end
        return [
            ((i, i + 1), other_bases_sets[sequence[i : i + 1]])
            for i in range(start, end)
        ]
        # return [(i, set(sequence[i])) for i in range(start, end)]

    def short_label(self):
        return "change"

