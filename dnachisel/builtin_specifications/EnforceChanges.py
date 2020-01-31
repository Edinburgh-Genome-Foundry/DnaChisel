"""Implementation of EnforceChanges."""

import numpy as np

from ..Specification import Specification, SpecEvaluation
from ..biotools import (
    sequences_differences_array,
    group_nearby_indices,
    OTHER_BASES,
)
from ..Location import Location

other_bases_sets = {
    base: other_bases for (base, other_bases) in OTHER_BASES.items()
}


class EnforceChanges(Specification):
    """Specify that some locations of the sequence should be changed/different.

    The input for this specification is a bit complex as it is different for
    objectives and constraints:

    - **For constraints**, use EnforceChanges() to force all nucleotides to
      be different from the reference (by default, the reference is the
      starting sequence). Use EnforceChanges(minimum=3) to enforce at least
      3 nucleotide changes, or EnforceChanges(minumum_percent=70) to enforce
      at least 70% change. Note that having a minimum makes computations much
      slower.
    - **For objectives**, use EnforceChanges() to maximize the number of
      nucleotide changes, or use EnforceChanges(amount=10) to aim at a
      10-nucleotides change (no more, no less), or use
      EnforceChanges(minumum_amount=70) to aim at 70% change, no more no less.

    As a genbank annotation:

    - ``@change``, ``@change(minimum=40%)``, ``@change(minimum=3)`` will
      enforce respectivel 100% different nucleotides, 40%+ different, and 3+
      nucleotides different.
    - ``~change``, ``~change(40%)``, ``~change(3)`` to aim at, respectively,
      as close as possible to 100% change, 40% changes, or a target of 3
      nucleotides changes.

    Parameters
    ----------
    location
      Location object indicating the position of the segment that must be
      different from the original sequence. Alternatively,
      indices can be provided. If neither is provided, the assumed location
      is the whole sequence.

    amount
       Number of differences desired in the final sequence when the spec. is
       used as an objective. An amount_percent can be provided instead.

    amount_percent
      Same as ``amount``, except the desired number of differences is decided
      at initialization, as a percentage of the number of nucleotides in the
      region covered.

    minimum
      Minimal number of differences enforced in the final sequence when the
      spec. is used as an constraint. A minimum_percent can be provided
      instead.

    minimum_percent
      Same as ``minimum``, except the desired number of differences is decided
      at initialization, as a percentage of the number of nucleotides in the
      region covered.

    indices
      List of indices that must be different from the original sequence.

    reference
      At the moment, this is rather an internal variable. Do not use unless
      you're not afraid of side effects.

    """

    localization_interval_length = 7  # used when optimizing
    best_possible_score = 0
    shorthand_name = "change"

    def __init__(
        self,
        amount=None,
        amount_percent=None,
        minimum=None,
        minimum_percent=None,
        location=None,
        indices=None,
        reference=None,
        boost=1.0,
    ):
        """Initialize."""
        # raise NotImplementedError("This class is not yet implemented")
        # if location is None and (indices is not None):
        #     location = (min(indices), max(indices) + 1)
        self.location = Location.from_data(location)
        if (self.location is not None) and self.location.strand == -1:
            self.location.strand = 1
        self.indices = np.array(indices) if (indices is not None) else None
        self.reference = reference
        # self.passive_objective = passive_objective
        self.amount = amount
        self.amount_percent = amount_percent
        self.minimum = minimum
        self.minimum_percent = minimum_percent
        if isinstance(amount, str) and amount.endswith("%"):
            self.amount = None
            self.amount_percent = float(amount[:-1])
        if isinstance(minimum, str) and minimum.endswith("%"):
            self.minimum = None
            self.minimum_percent = float(minimum[:-1])
        self.boost = boost

    def extract_subsequence(self, sequence):
        """Extract a subsequence from the location or indices.

        Used to initialize the function when the sequence is provided.

        """
        if (self.location is None) and (self.indices is None):
            return sequence
        elif self.indices is not None:
            return "".join(np.array(list(sequence))[self.indices])
        else:  # self.location is not None:
            return self.location.extract_sequence(sequence)

    def initialized_on_problem(self, problem, role=None):
        """Find out what sequence it is that we are supposed to conserve."""

        result = self._copy_with_full_span_if_no_location(problem)
        L = len(result.location if result.indices is None else result.indices)

        # FIRST DEAL WITH THE AMOUNTS OF CHANGE

        both_amounts_none = (self.amount is None) and (
            self.amount_percent is None
        )
        both_minimums_none = (self.minimum is None) and (
            self.minimum_percent is None
        )

        # Only if minimum_percent=100 is the constraint nucleotide-enforced
        # So this variable is updated below if necessary
        self.enforced_by_nucleotide_restrictions = False

        if role == "constraint":
            if not both_amounts_none:
                raise ValueError(
                    "Use minimum(_percent) parameter when using EnforceChanges"
                    " as a constraint"
                )
            if both_minimums_none:
                result.minimum_percent = 100
            if result.minimum_percent is not None:
                if result.minimum_percent == 100:
                    self.enforced_by_nucleotide_restrictions = 100
                result.minimum = np.ceil(result.minimum_percent * L / 100.0)

        elif role == "objective":
            if not both_minimums_none:
                raise ValueError(
                    "Use amount(_percent) parameter when using EnforceChanges "
                    " as an objective"
                )
            if both_amounts_none:
                result.amount_percent = 100
            if result.amount_percent is not None:
                result.amount = result.amount_percent * L / 100.0

        # THEN DEAL WITH THE REFERENCE

        if (self.reference is None) and (self.indices is not None):
            result = result.copy_with_changes()
            result.reference = self.extract_subsequence(problem.sequence)
        elif self.reference is None or (
            (self.indices is None) and (len(self.reference) < L)
        ):
            result = result.copy_with_changes()
            result.reference = self.extract_subsequence(problem.sequence)
        return result

    def evaluate(self, problem):
        """Return a score equal to -number_of_equalities.

        Locations are "binned" equality regions. Each bin has a length
        in nucleotides equal to ``localization_interval_length`.`
        """

        # FIND THE INDICES WHERE THE SEQUENCE IS UNCHANGED

        # Note: at this stage any minimum_percent or amount_percent have been
        # transformed into abolsute self.minimum and self.amount.

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

        def indices_to_intervals(indices):
            intervals = group_nearby_indices(
                indices, max_group_spread=self.localization_interval_length
            )
            return [(interval[0], interval[-1] + 1) for interval in intervals]

        if self.indices is not None:
            n_indices = len(self.indices)
        else:
            n_indices = len(self.location)
        n_differences = n_indices - len(equalities)
        if self.minimum is not None:
            score = n_differences - self.minimum
            intervals = indices_to_intervals(equalities)
        else:
            score = -abs(n_differences - self.amount)
            if n_differences <= self.amount:
                intervals = indices_to_intervals(equalities)
            else:
                differences = [
                    i for i in self.location.indices if i not in equalities
                ]
                intervals = indices_to_intervals(differences)
        locations = (
            [self.location]
            if (self.minimum is not None)
            else [Location(start, end, 1) for start, end in intervals]
        )
        return SpecEvaluation(self, problem, score=score, locations=locations)

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the spec to the overlap of its location and the new.
        """

        # If minimum_percent=100: everything is covered by mutation restriction
        # If minimum_percent<100: no smart strategy return self.
        # If amount_percent<100: no smart strategy return self.
        # If amount_percent==100: easy to localize, it is always locally 100%

        if self.amount_percent != 100:
            return self

        # Now the only case where localization makes sense: amount_percent=100

        start, end = location.start, location.end
        if self.indices is not None:
            pos = ((start <= self.indices) & (self.indices < end)).nonzero()[0]
            new_indices = self.indices[pos]
            new_reference = "".join(np.array(list(self.reference))[pos])
            return self.copy_with_changes(
                indices=new_indices,
                reference=new_reference,
                minimum=None if self.minimum is None else len(new_indices),
                amount=None if self.amount is None else len(new_indices),
            )
        else:
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return None
            else:
                new_constraint = self.copy_with_changes(
                    location=new_location,
                    minimum=None if self.minimum is None else len(location),
                    amount=None if self.amount is None else len(location),
                )
                relative_location = new_location + (-self.location.start)
                new_constraint.reference = relative_location.extract_sequence(
                    self.reference
                )
                return new_constraint

    def restrict_nucleotides(self, sequence, location=None):
        """When localizing, forbid any nucleotide but the one already there."""

        # IF IT IS NOT AN ALL-DIFFERENT CONSTRAINT, DO NOTHING

        if self.minimum_percent != 100:
            return []

        # ELSE, IMPOSE THAT EACH NUCL. IS ONE OF THE 3 NON-INITIAL ONES

        if location is not None:
            start = max(location.start, self.location.start)
            end = min(location.end, self.location.end)
        else:
            start, end = self.location.start, self.location.end

        if self.indices is not None:
            return [
                ((i, i + 1), other_bases_sets[sequence[i : i + 1]])
                for i in self.indices
                if start <= i < end
            ]
        else:
            return [
                ((i, i + 1), other_bases_sets[sequence[i : i + 1]])
                for i in range(start, end)
            ]

    def short_label(self):
        return "change"
    
    def breach_label(self):
        return "unchanged"
