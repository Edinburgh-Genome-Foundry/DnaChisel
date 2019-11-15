"""Implementation of AvoidMatches."""

from ..Specification import Specification, SpecEvaluation
from ..biotools.bowtie import (
    find_all_bowtie_matches,
    create_bowtie_index_from_sequences,
)
from ..Location import Location
import tempfile
import shutil


class AvoidMatches(Specification):
    """Enforce that the given pattern is absent in the sequence.

    Uses Bowtie.

    Parameters
    ----------

    bowtie_index
      Path to a local bowtie2 index.

    match_length
      Length of the matches to avoid

    mismatches
      Number of single-nucleotide mismatches allowed inside each match. Only
      0-3 mismatches is supported (this is what Bowtie supports.)

    location
      Location of the sequence on which the specification applies
    """

    priority = -2
    best_possible_score = 0
    blasts_paths = {}
    shorthand_name = "no_match"
    bowtie_indexes_paths = None

    def __init__(
        self,
        match_length,
        bowtie_index=None,
        sequences=None,
        mismatches=0,
        location=None,
        boost=1,
    ):
        """Initialize."""
        self.bowtie_index = bowtie_index
        self.match_length = match_length
        self.mismatches = mismatches
        self.sequences = sequences
        self.location = Location.from_data(location)
        self.boost = boost
        self._tmp_data_dir = None

        if self.sequences is not None:
            self._tmp_data_dir = tempfile.mkdtemp(prefix="chisel_AvoidMatch")
            self.bowtie_index = create_bowtie_index_from_sequences(
                sequences=sequences, path=self._tmp_data_dir
            )
        self.is_localized = False

    def initialized_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        """Find matches and count them negatively"""
        location = self.location
        if location is None:
            location = Location(0, len(problem.sequence))
        sequence = location.extract_sequence(problem.sequence)
        bowtie_index = self.bowtie_index
        if self.bowtie_indexes_paths is not None:
            if self.bowtie_index in self.bowtie_indexes_dict:
                bowtie_index = self.bowtie_indexes_paths[self.bowtie_index]
        # matches will be of the form [ ((start, end), n_mismatches), (...)]
        matches = find_all_bowtie_matches(
            sequence=sequence,
            bowtie_index_path=bowtie_index,
            match_length=self.match_length,
            max_mismatches=self.mismatches,
        )
        # the score is designed to be -1 per exact match, and improve when
        # there are mismatches more mismatches,
        score = -sum([1.0 / (1 + mis) for _, mis in matches])
        locations = [Location(start, end) for (start, end), _ in matches]

        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=locations,
            message="Failed - %s matches at %s" % (len(locations), locations),
        )

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the evaluation."""
        new_location = self.location.overlap_region(location)
        if new_location is None:
            return None
        new_location = location.extended(
            self.match_length - 1, right=with_righthand
        )
        return self.copy_with_changes(location=new_location, is_localized=True)

    def feature_label_parameters(self):
        return [self.blast_db]
    
    def remove_temp_directory(self):
        if self._tmp_data_dir is not None:
            shutil.rmtree(self._tmp_data_dir)