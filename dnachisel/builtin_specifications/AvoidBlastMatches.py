"""Implementation of AvoidBlastMatches."""

from ..Specification import Specification
from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from dnachisel.biotools import blast_sequence, group_nearby_segments
from dnachisel.Location import Location

class AvoidBlastMatches(Specification):
    """Enforce that the given pattern is absent in the sequence.

    Uses NCBI Blast+. Only local BLAST is supported/tested as for now

    Parameters
    ----------

    blast_db
      Path to a local BLAST database. These databases can be obtained with
      NCBI's `makeblastdb`. Omit the extension, e.g. `ecoli_db/ecoli_db`.

    word_size
      Word size used by the BLAST algorithm

    perc_identity
      Minimal percentage of identity for BLAST matches. 100 means that only
      perfect matches are considered.

    num_alignments
      Number alignments

    num_threads
      Number of threads/CPU cores to use for the BLAST algorithm.

    min_align_length
      Minimal length that an alignment should have to be considered.
    """
    priority = -2
    best_possible_score = 0

    def __init__(self, blast_db=None, sequences=None, word_size=4,
                 perc_identity=100, num_alignments=100000, num_threads=3,
                 min_align_length=20, ungapped=True, e_value=1e80,
                 culling_limit=1, location=None):
        """Initialize."""
        if isinstance(location, tuple):
            location = Location.from_tuple(location)
        self.blast_db = blast_db
        self.sequences = sequences
        self.word_size = word_size
        self.perc_identity = perc_identity
        self.num_alignments = num_alignments
        self.num_threads = num_threads
        self.min_align_length = min_align_length
        self.location = location
        self.e_value = e_value
        self.ungapped = ungapped
        self.culling_limit = culling_limit

    def initialize_on_problem(self, problem, role):
        """Find out what sequence it is that we are supposed to conserve."""
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            return self.copy_with_changes(location=location)
        else:
            return self

    def evaluate(self, problem):
        """Score as (-total number of blast identities in matches)."""
        location = self.location
        if location is None:
            location = Location(0, len(problem.sequence))
        sequence = location.extract_sequence(problem.sequence)

        blast_record = blast_sequence(
            sequence, blast_db=self.blast_db,
            subject_sequences=self.sequences,
            word_size=self.word_size,
            perc_identity=self.perc_identity,
            num_alignments=self.num_alignments,
            num_threads=self.num_threads,
            ungapped=self.ungapped,
            e_value=self.e_value,
            culling_limit=self.culling_limit
        )

        if isinstance(blast_record, list):
            alignments = [
                alignment
                for rec in blast_record
                for alignment in rec.alignments
            ]
        else:
            alignments = blast_record.alignments

        query_hits = [
            (
                min(hit.query_start, hit.query_end) + location.start - 1,
                max(hit.query_start, hit.query_end) + location.start,
                1 - 2 * (hit.query_start > hit.query_end),
                hit.identities
            )
            for alignment in alignments
            for hit in alignment.hsps
        ]

        locations = sorted([
            (start, end, ids)
            for (start, end, strand, ids) in query_hits
            if (end - start) >= self.min_align_length
        ])
        # locations = [
        #     (r[0][0], r[-1][-1])
        #     for r in group_nearby_segments(locations, max_start_spread=2)
        # ]

        score = -sum([ids for start, end, ids in locations])
        locations = [Location(start, end) for start, end, ids in locations]

        if locations == []:
            return SpecEvaluation(self, problem, score=1,
                                  message="Passed: no BLAST match found")



        return SpecEvaluation(
            self, problem, score=score, locations=locations,
            message="Failed - matches at %s" % locations)

    def localized(self, location, problem=None, with_righthand=True):
        """Localize the evaluation."""
        new_location = self.location.overlap_region(location)
        if new_location is None:
            return VoidSpecification(parent_specification=self)

        new_location = location.extended(self.min_align_length - 1,
                                         right=with_righthand)
        return self.copy_with_changes(location=new_location)

    def feature_label_parameters(self):
        return [self.blast_db]
