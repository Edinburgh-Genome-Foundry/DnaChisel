"""Collection of useful pre-defined objectives/constraints for DnaChisel."""
from collections import Counter, defaultdict
import itertools

import numpy as np

from ..biotools.biotables import (CODON_USAGE_TABLES, CODONS_SEQUENCES,
                                  CODONS_TRANSLATIONS)
from ..biotools.biotools import (gc_content, reverse_complement,
                                 blast_sequence, translate)
from .Objective import (Objective, PatternObjective, TerminalObjective,
                        ObjectiveEvaluation, VoidObjective)
from ..Location import Location


class AvoidBlastMatches(Objective):
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

    def __init__(self, blast_db, word_size=4, perc_identity=100,
                 num_alignments=1000, num_threads=3, min_align_length=20,
                 location=None):
        self.blast_db = blast_db
        self.word_size = word_size
        self.perc_identity = perc_identity
        self.num_alignments = num_alignments
        self.num_threads = num_threads
        self.min_align_length = min_align_length
        self.location = location

    def evaluate(self, problem):
        """Return (-M) as a score, where M is the number of BLAST matches found
        in the BLAST database."""
        location = self.location
        if location is None:
            location = Location(0, len(problem.sequence))
        sequence = location.extract_sequence(problem.sequence)
        blast_record = blast_sequence(
            sequence, blast_db=self.blast_db,
            word_size=self.word_size,
            perc_identity=self.perc_identity,
            num_alignments=self.num_alignments,
            num_threads=self.num_threads
        )
        query_locations = [
            Location(min(hit.query_start, hit.query_end),
                     max(hit.query_start, hit.query_end),
                     1 - 2 * (hit.query_start > hit.query_end))
            for alignment in blast_record.alignments
            for hit in alignment.hsps
        ]
        locations = sorted([
            loc for loc in query_locations
            if len(location) >= self.min_align_length
        ])
        if locations == []:
            return ObjectiveEvaluation(self, problem, score=1,
                                       message="Passed: no BLAST match found")

        return ObjectiveEvaluation(
            self, problem, score=-len(locations), locations=locations,
            message="Failed - matches at %s" % locations)

    def localized(self, location):
        """Localize the evaluation."""
        if self.location is not None:
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return VoidObjective(parent_objective=self)
        else:
            new_location = location.extended(self.min_align_length)

        return self.copy_with_changes(location=new_location)

    def __repr__(self):
        return "NoBlastMatchesObjective%s(%s, %d+ bp, perc %d+)" % (
            self.location, self.blast_db, self.min_align_length,
            self.perc_identity
        )


class AvoidIDTHairpins(Objective):
    """Avoid Hairpin patterns as defined by the IDT guidelines.

    A hairpin is defined by a sequence segment which has a reverse complement
    "nearby" in a given window.

    Parameters
    ----------

    stem_size
      Size of the stem of a hairpin, i.e. the length of the sequence which
      should have a reverse complement nearby to be considered a hairpin.

    hairpin_window
      The window in which the stem's reverse complement should be searched for.

    boost
      Multiplicative factor, importance of this objective in a multi-objective
      optimization.
    """

    best_possible_score = 0

    def __init__(self, stem_size=20, hairpin_window=200, boost=1.0):

        self.stem_size = stem_size
        self.hairpin_window = hairpin_window
        self.boost = boost

    def evaluate(self, problem):
        sequence = problem.sequence
        reverse = reverse_complement(sequence)
        locations = []
        for i in range(len(sequence) - self.hairpin_window):
            word = sequence[i:i + self.stem_size]
            rest = reverse[-(i + self.hairpin_window):-(i + self.stem_size)]
            if word in rest:
                locations.append([i, i + self.hairpin_window])
        score = -len(locations)

        locations = sorted([Location(l[0], l[1]) for l in locations])

        return ObjectiveEvaluation(self, problem, score, locations=locations)

    def localized(self, location):
        # TODO: I'm pretty sure this can be localized
        return self

    def __repr__(self):
        return "NoHairpinsIDTObjective(size=%d, window=%d)" % \
            (self.stem_size, self.hairpin_window)


class AvoidNonuniqueSegments(Objective):
    """Avoid sub-sequence which have repeats elsewhere in the sequence.

    Parameters
    ----------

    min_length
      Minimal length of sequences to be considered repeats

    location
      Segment of the sequence in which to look for repeats. If None, repeats
      are searched in the full sequence.

    include_reverse_complement
      If True, the sequence repeats are also searched for in the reverse
      complement of the sequence (or sub sequence if `location` is not None).

    Examples
    --------

    >>> from dnachisel import *
    >>> sequence = random_dna_sequence(50000)
    >>> constraint= AvoidNonuniqueSegments(10, include_reverse_complement=True)
    >>> problem = DnaOptimizationProblem(sequence, constraints= [contraint])
    >>> print (problem.constraints_summary())
    """

    def __init__(self, min_length, location=None,
                 include_reverse_complement=False):
        self.min_length = min_length
        self.location = location
        self.include_reverse_complement = include_reverse_complement

    def evaluate(self, problem):
        """Return 1 if the sequence has no repeats, else -N where N is the
        number of non-unique segments in the sequence."""
        if self.location is not None:
            location = self.location
        else:
            location = Location(0, len(problem.sequence))
        sequence = location.extract_sequence(problem.sequence)
        rev_complement = reverse_complement(sequence)
        kmers_locations = defaultdict(lambda: [])
        for i in range(len(sequence) - self.min_length):
            start, end = i, i + self.min_length
            kmers_locations[sequence[start:end]].append((start, end))
        if self.include_reverse_complement:
            for i in range(len(sequence) - self.min_length):
                start, end = i, i + self.min_length
                kmers_locations[rev_complement[start:end]].append(
                    (len(sequence) - end, len(sequence) - start)
                )

        locations = sorted([
            Location(*min(positions_list, key=lambda p: p[0]))
            for positions_list in kmers_locations.values()
            if len(positions_list) > 1
        ])

        if locations == []:
            return ObjectiveEvaluation(
                self, problem, score=1,
                message="Passed: no nonunique %d-mer found." % self.min_length)

        return ObjectiveEvaluation(
            self, problem, score=-len(locations),
            locations=locations,
            message="Failed, the following positions are the first occurences"
                    "of non-unique segments %s" % locations)

    def __repr__(self):
        return "NoNonuniqueKmers(%d)" % (self.min_length)


class AvoidPattern(PatternObjective):
    """Enforce that the given pattern is absent in the sequence."""
    best_possible_score = 0

    def evaluate(self, problem):
        locations = self.pattern.find_matches(problem.sequence, self.location)
        score = -len(locations)
        if score == 0:
            message = "Passed. Pattern not found !"
        else:
            message = "Failed. Pattern found at positions %s" % locations
        return ObjectiveEvaluation(
            self, problem, score, locations=locations, message=message
        )

    def __repr__(self):
        return "AvoidPattern(%s, %s)" % (self.pattern, self.location)



class CodonOptimize(Objective):
    """Objective to codon-optimize a coding sequence for a particular organism.

    Several codon-optimization policies exist. At the moment this Objective
    implements a method in which codons are replaced by the most frequent
    codon in the species.

    (as long as this doesn't break any Objective or lowers the global
    optimization objective)

    Supported organisms are ``E. coli``, ``S. cerevisiae``, ``H. Sapiens``,
    ``C. elegans``, ``D. melanogaster``, ``B. subtilis``.

    Parameters
    ----------

    organism
      Name of the organism to codon-optimize for. Supported organisms are
      ``E. coli``, ``S. cerevisiae``, ``H. Sapiens``, ``C. elegans``,
      ``D. melanogaster``, ``B. subtilis``.
      Note that the organism can be omited if a ``codon_usage_table`` is
      provided instead

    location
      Pair (start, end) indicating the position of the gene to codon-optimize.
      If not provided, the whole sequence is considered as the gene. Make
      sure the length of the sequence in the location is a multiple of 3.
      The location strand is either 1 if the gene is encoded on the (+) strand,
      or -1 for antisense.

    codon_usage_table
      A dict of the form ``{"TAC": 0.112, "CCT": 0.68}`` giving the RSCU table
      (relative usage of each codon). Only provide if no ``organism`` name
      was provided.

    Examples
    --------

    >>> objective = CodonOptimizationObjective(
    >>>     organism = "E. coli",
    >>>     location = (150, 300), # coordinates of a gene
    >>>     strand = -1
    >>> )


    """

    def __init__(self, organism=None, location=None,
                 codon_usage_table=None, boost=1.0):
        self.boost = boost
        self.location = location
        self.organism = organism
        if organism is not None:
            codon_usage_table = CODON_USAGE_TABLES[self.organism]
        if codon_usage_table is None:
            raise ValueError("Provide either an organism name or a codon "
                             "usage table")
        self.codon_usage_table = codon_usage_table

    def evaluate(self, problem):

        location = (self.location if self.location is not None
                    else Location(0, len(problem.sequence)))
        subsequence = (location.extract_sequence(problem.sequence)
                       .replace("T", "U"))
        length = len(subsequence)
        if (length % 3):
            raise ValueError("CodonOptimizationObjective on a window/sequence"
                             "with size %d not multiple of 3)" % length)
        score = sum([
            self.codon_usage_table[subsequence[3 * i:3 * (i + 1)]]
            for i in range(int(length / 3))
        ])
        return ObjectiveEvaluation(
            self, problem, score, locations=[location],
            message="Codon opt. on window %s scored %.02E" %
                    (location, score)
        )

    def __str__(self):
        return "CodonOptimize(%s, %s)" % (str(self.location), self.organism)

    def __repr__(self):
        return str(self)


class DoNotModify(Objective):
    """Specify that some locations of the sequence should not be changed.

    ``DoNotModify`` Objectives are used to constrain the mutations space
    of DNA OptimizationProblem.

    Parameters
    ----------

    location
      Location object indicating the position of the segment that
      must be left unchanged.
    """

    best_possible_score = 1

    def __init__(self, location=None, indices=None, boost=1.0):
        self.location = location
        self.indices = np.array(indices)
        self.boost = boost

    def evaluate(self, problem):
        sequence = problem.sequence
        original = problem.sequence_before
        if (self.location is None) and (self.indices is None):
            return ObjectiveEvaluation(sequence == original,
                                       locations=[self.location])
        elif self.location is not None:
            subseq, suboriginal = [self.location.extract_sequence(s)
                                   for s in (sequence, original)]
            score = 1 if (subseq == suboriginal) else -1
            return ObjectiveEvaluation(self, problem, score,
                                       locations=[self.location])
        else:
            sequence = np.fromstring(sequence, dtype="uint8")
            original = np.fromstring(original, dtype="uint8")
            if (sequence[self.indices] == original[self.indices]).min():
                score = 1
            else:
                score = -1

            return ObjectiveEvaluation(self, problem, score,
                                       locations=[self.location])

    def localized(self, location):
        """Localize the DoNotModify to the overlap of its location and the new.
        """
        if self.location is not None:
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return VoidObjective(parent_objective=self)
            return self.copy_with_changes(location=new_location)
        else:
            start, end = location.start, location.end
            inds = self.indices
            new_indices = inds[(start <= inds) & (inds <= end)]
            return self.copy_with_changes(indices=new_indices)

    def __repr__(self):
        return "DoNotModify(%s)" % str(self.location)


class EnforceGCContent(Objective):
    """Objective on the local or global proportion of G/C nucleotides.

    Examples
    --------

    >>> # Enforce global GC content between 40 and 70 percent.
    >>> Objective = GCContentObjective(0.4, 0.7)
    >>> # Enforce 30-80 percent local GC content over 50-nucleotides windows
    >>> Objective = GCContentObjective(0.3, 0.8, gc_window=50)


    Parameters
    ----------

    gc_min
      Minimal proportion of G-C (e.g. ``0.35``)

    gc_max
      Maximal proportion of G-C (e.g. ``0.75``)

    gc_window
      Length of the sliding window, in nucleotides, for local GC content.
      If not provided, the global GC content of the whole sequence is
      considered

    location
      Location objet indicating that the Objective only applies to a
      subsegment of the sequence. Make sure it is bigger than ``gc_window``
      if both parameters are provided

    """

    def __init__(self, gc_min=0, gc_max=1.0, gc_objective=None,
                 gc_window=None, location=None, boost=1.0):
        if gc_objective is not None:
            gc_min = gc_max = gc_objective
        self.gc_objective = gc_objective
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.gc_window = gc_window
        self.location = location
        self.boost = boost

    def evaluate(self, problem):
        location = (self.location if self.location is not None
                    else Location(0, len(problem.sequence)))
        wstart, wend = location.start, location.end
        sequence = location.extract_sequence(problem.sequence)
        gc = gc_content(sequence, self.gc_window)
        breaches = (np.maximum(0, self.gc_min - gc) +
                    np.maximum(0, gc - self.gc_max))
        score = - (breaches.sum())
        breaches_starts = (breaches > 0).nonzero()[0]

        if len(breaches_starts) == 0:
            breaches_locations = []
        elif len(breaches_starts) == 1:
            if self.gc_window is not None:
                start = breaches_starts[0]
                breaches_locations = [[start, start + self.gc_window]]
            else:
                breaches_locations = [[wstart, wend]]
        else:
            breaches_locations = []
            current_start = breaches_starts[0]
            last_end = current_start + self.gc_window
            for i in breaches_starts[1:]:
                if (i > last_end + self.gc_window):
                    breaches_locations.append([
                        wstart + current_start, wstart + last_end]
                    )
                    current_start = i
                    last_end = i + self.gc_window

                else:
                    last_end = i + self.gc_window
            breaches_locations.append(
                [wstart + current_start, wstart + last_end])

        if breaches_locations == []:
            message = "Passed !"
        else:
            breaches_locations = [Location(*loc) for loc in breaches_locations]
            message = ("Failed: GC content out of bound on segments " +
                       ", ".join([str(l) for l in breaches_locations]))
        return ObjectiveEvaluation(self, problem, score, breaches_locations,
                                   message=message)

    def localized(self, location):
        """Localize the GC content evaluation

        For a location, the GC content evaluation will be restricted
        to [start - gc_window, end + gc_window]
        """
        if self.location is not None:
            if self.gc_window is None:
                return self
            new_location = self.location.overlap_region(location)
            if new_location is None:
                return VoidObjective(parent_objective=self)
            else:
                extension = 0 if self.gc_window is None else self.gc_window - 1
                extended_location = location.extended(extension)

                new_location = self.location.overlap_region(extended_location)
        else:
            if self.gc_window is not None:
                new_location = location.extended(self.gc_window + 1)
            else:
                new_location = None
        return self.copy_with_changes(location=new_location)

    def __repr__(self):
        return (
            "EnforceGCContent(min %.02f, max %.02f, gc_win %s, location %s)" %
            (self.gc_min, self.gc_max,
             "global" if (self.gc_window is None) else self.gc_window,
             self.location))


class EnforcePattern(PatternObjective):
    """Enforce that the given pattern is present in the sequence.

    Parameters
    ----------

    pattern
      A SequencePattern or DnaNotationPattern

    dna_pattern
      A string of ATGC that will be converted automatically to a DNA pattern

    enzyme
      Enzyme name, can be provided instead of pattern or dna_pattern

    location
      Location of the DNA segment on which to enforce the pattern e.g.
      ``Location(10, 45, 1)``
    """
    best_possible_score = 0
    shrink_when_localized = False

    def __init__(self, pattern=None, dna_pattern=None, enzyme=None,
                 location=None, occurences=1, boost=1.0):
        PatternObjective.__init__(self, pattern=pattern, location=location,
                                  dna_pattern=dna_pattern,
                                  enzyme=enzyme)
        self.occurences = occurences
        self.boost = boost

    def evaluate(self, problem):
        location = (self.location if (self.location is not None) else
                    Location(0, len(problem.sequence)))
        locations = self.pattern.find_matches(problem.sequence, location)
        score = -abs(len(locations) - self.occurences)

        if score == 0:
            message = "Passed. Pattern found at positions %s" % locations
        else:
            if self.occurences == 0:
                message = "Failed. Pattern not found."
            else:
                message = ("Failed. Pattern found %d times instead of %d"
                           " wanted at positions %s") % (len(locations),
                                                         self.occurences,
                                                         location)
        return ObjectiveEvaluation(
            self, problem, score, message=message,
            locations=None if location is None else [location],
        )

    def __repr__(self):
        return "EnforcePattern(%s, %s)" % (self.pattern, self.location)


class EnforceRegionsCompatibility(Objective):
    max_possible_score = 0

    def __init__(self, regions, compatibility_condition, boost=1.0):
        self.regions = regions
        self.compatibility_condition = compatibility_condition
        self.boost = boost

    def evaluate(self, problem):
        incompatible_regions_pairs = []
        for (r1, r2) in itertools.combinations(self.regions, 2):
            if not self.compatibility_condition(r1, r2, problem):
                incompatible_regions_pairs.append((r1, r2))

        all_regions_with_incompatibility = [
            region
            for incompatibles_pair in incompatible_regions_pairs
            for region in incompatibles_pair
        ]
        counter = Counter(all_regions_with_incompatibility)
        all_regions_with_incompatibility = sorted(
            list(set(all_regions_with_incompatibility)),
            key=counter.get
        )
        all_regions_with_incompatibility = [
            Location(*r) for r in all_regions_with_incompatibility
        ]

        score = -len(incompatible_regions_pairs)
        if score == 0:
            message = "All compatible !"
        else:
            message = "Found the following imcompatibilities: %s" % (
                incompatible_regions_pairs
            )
        return ObjectiveEvaluation(
            self, problem,
            score=score,
            locations=all_regions_with_incompatibility,
            message=message
        )

    def localized(self, window):
        # FIXME: weird stuff here
        wstart, wend = window
        included_regions = [
            (a, b) for (a, b) in self.regions
            if wstart <= a <= b <= wend
        ]

        def evaluate(problem):
            """Objective evaluation"""
            # compute incompatibilities but exclude to
            # consider pairs of regions that are both
            # outside the current localization window
            # FIXME: weird stuff here
            incompatible_regions = [
                region
                for region in self.regions
                for included_region in included_regions
                if (region != included_region) and
                not self.compatibility_condition(
                    included_region, region, problem
                )
            ]
            score = -len(incompatible_regions)
            return ObjectiveEvaluation(
                self, problem,
                score=score
            )
        return Objective(evaluate, boost=self.boost)

    def __repr__(self):
        return "CompatSeq(%s...)" % str(self.regions[0])


class EnforceTerminalGCContent(TerminalObjective):
    """Enforce bounds for the GC content at the sequence's terminal ends.

    Parameters
    ----------

    window_size
      Size in basepair of the two terminal ends to consider

    gc_min
      A float between 0 and 1, minimal proportion of GC that the ends should
      contain

    gc_max
      Float between 0 and 1, maximal proportion of GC that the ends should
      contain

    boost
      Multiplicatory factor applied to this objective.
    """

    def __init__(self, window_size, gc_min=0, gc_max=1, boost=1.0):
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.window_size = window_size
        self.boost = boost

    def evaluate_end(self, sequence):
        return (self.gc_min < gc_content(sequence) < self.gc_max)

    def __repr__(self):
        return "Terminal(%.02f < gc < %.02f, window: %d)" % \
            (self.gc_min, self.gc_max, self.window_size)

class EnforceTranslation(Objective):
    """Enforce that the DNA segment sequence translates to a specific
    amino-acid sequence.

    This class enforces the standard translation, but it is also possible to
    change the class' `codons_sequences` and `codons_translations`
    dictionnaries for more exotic kind of translations

    Parameters
    -----------

    location
      A pair (start, end) indicating the segment that is a coding sequence

    strand
      Set to 1 (default) if the gene is read in direct sense, -1 for antisense

    translation
      String representing the protein sequence that the DNA segment should
      translate to, eg. "MKY...LL*" ("*" stands for stop codon).
      This parameter can be omited if the ``sequence`` parameter is provided

    Examples
    --------

    >>> from dnachisel import *
    >>> sequence = some_dna_sequence # with a gene in segment 150-300
    >>> Objective = EnforceTranslationObjective(
    >>>     window=(150,300),
    >>>     strand = 1,
    >>>     translation= translate(sequence[150:300]) # "MKKLYQ...YNL*"
    >>> )
    >>> # OR EQUIVALENT IF THE GENE ALREADY ENCODES THE RIGHT PROTEIN:
    >>> Objective = EnforceTranslationObjective(
    >>>     window=(150,300),
    >>>     strand = 1,
    >>>     sequence = sequence
    >>> )
    """

    best_possible_score = 0
    codons_sequences = CODONS_SEQUENCES
    codons_translations = "Bacterial"

    def __init__(self, location=None, translation=None, boost=1.0):
        self.translation = translation
        self.set_location(location)
        self.boost = boost

        self.initialize_translation_from_problem = (translation is None)
        self.initialize_location_from_problem = (location is None)

    def set_location(self, location):
        if location is not None:
            if len(location) % 3:
                raise ValueError(
                    "Location length in Codon Objectives should be a 3x. "
                    "Location %s has length %d" % (location, len(location))
                )

            if ((self.translation is not None) and
                    (len(location) != 3 * len(self.translation))):
                raise ValueError(
                    ("Window size (%d bp) incompatible with translation "
                     "(%d aa)") % (len(location), len(self.translation))
                )
        self.location = location

    def set_translation(self, translation):
        if (translation is not None) and (self.location is not None):
            if (len(self.location) != 3 * len(self.translation)):
                raise ValueError(
                    ("Window size (%d bp) incompatible with translation "
                     "(%d aa)") % (len(self.location), len(self.translation))
                )
        self.translation = translation

    def initialize_problem(self, problem, role):
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            result = self.copy_with_changes()
            result.set_location(location)
        else:
            result = self
        if self.translation is None:
            subsequence = result.location.extract_sequence(problem.sequence)
            translation = translate(subsequence, self.codons_translations)

            result = self.copy_with_changes(translation=translation)
        return result

    def evaluate(self, problem):
        location = (self.location if self.location is not None else
                    Location(0, len(problem.sequence)))
        subsequence = location.extract_sequence(problem.sequence)
        translation = translate(subsequence, self.codons_translations)
        errors = [
            ind
            for ind in range(len(translation))
            if translation[ind] != self.translation[ind]
        ]
        errors_locations = [
            Location(3 * ind, 3 * (ind + 1)) if self.location.strand >= 0 else
            Location(self.location.end - 3 * (ind + 1),
                     self.location.end - 3 * ind)
            for ind in errors
        ]
        success = (len(errors) == 0)
        return ObjectiveEvaluation(self, problem, score=-len(errors),
                                   locations=errors_locations,
                                   # self.location],
                                   message="All OK." if success else
                                   "Wrong translation at indices %s"
                                   % errors)


    def localized(self, location):
        """"""
        if self.location is not None:
            overlap = self.location.overlap_region(location)
            if overlap is None:
                return VoidObjective(parent_objective=self)
            else:
                # return self
                o_start, o_end = overlap.start, overlap.end
                w_start, w_end = self.location.start, self.location.end

                if self.location.strand == 1:
                    start_codon = int((o_start - w_start) / 3)
                    end_codon = int((o_end - w_start) / 3)
                    new_location = Location(
                        start=w_start + 3 * start_codon,
                        end=min(w_end, w_start + 3 * (end_codon + 1)),
                        strand=self.location.strand
                    )
                else:
                    start_codon = int((w_end - o_end) / 3)
                    end_codon = int((w_end - o_start) / 3)
                    new_location = Location(
                        start=max(w_start, w_end - 3 * (end_codon + 1)),
                        end=w_end - 3 * start_codon,
                        strand=self.location.strand
                    )

                new_translation = self.translation[start_codon:end_codon + 1]
                return self.__class__(new_location,
                                      translation=new_translation,
                                      boost=self.boost)
        return self

    def __repr__(self):
        return "EnforceTranslation(%s)" % str(self.location)


class AvoidStopCodon(EnforceTranslation):
    """Objective: do not introduce any new stop codon in that frame."""
    codons_translations = {
        codon: "*" if (translation == '*') else "_"
        for codon, translation in CODONS_TRANSLATIONS.items()
    }
    codons_sequences = None

    def __str__(self):
        return "AvoidStopCodon(%s)" % self.location

    def __init__(self, location, translation=None, boost=1.0):
        if (len(location) % 3) != 0:
            raise ValueError("Loc. %s for AvoidStopCodon is not 3x" % location)
        self.translation = '_' * int(len(location) / 3)
        self.boost = boost
        self.location = location


class MinimizeDifferences(Objective):
    """Objective to minimize the differences to a given sequence.


    This can be used to enforce "conservative" optimization, in which we try
    to minimize the changes from the original sequence

    Parameters
    ----------

    location
      Pair (start, end) indicating the segment of the sequence. If none
      provided, the whole sequence is considered.

    target_sequence
      The DNA sequence that the problem' sequence (or subsequence) should equal.
      Can be omitted if ``sequence_before`` is provided instead

    sequence_before
      A DNA sequence (will generally be the problem' sequence itself) with
      already the right sequence at the given ``location``. Only provide if
      you are not providing a ``target_sequence``


    Examples
    --------

    >>> from dnachisel import *
    >>> sequence = random_dna_sequence(length=10000)
    >>> # Fix the sequence's local gc content while minimizing changes.
    >>> problem = DnaOptimizationProblem(
    >>>     sequence = sequence,
    >>>     Objectives = [GCContentObjective(0.3,0.6, gc_location=50)],
    >>>     objective = [MinimizeDifferencesObjective(
    >>>                     sequence_before=sequence)]
    >>> )
    >>> problem.solve_all_Objectives_one_by_one()
    >>> problem.maximize_all_objectives_one_by_one()

    """

    best_possible_score = 0

    def __init__(self, location=None, target_sequence=None, boost=1.0):
        self.boost = boost
        self.location = location
        self.target_sequence = target_sequence
        self.initialize_sequence_from_problem = (target_sequence is None)
        self.reference_sequence = target_sequence

    def initialize_problem(self, problem, role):
        if not self.initialize_sequence_from_problem:
            return self
        if self.location is None:
            self.location = Location(0, len(problem.sequence))
        reference_sequence = self.location.extract_sequence(problem.sequence)
        if self.location.strand == -1:
            reference_sequence = reverse_complement(reference_sequence)
        return self.copy_with_changes(
            reference_sequence=reference_sequence,
            location=Location(self.location.start, self.location.end, 1)
        )

    def evaluate(self, problem):
        subsequence = self.location.extract_sequence(problem.sequence)
        locations = [
            Location(self.location.start + i,
                     self.location.start + i, 1)
            for i in range(len(subsequence))
            if subsequence[i] != self.reference_sequence[i]
        ]

        diffs = len(locations)
        #sequences_differences(subsequence, self.reference_sequence)
        return ObjectiveEvaluation(
            self, problem, score=-diffs, locations=locations,
            message="Found %d differences with target sequence" % diffs
        )

    def localized(self, location):
        new_location = location.overlap_region(self.location)
        if new_location is None:
            return VoidObjective(parent_objective=self)
        target_start = new_location.start - self.location.start
        target_end = new_location.end - self.location.end
        new_target = self.reference_sequence[target_start:target_end]
        return MinimizeDifferences(location=new_location, boost=self.boost,
                                   target_sequence=new_target)

    def __str__(self):
        return "MinimizeDifferencesObj(%s, %s...)" % (
            "global" if self.location is None else str(self.location),
            self.sequence[:7]
        )


class SequenceLengthBounds(Objective):
    """Checks that the sequence length is between bounds.

    Quite an uncommon objective as it can't really be solved or optimized.
    But practical at times, as part of a list of constraints to verify.

    Parameters
    ----------

    min_length
      Minimal allowed sequence length in nucleotides

    max_length
      Maximal allowed sequence length in nucleotides. None means no bound.
    """
    best_possible_score = 0

    def __init__(self, min_length=0, max_length=None):
        self.min_length = min_length
        self.max_length = max_length

    def evaluate(self, problem):
        """Return 0 if the sequence length is between the bounds, else -1"""
        L, mini, maxi = len(problem.sequence), self.min_length, self.max_length
        if maxi is None:
            score = (L >= mini)
        else:
            score = (mini <= L <= maxi)
        return ObjectiveEvaluation(self, problem, score - 1)

    def __repr__(self):
        return "Length(%d < L < %d)" % (self.min_length, self.max_length)
