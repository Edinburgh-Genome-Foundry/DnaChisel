"Implement EnforceTranslation."

from .CodonSpecification import CodonSpecification
# from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from dnachisel.biotools import CODONS_SEQUENCES, translate, reverse_complement
from dnachisel.Location import Location


class EnforceTranslation(CodonSpecification):
    """Enforce that the DNA segment sequence translates to a specific
    amino-acid sequence.

    This class enforces the standard translation, but it is also possible to
    change the class' `codons_sequences` and `codons_translations`
    dictionnaries for more exotic kind of translations

    Parameters
    -----------

    location
      Either a DnaChisel Location or a tuple of the form (start, end, strand)
      or just (start, end), with strand defaulting to +1, indicating the
      position of the gene to codon-optimize. If not provided, the whole
      sequence is considered as the gene. The location should have a length
      that is a multiple of 3. The location strand is either 1 if the gene is
      encoded on the (+) strand, or -1 for antisense.

    translation
      String representing the protein sequence that the DNA segment should
      translate to, eg. "MKY...LL*" ("*" stands for stop codon).
      Can be omitted if the sequence initially encodes the right sequence.

    """

    best_possible_score = 0
    codons_sequences = CODONS_SEQUENCES
    enforced_by_nucleotide_restrictions = True
    codons_translations = "Bacterial"

    def __init__(self, location=None, translation=None, boost=1.0):
        """Initialize."""
        self.translation = translation
        if isinstance(location, tuple):
            location = Location.from_tuple(location, default_strand=+1)
        if (location is not None) and (location.strand not in [-1, 1]):
            location = Location(location.start, location.end, 1)
        self.set_location(location)
        self.boost = boost

        self.initialize_translation_from_problem = (translation is None)
        self.initialize_location_from_problem = (location is None)

    def set_location(self, location):
        """Check that the location length is valid before setting it."""
        if location is not None:
            if len(location) % 3:
                raise ValueError(
                    "Location length in Codon Specifications should be a 3x. "
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
        """Check that the translation length is valid before setting it."""
        if (translation is not None) and (self.location is not None):
            if (len(self.location) != 3 * len(self.translation)):
                raise ValueError(
                    ("Window size (%d bp) incompatible with translation "
                     "(%d aa)") % (len(self.location), len(self.translation))
                )
        self.translation = translation

    def initialize_on_problem(self, problem, role):
        """Get translation from the sequence if it is not already set."""
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            result = self.copy_with_changes()
            result.set_location(location)
        else:
            result = self
        if result.translation is None:
            subsequence = result.location.extract_sequence(problem.sequence)
            translation = translate(subsequence, self.codons_translations)

            result = result.copy_with_changes(translation=translation)
        return result

    def evaluate(self, problem):
        """Score is the number of wrong-translation codons."""
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
            Location(start=self.location.end - 3 * (ind + 1),
                     end=self.location.end - 3 * ind,
                     strand=-1)
            for ind in errors
        ]
        success = (len(errors) == 0)
        return SpecEvaluation(self, problem, score=-len(errors),
                              locations=errors_locations,
                              message="All OK." if success else
                              "Wrong translation at indices %s" % errors)

    def localized_on_window(self, new_location, start_codon, end_codon):
        new_translation = self.translation[start_codon:end_codon]
        return self.__class__(new_location,
                              translation=new_translation,
                              boost=self.boost)

    def restrict_nucleotides(self, sequence, location=None):
        if self.codons_sequences is None:
            return []

        strand = self.location.strand
        start = self.location.start
        end = self.location.end

        if strand == 1:

            return [
                ((i, i + 3), set(self.codons_sequences[
                    self.translation[int((i - start) / 3)]
                ]))
                for i in range(start, end, 3)
            ]
        else:
            return [
                ((i, i + 3), set(
                    reverse_complement(n)
                    for n in self.codons_sequences[
                        self.translation[-int((i - start) / 3) - 1]
                    ]
                ))
                for i in range(start, end, 3)
            ]

    def __repr__(self):
        return "EnforceTranslation(%s)" % str(self.location)

    def __str__(self):
        return "EnforceTranslation(%s)" % str(self.location)
