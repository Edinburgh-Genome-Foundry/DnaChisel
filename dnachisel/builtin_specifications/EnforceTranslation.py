"Implement EnforceTranslation."

from .CodonSpecification import CodonSpecification

# from .VoidSpecification import VoidSpecification
from ..SpecEvaluation import SpecEvaluation
from dnachisel.biotools import (
    CODONS_SEQUENCES,
    translate,
    reverse_complement,
    get_backtranslation_table,
)
from dnachisel.Location import Location

from Bio.Data import CodonTable


class EnforceTranslation(CodonSpecification):
    """Enforce a specific amino-acid sequence translation.

    This class enforces the standard translation, but it is also possible to
    change the class' `codons_sequences` and `codons_translations`
    dictionnaries for more exotic kind of translations

    Shorthand for annotations: "cds".

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
    
    table

    """

    best_possible_score = 0
    # codons_sequences = CODONS_SEQUENCES
    enforced_by_nucleotide_restrictions = True
    shorthand_name = "cds"
    default_genetic_table = "Standard"

    def __init__(
        self,
        genetic_table="default",
        has_start_codon=False,
        translation=None,
        location=None,
        boost=1.0,
    ):
        """Initialize."""
        self.translation = translation
        if isinstance(location, tuple):
            location = Location.from_tuple(location, default_strand=+1)
        if (location is not None) and (location.strand not in [-1, 1]):
            location = Location(location.start, location.end, 1)
        self.set_location(location)
        self.boost = boost
        if genetic_table == "default":
            genetic_table = self.default_genetic_table
        self.genetic_table = genetic_table
        self.has_start_codon = has_start_codon

        self.initialize_translation_from_problem = translation is None
        self.initialize_location_from_problem = location is None
        self.backtranslation_table = get_backtranslation_table(genetic_table)

    def set_location(self, location):
        """Check that the location length is valid before setting it."""
        if location is not None:
            if len(location) % 3:
                raise ValueError(
                    "Location length in Codon Specifications should be a 3x. "
                    "Location %s has length %d" % (location, len(location))
                )

            if (self.translation is not None) and (
                len(location) != 3 * len(self.translation)
            ):
                raise ValueError(
                    (
                        "Window size (%d bp) incompatible with translation "
                        "(%d aa)"
                    )
                    % (len(location), len(self.translation))
                )
        self.location = location

    def set_translation(self, translation):
        """Check that the translation length is valid before setting it."""
        if (translation is not None) and (self.location is not None):
            if len(self.location) != 3 * len(self.translation):
                raise ValueError(
                    (
                        "Window size (%d bp) incompatible with translation "
                        "(%d aa)"
                    )
                    % (len(self.location), len(self.translation))
                )
        self.translation = translation

    def initialized_on_problem(self, problem, role):
        """Get translation from the sequence if it is not already set."""
        if self.location is None:
            location = Location(0, len(problem.sequence), 1)
            result = self.copy_with_changes()
            result.set_location(location)
        else:
            result = self
        if result.translation is None:
            subsequence = result.location.extract_sequence(problem.sequence)
            translation = translate(
                subsequence,
                table=self.genetic_table,
                assume_start_codon=self.has_start_codon,
            )

            result = result.copy_with_changes(translation=translation)
        return result

    def evaluate(self, problem):
        """Score is the number of wrong-translation codons."""
        # print ("here")
        location = (
            self.location
            if self.location is not None
            else Location(0, len(problem.sequence))
        )
        subsequence = location.extract_sequence(problem.sequence)
        translation = translate(subsequence, table=self.genetic_table)
        errors_locations = [
            self.codon_index_to_location(index)
            for (index, amino_acid) in enumerate(translation)
            if amino_acid != self.translation[index]
        ]
        # errors = [
        #     ind
        #     for ind in range(len(translation))
        #     if translation[ind] != self.translation[ind]
        # ]
        # errors_locations = [
        #     Location(3 * ind, 3 * (ind + 1))
        #     if self.location.strand >= 0
        #     else Location(
        #         start=self.location.end - 3 * (ind + 1),
        #         end=self.location.end - 3 * ind,
        #         strand=-1,
        #     )
        #     for ind in errors
        # ]
        return SpecEvaluation(
            self,
            problem,
            score=-len(errors_locations),
            locations=errors_locations,
            message="All OK."
            if len(errors_locations) == 0
            else "Wrong translation at locations %s" % errors_locations,
        )

    def localized_on_window(self, new_location, start_codon, end_codon):
        new_translation = self.translation[start_codon:end_codon]
        location_is_at_start = (
            self.location.strand == -1
            and new_location.end >= self.location.end
        ) or (new_location.start <= self.location.start)
        return self.__class__(
            location=new_location,
            translation=new_translation,
            boost=self.boost,
            genetic_table=self.genetic_table,
            has_start_codon=self.has_start_codon and location_is_at_start,
        )

    def restrict_nucleotides(self, sequence, location=None):
        if self.backtranslation_table is None:
            return []

        def get_first_codon_choices(first_codon):
            if not self.has_start_codon:
                return self.backtranslation_table[self.translation[0]]
            start_codons = self.backtranslation_table["START"]
            if first_codon not in start_codons:
                raise ValueError(
                    (
                        "Spec. %s starts with %s, not a start codon, yet "
                        "it has parameter has_start_codon set to True"
                    )
                    % (self.label(use_short_form=True), first_codon)
                )
            standard_amino_acid = translate(
                first_codon, table=self.genetic_table
            )
            synonyms = self.backtranslation_table[standard_amino_acid]
            return [codon for codon in synonyms if codon in start_codons]

        first_codon_location = self.codon_index_to_location(0)
        first_codon = first_codon_location.extract_sequence(sequence)
        choices = [
            (first_codon_location, get_first_codon_choices(first_codon))
        ] + [
            (self.codon_index_to_location(i), self.backtranslation_table[aa])
            for i, aa in list(enumerate(self.translation))[1:]
        ]
        # print (choices)
        # print (sorted([standardize_choice(choice) for choice in choices]))
        def standardize_choice(choice):
            location, choices_list = choice
            if location.strand == -1:
                choices_list = [reverse_complement(c) for c in choices_list]
            return (location.to_tuple()[:2], choices_list)

        return sorted([standardize_choice(choice) for choice in choices])

        # if strand == 1:
        #     if self.has_start_codon:
        #         first_codon = sequence[start : start + 3]
        #         first_choice = get_first_codon_choices(first_codon)
        #     else:
        #         first_choice = self.backtranslation_table[self.translation[0]]
        #     return [((start, start + 3), first_choice)] + [
        #         (
        #             (i, i + 3),
        #             set(
        #                 self.backtranslation_table[
        #                     self.translation[int((i - start) / 3)]
        #                 ]
        #             ),
        #         )
        #         for i in range(start + 3, end, 3)
        #     ]
        # else:
        #     if self.has_start_codon:
        #         first_codon = reverse_complement(sequence[end - 3 : end])
        #         last_choice = [
        #             codon for codon in get_first_codon_choices(first_codon)
        #         ]
        #     else:
        #         last_choice = self.backtranslation_table[self.translation[0]]
        #     last_choice = [reverse_complement(c) for c in last_choice]
        #     return [
        #         (
        #             (i, i + 3),
        #             set(
        #                 reverse_complement(n)
        #                 for n in self.backtranslation_table[
        #                     self.translation[-int((i - start) / 3) - 1]
        #                 ]
        #             ),
        #         )
        #         for i in range(start, end - 3, 3)
        #     ] + [((end - 3, end), last_choice)]

    def __repr__(self):
        return "EnforceTranslation(%s)" % str(self.location)

    def __str__(self):
        return "EnforceTranslation(%s)" % str(self.location)

    def short_label(self):
        return "cds"
