"Implement EnforceTranslation."


from ..Specification import SpecEvaluation
from ..biotools import (
    translate,
    reverse_complement,
    get_backtranslation_table,
)
from ..Location import Location
from .CodonSpecification import CodonSpecification


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

    genetic_table
      Either "Standard", "Bacterial", or any other Biopython genetic table name
      (see dnachisel.biotools.CODON_TABLE_NAMES for a list of accepted names).

    start_codon
      Signals that the first codon is a start codon and provides a policy for
      changing it. It is very important that this parameter be set when dealing
      with full coding sequences in organisms with different start codons, e.g.
      GTG in bacteria. If it is not set, then a GTG (start codon or Valine)
      could be changed into GTA (also valine but NOT a start codon).
      Can be False (the first codon is not considered a start codon, but can
      still naturally code for methionine), ``keep`` (freezes the original
      sub-sequence at this location, or a codon e.g. ``ATG`` or ``GTG`` or
      ``["ATG", "GTG"]``.

    translation
      String representing the protein sequence that the DNA segment should
      translate to, eg. "MKY...LL*" ("*" stands for stop codon).
      Can be omitted if the sequence initially encodes the right sequence.

    boost
      Score multiplicator (=weight) for when the specification is used as an
      optimization objective alongside competing objectives.
    """

    best_possible_score = 0
    enforced_by_nucleotide_restrictions = True
    shorthand_name = "cds"
    default_genetic_table = "Standard"

    def __init__(
        self,
        genetic_table="default",
        start_codon=None,
        translation=None,
        location=None,
        boost=1.0,
    ):
        """Initialize."""
        self.translation = translation
        location = Location.from_data(location)
        if (location is not None) and (location.strand not in [-1, 1]):
            location = Location(location.start, location.end, 1)
        self.set_location(location)
        self.start_codon = start_codon
        self.boost = boost
        if genetic_table == "default":
            genetic_table = self.default_genetic_table
        self.genetic_table = genetic_table

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

    def initialized_on_problem(self, problem, role):
        """Get translation from the sequence if it is not already set."""
        result = self._copy_with_full_span_if_no_location(problem)
        if result.translation is None:
            subsequence = result.location.extract_sequence(problem.sequence)
            translation = translate(
                subsequence,
                table=result.genetic_table,
                assume_start_codon=result.start_codon is not None,
            )
            result = result.copy_with_changes(translation=translation)
        if len(result.location) != 3 * len(result.translation):
            raise ValueError(
                (
                    "Window size (%d bp) incompatible with translation "
                    "(%d aa)"
                )
                % (len(result.location), len(result.translation))
            )
        if (result.start_codon is not None) and result.translation[0] != "M":
            raise ValueError(
                (
                    "Spec. %s specificies a start_codon parameter, but the "
                    "translation at this location does not start with Met. "
                    "Maybe you should provide a 'translation' parameter or "
                    "set up a different Genetic Code."
                )
                % result.label(use_short_form=True)
            )

        return result

    def evaluate(self, problem):
        """Score is the number of wrong-translation codons."""
        # Note: this method is actually very little used as this specification
        # class sets the enforced_by_nucleotide_restrictions attribute.
        location = (
            self.location
            if self.location is not None
            else Location(0, len(problem.sequence), 1)
        )
        subsequence = location.extract_sequence(problem.sequence)
        translation = translate(
            subsequence,
            table=self.genetic_table,
            assume_start_codon=self.start_codon is not None,
        )
        errors_locations = [
            self.codon_index_to_location(index)
            for (index, amino_acid) in enumerate(translation)
            if amino_acid != self.translation[index]
        ]
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
        if self.location.strand == -1:
            location_is_at_start = new_location.end >= self.location.end
        else:
            location_is_at_start = new_location.start <= self.location.start
        return self.__class__(
            location=new_location,
            translation=new_translation,
            boost=self.boost,
            genetic_table=self.genetic_table,
            start_codon=self.start_codon if location_is_at_start else None
            # has_start_codon=self.has_start_codon and location_is_at_start,
        )

    def restrict_nucleotides(self, sequence, location=None):
        if self.backtranslation_table is None:
            return []

        def get_first_codon_choices(first_codon):
            if self.start_codon is None:
                return self.backtranslation_table[self.translation[0]]
            elif isinstance(self.start_codon, (list, tuple)):
                return list(self.start_codon)
            elif self.start_codon == "keep":
                return [first_codon]
            else:
                return [self.start_codon]  # "ATG"
        first_codon_location = self.codon_index_to_location(0)
        first_codon = first_codon_location.extract_sequence(sequence)
        choices = [
            (first_codon_location, get_first_codon_choices(first_codon))
        ] + [
            (self.codon_index_to_location(i), self.backtranslation_table[aa])
            for i, aa in list(enumerate(self.translation))[1:]
        ]

        def standardize_choice(choice):
            location, choices_list = choice
            if location.strand == -1:
                choices_list = [reverse_complement(c) for c in choices_list]
            return (location.to_tuple()[:2], choices_list)

        return sorted([standardize_choice(choice) for choice in choices])

    def __repr__(self):
        return "EnforceTranslation(%s)" % str(self.location)

    def __str__(self):
        return "EnforceTranslation(%s)" % str(self.location)

    def short_label(self):
        return "cds"
    
    def breach_label(self):
        return "protein sequence changed"
