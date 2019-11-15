
from ..biotools import translate
from ..Location import Location
from ..Specification import SpecEvaluation
from .CodonSpecification import CodonSpecification

from Bio.Data import CodonTable


def get_backtranslation_table(table_name="Standard"):
    table = CodonTable.unambiguous_dna_by_name[table_name]
    back_translation_table = {}
    for codon, amino_acid in table.forward_table.items():
        if amino_acid not in back_translation_table:
            back_translation_table[amino_acid] = []
        back_translation_table[amino_acid].append(codon)
    back_translation_table["*"] = table.stop_codons
    back_translation_table["START"] = table.start_codons
    return back_translation_table


def get_genetic_table(table_name="Standard"):
    table = CodonTable.unambiguous_dna_by_name[table_name]
    genetic_table = {stop_codon: "*" for stop_codon in table.stop_codons}
    genetic_table.update(table.forward_table)
    return genetic_table


class AvoidStopCodons(CodonSpecification):
    """Do not introduce any new stop codon in that frame.

    This can be used for research purposes, to avoid breaking a reading frame
    when editing it with quasi-synonymous mutations.
    """

    def __init__(self, genetic_table="Standard", location=None, boost=1.0):
        self.genetic_table = genetic_table
        self.boost = boost
        self.location = Location.from_data(location)

    def initialized_on_problem(self, problem, role):
        """Get translation from the sequence if it is not already set."""
        return self._copy_with_full_span_if_no_location(problem)

    def localized_on_window(self, new_location, start_codon, end_codon):
        return self.copy_with_changes(location=new_location)

    def evaluate(self, problem):
        location = (
            self.location
            if self.location is not None
            else Location(0, len(problem.sequence))
        )
        subsequence = location.extract_sequence(problem.sequence)
        translation = translate(subsequence, table=self.genetic_table)
        errors_locations = [
            self.codon_index_to_location(index)
            for index in range(len(translation))
            if translation[index] == "*"
        ]
        return SpecEvaluation(
            self,
            problem,
            score=-len(errors_locations),
            locations=errors_locations,
            message="All OK."
            if len(errors_locations) == 0
            else "Stop codons found at indices %s" % errors_locations,
        )

        

    def __str__(self):
        """Represent."""
        return "AvoidStopCodons(%s)" % self.location

    def __str__(self):
        """Represent."""
        return "AvoidStopCodons(%s)" % self.location
