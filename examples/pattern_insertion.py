"""Pattern insertion in DnaChisel.

"""
from dnachisel import *

protein_sequence = random_protein_sequence(length=300, seed=123)
dna_sequence = reverse_translate(protein_sequence)

canvas = DnaCanvas(
    sequence=dna_sequence,
    constraints = [
        NoPatternConstraint(homopolymer_pattern("A",5)),
        NoPatternConstraint(enzyme_pattern("HindIII")),
        EnforceTranslationConstraint(window=[0, len(dna_sequence)],
                                     translation=protein_sequence),
    ]
)

pattern =enzyme_pattern("BsaI")
matches = pattern.find_matches(dna_sequence)
print ("\nBEFORE PATTERN INSERTION AND SOLVING\n")
print ("BsaI matches in the original sequence: %s" % matches)
canvas.print_constraints_summary()
canvas.include_pattern_by_successive_tries(enzyme_pattern("BsaI"))
print ("AFTER PATTERN INSERTION AND SOLVING")
canvas.print_constraints_summary()
