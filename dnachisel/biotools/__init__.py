from .biotools import (
    blast_sequence,
    gc_content,
    random_dna_sequence,
    random_protein_sequence,
    reverse_complement,
    reverse_translate,
    sequences_differences,
    translate,
    windows_overlap,
    change_biopython_record_sequence,
    sequence_to_biopython_record

)

from .features_annotations import (
    annotate_record,
    annotate_differences,
    annotate_pattern_occurrences
)

from .patterns import (
    DnaNotationPattern,
    homopolymer_pattern,
    enzyme_pattern,
    repeated_kmers
)

from .biotables import CODONS_SEQUENCES, CODONS_TRANSLATIONS
