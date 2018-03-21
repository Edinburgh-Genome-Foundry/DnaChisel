"""Biologically-related useful methods."""

from .biotools import (
    blast_sequence,
    change_biopython_record_sequence,
    complement,
    crop_record,
    dna_pattern_to_regexpr,
    find_specification_in_feature,
    gc_content,
    group_nearby_indices,
    group_nearby_segments,
    is_palyndromic,
    list_common_enzymes,
    load_record,
    random_dna_sequence,
    random_protein_sequence,
    reverse_complement,
    reverse_translate,
    sequences_differences,
    sequences_differences_array,
    sequences_differences_segments,
    sequence_to_biopython_record,
    subdivide_window,
    translate,
    windows_overlap,
    codons_frequencies_and_positions
)

from .features_annotations import (
    annotate_record,
    annotate_differences,
    annotate_pattern_occurrences
)

from .biotables import (CODONS_SEQUENCES, CODONS_TRANSLATIONS,
                        CODON_USAGE_TABLES, IUPAC_NOTATION,
                        NUCLEOTIDE_TO_REGEXPR, CODON_USAGE_BY_AA)
