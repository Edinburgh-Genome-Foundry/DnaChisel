"""Biologically-related useful methods."""

from .blast_sequence import blast_sequence
from .genbank_operations import (
    change_biopython_record_sequence,
    find_specification_in_feature,
    load_record,
    sequence_to_biopython_record,
    annotate_record,
    annotate_differences,
    annotate_pattern_occurrences,
)

from .sequences_operations import (
    complement,
    reverse_complement,
    reverse_translate,
    dna_pattern_to_regexpr,
    translate,
)

from .enzymes_operations import list_common_enzymes

from .gc_content import gc_content

from .indices_operations import (
    group_nearby_indices,
    group_nearby_segments,
    subdivide_window,
    windows_overlap,
)
from .sequences_differences import (
    sequences_differences,
    sequences_differences_array,
    sequences_differences_segments,
)

from .formatting_operations import (
    dict_to_pretty_string,
    round_all_numbers_in_dict,
    score_to_formatted_string
)

from .random_sequences import random_dna_sequence, random_protein_sequence


from .biotables import (
    CODONS_SEQUENCES,
    CODONS_TRANSLATIONS,
    IUPAC_NOTATION,
    NUCLEOTIDE_TO_REGEXPR,
)


__all__ = [
    'blast_sequence',
    'change_biopython_record_sequence',
    'find_specification_in_feature',
    'load_record',
    'sequence_to_biopython_record',
    'annotate_record',
    'annotate_differences',
    'annotate_pattern_occurrences',
    'complement',
    'reverse_complement',
    'reverse_translate',
    'dna_pattern_to_regexpr',
    'translate',
    'list_common_enzymes',
    'gc_content',
    'group_nearby_indices',
    'group_nearby_segments',
    'subdivide_window',
    'windows_overlap',
    'sequences_differences',
    'sequences_differences_array',
    'sequences_differences_segments',
    'dict_to_pretty_string',
    'round_all_numbers_in_dict',
    'score_to_formatted_string',
    'random_dna_sequence',
    'random_protein_sequence',
    'CODONS_SEQUENCES',
    'CODONS_TRANSLATIONS',
    'IUPAC_NOTATION',
    'NUCLEOTIDE_TO_REGEXPR'
]
