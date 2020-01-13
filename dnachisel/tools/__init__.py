from .biotables_tailor import (aa2codon_table, codon2aa_table, default_cai_table,
                        default_tai_tuller, sense_codons_list)
from .codontools import analyzeCodons, analyze_cai, get_alternate_codons,randomMutation

from .mathtools import pick_random,pick_random_tuple,hammingDistance

from .seqtools import analyze_ntcontent

from .range import Range,RangeSet

from .DBSQLite import DBSQLite


def appendLabelToDict(somedict, label):
    return {
        label + str(key):somedict[key]
        for key in somedict
    }