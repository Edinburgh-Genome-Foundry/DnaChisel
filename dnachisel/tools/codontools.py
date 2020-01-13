from .biotables_tailor import (aa2codon_table, codon2aa_table, default_cai_table,
                        default_tai_tuller, sense_codons_list)
from random import choice

from math import log,exp
def analyzeCodons(seq, positions=None, data_table=default_cai_table):
    '''
    given a sequence it returns a list with two elements: [ list_of_codons, list_of_codons_cai]
    
    seq (str):
      input dna sequence
    positions
    
    
    Usage:
    --
    >>> analyzeCodons("ATGCAGTAGCAGTGCAAGTTG")
    [['atg', 'cag', 'tag', 'cag', 'tgc', 'aag', 'ttg'], [1, 1, 1, 1, 1, 0.253, 0.02]]
    >>> analyzeCodons("ATGCAGTAGCAGTGCAAGTTG",[0,3,6])
    [['atg', 'cag', 'tag'], [1, 1, 1]]
    '''

    if positions == None:
        positions = range(0, len(seq), 3)

    seq = seq.lower()
    codons = []
    codons_cai = []
    for i in positions:
        codon = seq[i:i + 3]
        codons.append(codon)
        if codon in data_table:
            codons_cai.append(data_table[codon])
        else:
            codons_cai.append("NA")
    return [codons, codons_cai]


def get_alternate_codons(codon, data=default_tai_tuller, dist=0):
    """
    returns a alternate codon to codon
    
    data: 
        dictionary with a map between codons and tAI
    dist: 
        0   --> only synonymous codon
        1-3 --> only codon with 1-3 nt difference from original
    """
    if dist == 0:
        # return only syn codon
        return [(syn_cod, data[syn_cod])
                for syn_cod in aa2codon_table[codon2aa_table[codon]]
                if syn_cod != codon]
    else:

        def diff(str1, str2):
            nbr = 0
            for i in range(len(str1)):
                if str1[i] != str2[i]:
                    nbr += 1
            return nbr

        # return syn codon and codon 1 nt away
        return [(alt_cod, data[alt_cod]) for alt_cod in sense_codons_list
                if (alt_cod != codon and diff(codon, alt_cod) <= dist)]

def randomMutation(nucleotide):
    possible_mut = list(set('atcg') - set(nucleotide))
    
    return choice(possible_mut)

def analyze_cai(seq,cai_table = default_cai_table):
    seq = seq.lower()
    score = 0
    len_sq = 0
    for i in range(0, len(seq), 3):
        if seq[i:i + 3] in cai_table:
            score += log(cai_table[seq[i:i + 3]])
            len_sq += 1
    score /= len_sq
    return exp(score)