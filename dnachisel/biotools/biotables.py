"""This module provides useful biological tales as Python dictionnaries.
"""

from collections import defaultdict

def reverse_table(table):
    """ Return a dictionary {v1: [k1a, k1b,...]} where k1a, k1b are all
    the keys of table such that table[k1]=v1.
    """
    new_table = defaultdict(lambda: [])
    for (k, v) in table.items():
        new_table[v].append(k)
    return dict(new_table)


COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C"}

NUCLEOTIDE_TO_REGEXPR = {
    "A": "[A]",
    "B": "[BCGKSTY]",
    "C": "[C]",
    "D":  "[ADGKRTW]",
    "G": "[G]",
    "H": "[ACHMTWY]",
    "K": "[GKT]",
    "M": "[ACM]",
    "N": "[ABCDGHKMNRSTVWY]",
    "R": "[AGR]",
    "S": "[CGS]",
    "T": "[T]",
    "V": "[ACGMRSV]",
    "W": "[ATW]",
    "Y": "[CTY]",
}


AA_LONGNAMES = {
    "A": "Ala", "B": "Unk", "C": "Cys",
    "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile",
    "J": "Unk", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "O": "Unk",
    "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "U": "Unk",
    "V": "Val", "W": "Trp", "X": "Unk",
    "Y": "Tyr", "Z": "Unk", "*": "Stp"
}


CODON_TRANSLATIONS = {
    "TTT" : "F",   "TTC" : "F",   "TTA" : "L",   "TTG" : "L",
    "CTT" : "L",   "CTC" : "L",   "CTA" : "L",   "CTG" : "L",
    "ATT" : "I",   "ATC" : "I",   "ATA" : "I",   "ATG" : "M",
    "GTT" : "V",   "GTC" : "V",   "GTA" : "V",   "GTG" : "V",
    "TCT" : "S",   "TCC" : "S",   "TCA" : "S",   "TCG" : "S",
    "CCT" : "P",   "CCC" : "P",   "CCA" : "P",   "CCG" : "P",
    "ACT" : "T",   "ACC" : "T",   "ACA" : "T",   "ACG" : "T",
    "GCT" : "A",   "GCC" : "A",   "GCA" : "A",   "GCG" : "A",
    "TAT" : "Y",   "TAC" : "Y",   "TAA" : "*",   "TAG" : "*",
    "CAT" : "H",   "CAC" : "H",   "CAA" : "Q",   "CAG" : "Q",
    "AAT" : "N",   "AAC" : "N",   "AAA" : "K",   "AAG" : "K",
    "GAT" : "D",   "GAC" : "D",   "GAA" : "E",   "GAG" : "E",
    "TGT" : "C",   "TGC" : "C",   "TGA" : "*",   "TGG" : "W",
    "CGT" : "R",   "CGC" : "R",   "CGA" : "R",   "CGG" : "R",
    "AGT" : "S",   "AGC" : "S",   "AGA" : "R",   "AGG" : "R",
    "GGT" : "G",   "GGC" : "G",   "GGA" : "G",   "GGG" : "G"
}

CODONS_SEQUENCES = reverse_table(CODON_TRANSLATIONS)


# Values in this table are taken from the GeneDesign project:
# https://github.com/GeneDesign/GeneDesign

CODON_USAGE = {

    "S. cerevisiae": {


        "TTT" : 0.19,   "TTC" : 1.81,   "TTA" : 0.49,   "TTG" : 5.34,
        "CTT" : 0.02,   "CTC" : 0.00,   "CTA" : 0.15,   "CTG" : 0.02,
        "ATT" : 1.26,   "ATC" : 1.74,   "ATA" : 0.00,   "ATG" : 1.00,
        "GTT" : 2.07,   "GTC" : 1.91,   "GTA" : 0.00,   "GTG" : 0.02,
        "TCT" : 3.26,   "TCC" : 2.42,   "TCA" : 0.08,   "TCG" : 0.02,
        "CCT" : 0.21,   "CCC" : 0.02,   "CCA" : 3.77,   "CCG" : 0.00,
        "ACT" : 1.83,   "ACC" : 2.15,   "ACA" : 0.00,   "ACG" : 0.01,
        "GCT" : 3.09,   "GCC" : 0.89,   "GCA" : 0.03,   "GCG" : 0.00,
        "TAT" : 0.06,   "TAC" : 1.94,   "TAA" : 1.00,   "TAG" : 0.00,
        "CAT" : 0.32,   "CAC" : 1.68,   "CAA" : 1.98,   "CAG" : 0.02,
        "AAT" : 0.06,   "AAC" : 1.94,   "AAA" : 0.16,   "AAG" : 1.84,
        "GAT" : 0.70,   "GAC" : 1.30,   "GAA" : 1.98,   "GAG" : 0.02,
        "TGT" : 1.80,   "TGC" : 0.20,   "TGA" : 0.00,   "TGG" : 1.00,
        "CGT" : 0.63,   "CGC" : 0.00,   "CGA" : 0.00,   "CGG" : 0.00,
        "AGT" : 0.06,   "AGC" : 0.16,   "AGA" : 5.37,   "AGG" : 0.00,
        "GGT" : 3.92,   "GGC" : 0.06,   "GGA" : 0.00,   "GGG" : 0.02
    },

    "E. coli": {

        "TTT" : 0.34,   "TTC" : 1.66,   "TTA" : 0.06,   "TTG" : 0.07,
        "CTT" : 0.13,   "CTC" : 0.17,   "CTA" : 0.04,   "CTG" : 5.54,
        "ATT" : 0.48,   "ATC" : 2.51,   "ATA" : 0.01,   "ATG" : 1.00,
        "GTT" : 2.41,   "GTC" : 0.08,   "GTA" : 1.12,   "GTG" : 0.40,
        "TCT" : 2.81,   "TCC" : 2.07,   "TCA" : 0.06,   "TCG" : 0.00,
        "CCT" : 0.15,   "CCC" : 0.02,   "CCA" : 0.42,   "CCG" : 3.41,
        "ACT" : 1.87,   "ACC" : 1.91,   "ACA" : 0.10,   "ACG" : 0.12,
        "GCT" : 2.02,   "GCC" : 0.18,   "GCA" : 1.09,   "GCG" : 0.71,
        "TAT" : 0.38,   "TAC" : 1.63,   "TAA" : 1.00,   "TAG" : 1.00,
        "CAT" : 0.45,   "CAC" : 1.55,   "CAA" : 0.12,   "CAG" : 1.88,
        "AAT" : 0.02,   "AAC" : 1.98,   "AAA" : 1.63,   "AAG" : 0.37,
        "GAT" : 0.51,   "GAC" : 1.49,   "GAA" : 1.64,   "GAG" : 0.36,
        "TGT" : 0.60,   "TGC" : 1.40,   "TGA" : 1.00,   "TGG" : 1.00,
        "CGT" : 4.47,   "CGC" : 1.53,   "CGA" : 0.00,   "CGG" : 0.00,
        "AGT" : 0.13,   "AGC" : 0.93,   "AGA" : 0.00,   "AGG" : 0.00,
        "GGT" : 2.27,   "GGC" : 1.68,   "GGA" : 0.00,   "GGG" : 0.04
    },

    "H. Sapiens": {

        "TTT" : 0.27,   "TTC" : 1.73,   "TTA" : 0.05,   "TTG" : 0.31,
        "CTT" : 0.20,   "CTC" : 1.42,   "CTA" : 0.15,   "CTG" : 3.88,
        "ATT" : 0.45,   "ATC" : 2.43,   "ATA" : 0.12,   "ATG" : 1.00,
        "GTT" : 0.09,   "GTC" : 1.03,   "GTA" : 0.11,   "GTG" : 2.78,
        "TCT" : 0.45,   "TCC" : 2.09,   "TCA" : 0.26,   "TCG" : 0.68,
        "CCT" : 0.58,   "CCC" : 2.02,   "CCA" : 0.36,   "CCG" : 1.04,
        "ACT" : 0.36,   "ACC" : 2.37,   "ACA" : 0.36,   "ACG" : 0.92,
        "GCT" : 0.45,   "GCC" : 2.38,   "GCA" : 0.36,   "GCG" : 0.82,
        "TAT" : 0.34,   "TAC" : 1.66,   "TAA" : 1.00,   "TAG" : 1.00,
        "CAT" : 0.30,   "CAC" : 1.70,   "CAA" : 0.21,   "CAG" : 1.79,
        "AAT" : 0.33,   "AAC" : 1.67,   "AAA" : 0.34,   "AAG" : 1.66,
        "GAT" : 0.36,   "GAC" : 1.64,   "GAA" : 0.26,   "GAG" : 1.74,
        "TGT" : 0.42,   "TGC" : 1.58,   "TGA" : 1.00,   "TGG" : 1.00,
        "CGT" : 0.38,   "CGC" : 2.72,   "CGA" : 0.31,   "CGG" : 1.53,
        "AGT" : 0.31,   "AGC" : 2.22,   "AGA" : 0.22,   "AGG" : 0.84,
        "GGT" : 0.34,   "GGC" : 2.32,   "GGA" : 0.29,   "GGG" : 1.05
    },

    "C. elegans": {

        "TTT" : 0.72,   "TTC" : 1.28,   "TTA" : 0.45,   "TTG" : 1.35,
        "CTT" : 1.86,   "CTC" : 1.38,   "CTA" : 0.34,   "CTG" : 0.63,
        "ATT" : 1.52,   "ATC" : 1.23,   "ATA" : 0.25,   "ATG" : 1.00,
        "GTT" : 1.67,   "GTC" : 1.11,   "GTA" : 0.52,   "GTG" : 0.70,
        "TCT" : 1.47,   "TCC" : 0.98,   "TCA" : 1.44,   "TCG" : 0.83,
        "CCT" : 0.52,   "CCC" : 0.23,   "CCA" : 2.75,   "CCG" : 0.51,
        "ACT" : 1.34,   "ACC" : 1.02,   "ACA" : 1.15,   "ACG" : 0.49,
        "GCT" : 1.64,   "GCC" : 1.06,   "GCA" : 0.99,   "GCG" : 0.31,
        "TAT" : 0.97,   "TAC" : 1.03,   "TAA" : 1.00,   "TAG" : 1.00,
        "CAT" : 1.13,   "CAC" : 0.87,   "CAA" : 1.39,   "CAG" : 0.61,
        "AAT" : 1.10,   "AAC" : 0.90,   "AAA" : 0.84,   "AAG" : 1.16,
        "GAT" : 1.36,   "GAC" : 0.64,   "GAA" : 1.15,   "GAG" : 0.85,
        "TGT" : 1.14,   "TGC" : 0.86,   "TGA" : 1.00,   "TGG" : 1.00,
        "CGT" : 1.84,   "CGC" : 0.73,   "CGA" : 1.07,   "CGG" : 0.31,
        "AGT" : 0.76,   "AGC" : 0.52,   "AGA" : 1.79,   "AGG" : 0.26,
        "GGT" : 0.70,   "GGC" : 0.28,   "GGA" : 2.85,   "GGG" : 0.16
    },

    "D. melanogaster": {

        "TTT" : 0.12,   "TTC" : 1.88,   "TTA" : 0.03,   "TTG" : 0.69,
        "CTT" : 0.25,   "CTC" : 0.72,   "CTA" : 0.06,   "CTG" : 4.25,
        "ATT" : 0.74,   "ATC" : 2.26,   "ATA" : 0.00,   "ATG" : 1.00,
        "GTT" : 0.56,   "GTC" : 1.59,   "GTA" : 0.06,   "GTG" : 1.79,
        "TCT" : 0.87,   "TCC" : 2.74,   "TCA" : 0.04,   "TCG" : 1.17,
        "CCT" : 0.42,   "CCC" : 2.73,   "CCA" : 0.62,   "CCG" : 0.23,
        "ACT" : 0.65,   "ACC" : 3.04,   "ACA" : 0.10,   "ACG" : 0.21,
        "GCT" : 0.95,   "GCC" : 2.82,   "GCA" : 0.09,   "GCG" : 0.14,
        "TAT" : 0.23,   "TAC" : 1.77,   "TAA" : 1.00,   "TAG" : .00,
        "CAT" : 0.29,   "CAC" : 1.71,   "CAA" : 0.03,   "CAG" : 1.97,
        "AAT" : 0.13,   "AAC" : 1.87,   "AAA" : 0.06,   "AAG" : 1.94,
        "GAT" : 0.90,   "GAC" : 1.10,   "GAA" : 0.19,   "GAG" : 1.81,
        "TGT" : 0.07,   "TGC" : 1.93,   "TGA" : 1.00,   "TGG" : 1.00,
        "CGT" : 2.65,   "CGC" : 3.07,   "CGA" : 0.07,   "CGG" : 0.00,
        "AGT" : 0.04,   "AGC" : 1.13,   "AGA" : 0.00,   "AGG" : 0.21,
        "GGT" : 1.34,   "GGC" : 1.66,   "GGA" : 0.99,   "GGG" : 0.00
    },

    "B. subtilis": {

        "TTT" : 0.70,   "TTC" : 1.30,   "TTA" : 2.71,   "TTG" : 0.00,
        "CTT" : 2.13,   "CTC" : 0.00,   "CTA" : 1.16,   "CTG" : 0.00,
        "ATT" : 0.91,   "ATC" : 1.96,   "ATA" : 0.13,   "ATG" : 1.00,
        "GTT" : 1.88,   "GTC" : 0.25,   "GTA" : 1.38,   "GTG" : 0.50,
        "TCT" : 3.45,   "TCC" : 0.00,   "TCA" : 1.50,   "TCG" : 0.00,
        "CCT" : 2.29,   "CCC" : 0.00,   "CCA" : 1.14,   "CCG" : 0.57,
        "ACT" : 2.21,   "ACC" : 0.00,   "ACA" : 1.38,   "ACG" : 0.41,
        "GCT" : 2.94,   "GCC" : 0.08,   "GCA" : 0.60,   "GCG" : 0.38,
        "TAT" : 0.50,   "TAC" : 1.50,   "TAA" : 0.00,   "TAG" : 0.00,
        "CAT" : 2.00,   "CAC" : 0.00,   "CAA" : 1.71,   "CAG" : 0.29,
        "AAT" : 0.47,   "AAC" : 1.53,   "AAA" : 1.83,   "AAG" : 0.17,
        "GAT" : 0.53,   "GAC" : 1.47,   "GAA" : 1.40,   "GAG" : 0.60,
        "TGT" : 0.00,   "TGC" : 2.00,   "TGA" : 0.00,   "TGG" : 1.00,
        "CGT" : 3.11,   "CGC" : 1.78,   "CGA" : 1.00,   "CGG" : 0.00,
        "AGT" : 0.45,   "AGC" : 0.60,   "AGA" : 1.11,   "AGG" : 0.00,
        "GGT" : 1.38,   "GGC" : 0.97,   "GGA" : 1.66,   "GGG" : 0.00
    }
}

OTHER_BASES = {
    "A": ["T", "G", "C"],
    "T": ["A", "G", "C"],
    "G": ["T", "A", "C"],
    "C": ["T", "G", "A"],
}
