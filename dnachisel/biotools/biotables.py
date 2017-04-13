"""This module provides useful biological tales as Python dictionnaries.
"""

from collections import defaultdict
import os

def reverse_table(table):
    """ Return a dictionary {v1: [k1a, k1b,...]} where k1a, k1b are all
    the keys of table such that table[k1]=v1.
    """
    new_table = defaultdict(lambda: [])
    for (k, v) in table.items():
        new_table[v].append(k)
    return dict(new_table)


data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
codons_usage_dir = os.path.join(data_dir, "codon_usage_tables")

CODON_USAGE_TABLES = {}
for fname in os.listdir(codons_usage_dir):
    if not fname.endswith(".csv"):
        continue
    organism = "_".join(fname[:-4].split("_")[:-1])
    with open(os.path.join(codons_usage_dir, fname), "r") as f:
        CODON_USAGE_TABLES[organism] = {
            line.split(";")[1]: float(line.split(";")[2].strip("\n"))
            for line in f.readlines()[1:]
        }


COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C"}

NUCLEOTIDE_TO_REGEXPR = {
    "A": "A",
    "B": "[BCGKSTY]",
    "C": "C",
    "D":  "[ADGKRTW]",
    "G": "G",
    "H": "[ACHMTWY]",
    "K": "[GKT]",
    "M": "[ACM]",
    "N": "[ABCDGHKMNRSTVWY]",
    "R": "[AGR]",
    "S": "[CGS]",
    "T": "T",
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


CODONS_TRANSLATIONS = {
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

CODONS_SEQUENCES = reverse_table(CODONS_TRANSLATIONS)

OTHER_BASES = {
    "A": ["T", "G", "C"],
    "T": ["A", "G", "C"],
    "G": ["T", "A", "C"],
    "C": ["T", "G", "A"],
}
