"""This module provides useful biological tables as Python dictionnaries.
"""

from collections import defaultdict
import os


def reverse_table(table):
    """ Return a dictionary {v1: [k1a, k1b,...]} where k1a, k1b are all
    the keys of table such that table[k1]=v1.
    """
    new_table = defaultdict(lambda: [])
    for (k, v) in sorted(table.items()):
        new_table[v].append(k)
    return dict(new_table)


def dict_from_csv(filepath, sep=";"):
    """Read a CSV and store entries in a dict."""
    with open(filepath, "r") as f:
        return {
            line.split(sep)[0]: line.split(sep)[1].strip("\n")
            for line in f.readlines()
            if line not in ("", "\n")
        }


data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")

# TABLES DEFINITIONS START HERE

# COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C"}
with open(os.path.join(data_dir, "complements.csv"), "r") as f:
    COMPLEMENTS = dict([line.split(',')
                        for line in f.read().splitlines()])
OTHER_BASES = {
    "A": ["T", "G", "C"],
    "T": ["A", "G", "C"],
    "G": ["T", "A", "C"],
    "C": ["T", "G", "A"],
}

CODONS_TRANSLATIONS = dict_from_csv(
    os.path.join(data_dir, "codons_translations.csv"))
CODONS_SEQUENCES = reverse_table(CODONS_TRANSLATIONS)
CODONS_SYNONYMS = {
    codon: CODONS_SEQUENCES[CODONS_TRANSLATIONS[codon]]
    for codon in CODONS_TRANSLATIONS
}

AA_LONG_NAMES = dict_from_csv(os.path.join(data_dir, "aa_long_names.csv"))

NUCLEOTIDE_TO_REGEXPR = dict_from_csv(
    os.path.join(data_dir, "nucleotide_to_regexpr.csv"))



iupac_file = os.path.join(data_dir, "iupac_notation.csv")
IUPAC_NOTATION = {
    k: set(v)
    for k, v in dict_from_csv(iupac_file).items()
}
