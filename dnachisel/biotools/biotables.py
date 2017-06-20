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


def dict_from_csv(filepath, sep=";"):
    with open(filepath, "r") as f:
        return {
            line.split(sep)[0]: line.split(sep)[1].strip("\n")
            for line in f.readlines()
            if line not in ("", "\n")
        }

data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
codons_usage_dir = os.path.join(data_dir, "codon_usage_tables")

# TABLES DEFINITIONS START HERE

COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C"}
OTHER_BASES = {
    "A": ["T", "G", "C"],
    "T": ["A", "G", "C"],
    "G": ["T", "A", "C"],
    "C": ["T", "G", "A"],
}

CODON_USAGE_TABLES = {}
for fname in os.listdir(codons_usage_dir):
    if not fname.endswith(".csv"):
        continue
    organism = "_".join(fname[:-4].split("_")[:-1])
    with open(os.path.join(codons_usage_dir, fname), "r") as f:
        CODON_USAGE_TABLES[organism] = {
            line.split(";")[1]: float(line.split(";")[2].strip("\n"))
            for line in f.readlines()[1:]
            if line not in ("", "\n")
        }

AA_LONG_NAMES = dict_from_csv(os.path.join(data_dir, "aa_long_names.csv"))

NUCLEOTIDE_TO_REGEXPR = dict_from_csv(
    os.path.join(data_dir, "nucleotide_to_regexpr.csv"))

CODONS_TRANSLATIONS = dict_from_csv(
    os.path.join(data_dir, "codons_translations.csv"))

CODONS_SEQUENCES = reverse_table(CODONS_TRANSLATIONS)
