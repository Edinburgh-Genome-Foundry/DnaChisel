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
codons_usage_dir = os.path.join(data_dir, "codon_usage_tables")

# TABLES DEFINITIONS START HERE

COMPLEMENTS = {"A": "T", "T": "A", "C": "G", "G": "C"}
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

CODON_USAGE_TABLES = {}
for fname in os.listdir(codons_usage_dir):
    if not fname.endswith(".csv"):
        continue
    organism = "_".join(fname[:-4].split("_")[:-1])
    with open(os.path.join(codons_usage_dir, fname), "r") as f:
        CODON_USAGE_TABLES[organism] = table = {
            line.split(";")[1]: float(line.split(";")[2].strip("\n"))
            for line in f.readlines()[1:]
            if line not in ("", "\n")
        }
        for (codon, usage) in list(table.items()):
            table[codon.replace('U', 'T')] = usage
        table['best_frequencies'] = {
            aa: max([table[codon] for codon in CODONS_SEQUENCES[aa]])
            for aa in CODONS_SEQUENCES
        }

CODON_USAGE_BY_AA = {
    species: {
        aa: {
            codon: species_data[codon]
            for codon in codons
        }
        for aa, codons in CODONS_SEQUENCES.items()
    }
    for species, species_data in CODON_USAGE_TABLES.items()
}


AA_LONG_NAMES = dict_from_csv(os.path.join(data_dir, "aa_long_names.csv"))

NUCLEOTIDE_TO_REGEXPR = dict_from_csv(
    os.path.join(data_dir, "nucleotide_to_regexpr.csv"))



iupac_file = os.path.join(data_dir, "iupac_notation.csv")
IUPAC_NOTATION = {
    k: set(v)
    for k, v in dict_from_csv(iupac_file).items()
}
