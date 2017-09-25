import matplotlib
matplotlib.use("Agg")

from dnachisel.reports import plot_sequence_manufacturability_difficulties
from dnachisel import random_dna_sequence


def test_plot_sequence_manufacturability_difficulties():
    sequence = random_dna_sequence(10000)
    ax = plot_sequence_manufacturability_difficulties(sequence)
