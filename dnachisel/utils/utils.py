"""Misc. functions using DnaChisel which can be useful in other programs."""

from ..DnaOptimizationProblem import DnaOptimizationProblem
from ..biotools import random_dna_sequence

def random_compatible_dna_sequence(sequence_length, constraints, probas=None,
                                   **kwargs):
    sequence = random_dna_sequence(sequence_length, probas=probas)
    canvas = DnaOptimizationProblem(sequence, constraints=constraints)
    canvas.solve_all_constraints_one_by_one(**kwargs)
    return canvas.sequence
