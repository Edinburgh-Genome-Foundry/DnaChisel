"""Misc. functions using DnaChisel which can be useful in other programs."""

from ..DnaOptimizationProblem import DnaOptimizationProblem
from ..biotools import random_dna_sequence

def random_compatible_dna_sequence(sequence_length, constraints, probas=None,
                                   seed=None, **kwargs):
    sequence = random_dna_sequence(sequence_length, probas=probas, seed=seed)
    canvas = DnaOptimizationProblem(sequence, constraints=constraints)
    canvas.resolve_constraints(**kwargs)
    return canvas.sequence
