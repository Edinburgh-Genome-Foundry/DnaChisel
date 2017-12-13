"""Misc. functions using DnaChisel which can be useful in other programs."""

from ..DnaOptimizationProblem import DnaOptimizationProblem
from ..biotools import random_dna_sequence

def random_compatible_dna_sequence(sequence_length, constraints, probas=None,
                                   seed=None, max_random_iters=5000,
                                   logger='bar', **kwargs):

    sequence = random_dna_sequence(sequence_length, probas=probas, seed=seed)
    problem = DnaOptimizationProblem(sequence, constraints=constraints,
                                     logger=logger)
    problem.max_random_iters = max_random_iters
    problem.resolve_constraints(**kwargs)
    return problem.sequence
