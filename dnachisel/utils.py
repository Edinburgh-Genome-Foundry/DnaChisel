from .DnaCanvas import DnaCanvas
from .biotools import random_dna_sequence

def random_compatible_dna_sequence(sequence_length, constraints, **kwargs):
    canvas = DnaCanvas(random_dna_sequence(sequence_length))
    canvas.solve_all_constraints_one_by_one(**kwargs)
    return canvas.sequence
