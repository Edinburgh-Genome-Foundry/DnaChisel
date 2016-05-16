from .DnaCanvas import DnaCanvas
from .biotools import random_dna_sequence

def random_compatible_dna_sequence(sequence_length, constraints, probas=None,
                                   **kwargs):
    sequence = random_dna_sequence(sequence_length, probas=probas)
    canvas = DnaCanvas(sequence, constraints=constraints)
    canvas.solve_all_constraints_one_by_one(**kwargs)
    return canvas.sequence
