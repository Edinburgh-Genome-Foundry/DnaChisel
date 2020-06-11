"""Example of use of the AvoidPAttern specification"""

from dnachisel import (
    DnaOptimizationProblem,
    random_dna_sequence,
    reverse_complement,
    AvoidHairpins,
)
import numpy


def test_avoid_hairpin_basics():
    numpy.random.seed(123)
    random_sequences = [random_dna_sequence(30) for i in range(10)]

    full_sequence = "".join(
        [
            seq
            for sequence in random_sequences
            for seq in (
                random_dna_sequence(50),
                sequence,
                random_dna_sequence(50),
                reverse_complement(sequence),
                random_dna_sequence(50),
            )
        ]
    )

    problem = DnaOptimizationProblem(
        full_sequence, constraints=[AvoidHairpins()], logger=None
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()

def test_avoid_hairpins_on_extremities():
    # see https://github.com/Edinburgh-Genome-Foundry/DnaChisel/issues/37
    problem = DnaOptimizationProblem(
        sequence="attcaatgggggggggggggggggggggggggtagccta",
        constraints=[AvoidHairpins(stem_size=3, hairpin_window=8)] 
    )
    evaluation = problem.constraints_evaluations().evaluations[0]
    assert str(evaluation.locations) == "[0-6, 32-39]"
    problem.resolve_constraints()
    assert problem.all_constraints_pass()