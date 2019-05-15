"""Example of use of the AvoidPAttern specification"""

from dnachisel import (DnaOptimizationProblem, random_protein_sequence,
                       random_dna_sequence, Location, CodonOptimize,
                       reverse_translate, EnforceTranslation, biotools)
from python_codon_tables import get_codons_table
import numpy

def test_codon_optimize_bestcodon():
    numpy.random.seed(123)
    protein = random_protein_sequence(3000, seed=123)
    sequence = reverse_translate(protein)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species='e_coli')]
    )
    assert problem.objective_scores_sum() < 0
    problem.optimize()
    assert problem.objective_scores_sum() == 0

def test_codon_optimize_harmonized():
    numpy.random.seed(123)
    protein = random_protein_sequence(500, seed=123)
    sequence = reverse_translate(protein)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species='e_coli', mode='harmonized')]
    )
    assert (-700 < problem.objective_scores_sum() < -600)
    problem.optimize()
    assert (-350 < problem.objective_scores_sum())

def test_codon_optimize_harmonized_short_sequence():
    protein = "DDDKKKKKK"
    sequence = reverse_translate(protein)
    harmonization = CodonOptimize(species='b_subtilis', mode='harmonized')
    problem = DnaOptimizationProblem(
                sequence=sequence,
                constraints=[EnforceTranslation()],
                objectives=[harmonization]
            )
    assert problem.objective_scores_sum() < -7
    problem.optimize()
    assert -1 < problem.objective_scores_sum()



def test_codon_optimize_as_hard_constraint():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(2000, seed=123),
        constraints=[
            EnforceTranslation(location=Location(1000, 1300)),
            CodonOptimize(location=Location(1000, 1300), species='e_coli')
        ]
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()

def test_codon_optimize_with_custom_table():
    table = get_codons_table('b_subtilis')
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(1200, seed=123),
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(codon_usage_table=table)]
    )
    assert (problem.objective_scores_sum() < -10)
    problem.optimize()
    assert (problem.objective_scores_sum() == 0)
