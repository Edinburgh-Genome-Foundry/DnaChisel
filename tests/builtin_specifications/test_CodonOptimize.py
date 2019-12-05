"""Example of use of the AvoidPAttern specification"""

from dnachisel import (
    DnaOptimizationProblem,
    random_protein_sequence,
    random_dna_sequence,
    Location,
    CodonOptimize,
    reverse_translate,
    EnforceTranslation,
    biotools,
)
from python_codon_tables import get_codons_table
import numpy


def test_codon_optimize_bestcodon():
    numpy.random.seed(123)
    protein = random_protein_sequence(3000, seed=123)
    sequence = reverse_translate(protein)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(species="e_coli")],
        logger=None,
    )
    assert problem.objective_scores_sum() < 0
    problem.optimize()
    assert problem.objective_scores_sum() == 0


def test_codon_optimize_match_usage_random_sequence():
    numpy.random.seed(123)
    protein = random_protein_sequence(500, seed=123)
    sequence = reverse_translate(protein)
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[
            CodonOptimize(species="e_coli", method="match_codon_usage")
        ],
        logger=None,
    )
    assert -600 < problem.objective_scores_sum() < -550
    problem.optimize()
    print (problem.objective_scores_sum())
    assert -17 < problem.objective_scores_sum()


def test_codon_optimize_match_usage_gfp_sequence():
    sequence = (
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTG"
        "GTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGCGCGGC"
        "GAGGGCGAGGGCGATGCCACCAACGGCAAGCTGACCCTGAAGTTCATC"
    )
    spec = CodonOptimize(species="s_cerevisiae", method="match_codon_usage")
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[spec],
        logger=None,
    )
    assert problem.objective_scores_sum() < -61
    problem.optimize()
    assert problem.objective_scores_sum() > -16

    # Just for coverage, we run the compare_frequency function in text mode
    spec = problem.objectives[0]
    codons = spec.get_codons(problem)
    print(spec.compare_frequencies(codons, text_mode=True))


def test_codon_optimize_match_usage_short_sequence():
    numpy.random.seed(123)
    protein = "DDDKKKKKK"
    sequence = reverse_translate(protein)
    harmonization = CodonOptimize(
        species="b_subtilis", method="match_codon_usage"
    )
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[harmonization],
        logger=None,
    )
    assert problem.objective_scores_sum() < -5.5
    problem.optimize()
    assert -0.6 < problem.objective_scores_sum()
    print(problem.objective_scores_sum())
    assert problem.sequence == "GATGATGACAAGAAAAAGAAAAAAAAA"


def test_codon_optimize_harmonize_rca_short_sequence():
    protein = random_protein_sequence(500, seed=123)
    sequence = reverse_translate(protein)
    harmonization = CodonOptimize(
        species="h_sapiens", original_species="e_coli", method="harmonize_rca"
    )
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation()],
        objectives=[harmonization],
        logger=None,
    )
    assert problem.objective_scores_sum() < -123
    problem.optimize()
    assert -74 < problem.objective_scores_sum()


def test_codon_optimize_as_hard_constraint():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(2000, seed=123),
        constraints=[
            EnforceTranslation(location=Location(1000, 1300)),
            CodonOptimize(location=Location(1000, 1300), species="e_coli"),
        ],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_codon_optimize_with_custom_table():
    table = get_codons_table("b_subtilis")
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(1200, seed=123),
        constraints=[EnforceTranslation()],
        objectives=[CodonOptimize(codon_usage_table=table)],
        logger=None,
    )
    assert problem.objective_scores_sum() < -10
    problem.optimize()
    assert problem.objective_scores_sum() == 0
