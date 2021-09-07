"""Example of use of the AvoidPAttern specification"""
from io import StringIO
from dnachisel import (
    DnaOptimizationProblem,
    EnforceTranslation,
    random_dna_sequence,
    AvoidPattern,
    RepeatedKmerPattern,
    AvoidChanges,
    MotifPssmPattern,
    Location,
)
import numpy


def test_avoid_pattern_basics():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000, seed=123),
        constraints=[AvoidPattern("BsaI_site")],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_avoid_pattern_overlapping_locations():
    numpy.random.seed(123)
    problem = DnaOptimizationProblem(
        sequence="AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG",
        constraints=[AvoidPattern("NAN")],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert "A" not in problem.sequence[1:-1]


def test_avoid_repeated_small_kmers():
    problem = DnaOptimizationProblem(
        sequence="AGAAGAAGAAGAAGAAGATTTTTTTTTTTTTGGAGGAGGAGGACCCCCCCCCCCCGAGG",
        constraints=[AvoidPattern(RepeatedKmerPattern(3, 3))],
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_pattern_and_reverse():
    bsmbi = "CGTCTC"
    bsmbi_rev = "GAGACG"
    sequence = 10 * bsmbi + 25 * bsmbi_rev + 15 * bsmbi + 15 * bsmbi_rev
    problem = DnaOptimizationProblem(
        sequence,
        constraints=[AvoidPattern("BsmBI_site")],
        objectives=[AvoidChanges()],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    assert sum(problem.sequence_edits_as_array()) < 70


def test_AvoidPattern_on_strands():

    # Negative strand only
    sequence = "CATGCTATGC"
    problem = DnaOptimizationProblem(
        sequence, constraints=[AvoidPattern("CAT", location=(0, 10, -1))], logger=None,
    )
    problem.resolve_constraints()
    assert "CAT" in problem.sequence
    assert "ATG" not in problem.sequence

    # Negative strand only
    sequence = "CATGCTATGC"
    problem = DnaOptimizationProblem(
        sequence, constraints=[AvoidPattern("CAT", location=(0, 10, -1))], logger=None,
    )
    problem.resolve_constraints()
    assert "CAT" in problem.sequence
    assert "ATG" not in problem.sequence

    # Both strands
    sequence = "CATGCTATGC"
    problem = DnaOptimizationProblem(
        sequence, constraints=[AvoidPattern("CAT")], logger=None,
    )
    problem.resolve_constraints()
    assert "CAT" not in problem.sequence
    assert "ATG" not in problem.sequence


JASPAR_CONTENT = """
>MA0006.1	Ahr::Arnt
A  [     3      0      0      0      0      0 ]
C  [     8      0     23      0      0      0 ]
G  [     2     23      0     23      0     24 ]
T  [    11      1      1      1     24      0 ]
>MA0151.1	Arid3a
A  [    27      0      1     27     27     20 ]
C  [     0      0      9      0      0      0 ]
G  [     0      0      0      0      0      1 ]
T  [     0     27     17      0      0      6 ]
"""


def test_AvoidPattern_with_jaspar_motifs():
    stringio = StringIO(JASPAR_CONTENT)
    motif_patterns = MotifPssmPattern.list_from_file(
        stringio, file_format="jaspar", relative_threshold=0.9
    )
    problem = DnaOptimizationProblem(
        sequence="GGGGGGGGGGTGCGTGATTAAAGGGGG",
        constraints=[AvoidPattern(p) for p in motif_patterns],
    )
    assert 2 == len(problem.constraints_evaluations().all_locations())
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_AvoidPattern_with_regular_expression():
    sequence = (
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTG"
        "GTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGCGCGGC"
        "GAGGGCGAGGGCGATGCCACCAACGGCAAGCTGACCCTGAAGTTCATC"
    )
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation(), AvoidPattern(r"GGT(.*)GAT")],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()


def test_location_strand_gets_conserved():
    cst = AvoidPattern("AACAAAT", Location(4, 1624, -1))
    location = Location(9, 10)
    new_cst = cst.localized(location)
    assert new_cst.location.to_tuple() == (4, 16, -1)


def test_avoid_pattern_options():
    # Checks Github issue #53
    pattern = "C" * 4
    sequence = "A" * 6 + pattern

    # location=None
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidPattern(pattern, strand="from_location")],
        logger=None,
    )
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert pattern not in problem.sequence

    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidPattern(pattern, strand="both")],
        logger=None,
    )
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert pattern not in problem.sequence

    problem = DnaOptimizationProblem(
        sequence=sequence, constraints=[AvoidPattern(pattern, strand=-1)], logger=None,
    )
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert pattern in problem.sequence

    # location specified
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            AvoidPattern(pattern, location=Location(0, 10, -1), strand="from_location")
        ],
        logger=None,
    )
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    # sequence not changed because location strand is -1:
    assert pattern in problem.sequence

    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            AvoidPattern(pattern, location=Location(0, 10, -1), strand="both")
        ],
        logger=None,
    )
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    # sequence changed because strand option overwrites location:
    assert pattern not in problem.sequence

    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[AvoidPattern(pattern, location=Location(0, 10, 1), strand=-1)],
        logger=None,
    )
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    # sequence not changed because strand option overwrites location strand:
    assert pattern in problem.sequence
