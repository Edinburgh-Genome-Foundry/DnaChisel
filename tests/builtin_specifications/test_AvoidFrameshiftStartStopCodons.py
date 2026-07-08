from Bio.Data import CodonTable

from dnachisel import (
    AvoidFrameshiftStartStopCodons,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceTranslation,
    Location,
    translate,
)


def has_frameshift_start_stop(sequence, genetic_table="Standard"):
    table = CodonTable.unambiguous_dna_by_name[genetic_table]
    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)
    for frame in (1, 2):
        for index in range(frame, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]
            if (codon in start_codons) or (codon in stop_codons):
                return True
    return False


def test_AvoidFrameshiftStartStopCodons():
    sequence = "AATGAACTGCAAGCTGAA"
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            EnforceTranslation(),
            AvoidFrameshiftStartStopCodons(),
        ],
        logger=None,
    )
    assert not problem.all_constraints_pass()
    problem.resolve_constraints()
    assert problem.all_constraints_pass()
    assert not has_frameshift_start_stop(problem.sequence)
    assert translate(problem.sequence) == translate(sequence)


def test_AvoidFrameshiftStartStopCodons_localized_preserves_frame_origin():
    sequence = "ATGTTAAGCGGC"
    problem = DnaOptimizationProblem(sequence=sequence, constraints=[], logger=None)
    constraint = AvoidFrameshiftStartStopCodons(location=(0, len(sequence)))

    localized = constraint.localized(Location(3, 6), problem=problem)
    evaluation = localized.evaluate(problem)

    assert not evaluation.passes
    assert [(location.start, location.end) for location in evaluation.locations] == [(4, 7)]


def test_AvoidFrameshiftStartStopCodons_constrains_codon_optimization():
    amino_acid_sequence = "MLSGKKGETIKK"
    sequence = "ATGCTGAGCGGCAAAAAAGGCGAAACCATTAAAAAA"
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            EnforceTranslation(translation=amino_acid_sequence),
            AvoidFrameshiftStartStopCodons(),
        ],
        objectives=[CodonOptimize(species="e_coli")],
        logger=None,
    )

    problem.resolve_constraints()
    problem.optimize()

    assert problem.all_constraints_pass(autopass=False)
    assert not has_frameshift_start_stop(problem.sequence)
    assert translate(problem.sequence) == amino_acid_sequence
