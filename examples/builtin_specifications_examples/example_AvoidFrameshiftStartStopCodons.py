"""Demo for AvoidFrameshiftStartStopCodons."""

from Bio.Data import CodonTable

from dnachisel import (
    AvoidFrameshiftStartStopCodons,
    CodonOptimize,
    DnaOptimizationProblem,
    EnforceTranslation,
    translate,
)


def find_frameshift_start_stop(sequence, genetic_table="Standard"):
    table = CodonTable.unambiguous_dna_by_name[genetic_table]
    start_codons = set(table.start_codons)
    stop_codons = set(table.stop_codons)
    hits = []
    for frame in (1, 2):
        for index in range(frame, len(sequence) - 2, 3):
            codon = sequence[index : index + 3]
            if codon in start_codons or codon in stop_codons:
                hits.append((frame, index, codon))
    return hits


def run_without_spec(sequence, translation):
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[EnforceTranslation(translation=translation)],
        objectives=[CodonOptimize(species="e_coli")],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence


def run_with_spec(sequence):
    problem = DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            EnforceTranslation(),
            AvoidFrameshiftStartStopCodons(),
        ],
        logger=None,
    )
    problem.resolve_constraints()
    return problem.sequence


def main():
    sequence = "AATGAACTGCAAGCTGAA"
    print("Initial translation (+0):", translate(sequence))
    print("Frameshift hits:", find_frameshift_start_stop(sequence))

    fixed = run_with_spec(sequence)
    print("\nWith AvoidFrameshiftStartStopCodons")
    print("Optimized sequence:", fixed)
    print("Optimized translation (+0):", translate(fixed))
    print("Frameshift hits after optimization:", find_frameshift_start_stop(fixed))

    translation = "MKTSKQYQK*"
    start_seq = "ATGAAAACCAGCAAACAGTATCAGAAATAA"
    optimized = run_without_spec(start_seq, translation)
    print("\nWithout AvoidFrameshiftStartStopCodons")
    print("Optimized sequence:", optimized)
    print("Optimized translation (+0):", translate(optimized))
    print("Frameshift hits after optimization:", find_frameshift_start_stop(optimized))


if __name__ == "__main__":
    main()
