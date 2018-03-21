"""Misc. functions using DnaChisel which can be useful in other programs."""
import dnachisel as dc
from Bio import Restriction

def random_compatible_dna_sequence(sequence_length, constraints, probas=None,
                                   seed=None, max_random_iters=5000,
                                   logger='bar', **kwargs):

    sequence = dc.random_dna_sequence(
        sequence_length, probas=probas, seed=seed)
    problem = dc.DnaOptimizationProblem(sequence, constraints=constraints,
                                        logger=logger)
    problem.max_random_iters = max_random_iters
    problem.resolve_constraints(**kwargs)
    return problem.sequence

def make_restriction_part(part_length, left_overhang, right_overhang,
                          enzyme, forbidden_enzymes, assembly_enzyme='BsmBI'):
    l_left = len(left_overhang)
    l_right = len(right_overhang)
    left_overhang_location = (0, l_left)
    right_overhang_location = (l_left + part_length,
                               l_left + part_length + l_right)
    center_location = (l_left, l_left + part_length)
    core_sequence = (left_overhang + dc.random_dna_sequence(part_length)
                     + right_overhang)
    enforce_enzyme = dc.EnforcePatternOccurence(
        enzyme=enzyme, location=center_location)
    problem = dc.DnaOptimizationProblem(
        sequence=core_sequence,
        constraints=[
            dc.AvoidChanges(left_overhang_location),
            dc.AvoidChanges(right_overhang_location),
        ] + [enforce_enzyme] + [
            dc.AvoidPattern(enzyme=enzyme_name)
            for enzyme_name in forbidden_enzymes + [assembly_enzyme]
        ]
    )
    problem.resolve_constraints()
    core_sequence = dc.sequence_to_biopython_record(problem.sequence)
    for loc in [left_overhang_location, right_overhang_location]:
        dc.annotate_record(core_sequence, loc, 'overhang')
    site_location = enforce_enzyme.evaluate(problem).data['matches'][0]
    dc.annotate_record(core_sequence, site_location.to_tuple(), enzyme)
    assembly_site = Restriction.__dict__[assembly_enzyme].site
    flank = dc.sequence_to_biopython_record(assembly_site + 'A')
    dc.annotate_record(flank, label='flank')
    return flank + core_sequence + flank.reverse_complement()
