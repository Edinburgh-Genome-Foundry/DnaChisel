from dnachisel import *

sequence = random_dna_sequence(300) + "ATTATTATT" + random_dna_sequence(300)

canvas = DNACanvas(
    sequence=sequence,
    constraints=[
        NoPatternConstraint(DNAPattern("ATTATTATT")),
        GCPercentConstraint(0.3, 0.7, gc_window=50)
    ],
    objectives = [GCPercentObjective(0.4)]
)

print ("\n\n=== Status before optimization ===")
canvas.print_constraints_summary(failed_only=True)
canvas.print_objectives_summary()

canvas.solve_all_constraints_one_by_one()
canvas.maximize_all_objectives_one_by_one()

print ("\n\n=== Status after optimization ===")
canvas.print_constraints_summary(failed_only=True)
canvas.print_objectives_summary()

print ("\n\nFinal sequence: %s" % canvas.sequence)
