from dnachisel import *

canvas = DnaCanvas(
    sequence=random_dna_sequence(10000),
    constraints=[NoPatternConstraint(enzyme_pattern("BsaI")),
                 GCContentConstraint(0.3, 0.7, gc_window=50)],
    objectives = [GCContentObjective(0.4)]
)

print ("\n\n=== Status before optimization ===")
canvas.print_constraints_summary(failed_only=True)
canvas.print_objectives_summary()

canvas.solve_all_constraints_one_by_one()
canvas.maximize_all_objectives_one_by_one(max_random_iters=10000)

print ("\n\n=== Status after optimization ===")
canvas.print_constraints_summary()
canvas.print_objectives_summary()

# print ("\n\nFinal sequence: %s" % canvas.sequence)
