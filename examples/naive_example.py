from dnachisel import *
import numpy

# We setup the randomizer to always get the same sequence
numpy.random.seed(123)

canvas = DnaCanvas(sequence=random_dna_sequence(10000),
                   constraints=[AvoidPattern(enzyme_pattern("BsaI")),
                                EnforceGCContent(gc_min=0.3, gc_max=0.7,
                                                 gc_window=50)],
                   objectives=[EnforceGCContent(gc_objective=0.4)])

print ("\n\n=== Status before optimization ===")
print(canvas.constraints_summary(failed_only=True))
print(canvas.objectives_summary())

print ("Now solving constraints...")
canvas.solve_all_constraints_one_by_one()
print ("Done. Now optimizing objectives...")
canvas.maximize_all_objectives_one_by_one(max_random_iters=10000)

print ("\n\n=== Status after optimization ===")
print (canvas.constraints_summary())
print (canvas.objectives_summary())
