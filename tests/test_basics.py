import dnachisel as dc
import numpy

def test_avoid_pattern_gc_content():

    numpy.random.seed(123)

    problem = dc.DnaOptimizationProblem(
        sequence=dc.random_dna_sequence(10000),
        constraints=[
            dc.AvoidPattern(dc.enzyme_pattern("BsaI")),
            dc.EnforceGCContent(gc_min=0.3, gc_max=0.7, gc_window=50)
        ],
        objectives=[dc.EnforceGCContent(gc_objective=0.4)]
    )

    print ("\n\n=== Status before optimization ===")
    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

    problem.solve_constraints()
    assert problem.all_constraints_pass()
