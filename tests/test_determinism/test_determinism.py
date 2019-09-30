"""These tests enforce the contract that setting numpy.random.seed(xxx) will
determine the result of DNA chisel optimizations.

They can also be used to detect changes that break previous determinisms.
"""

import os
import dnachisel as dc
import numpy as np


def experiment_1(seed=123):
    """A DNA chisel optimization whose results produced the file
    test_determinism.py"""
    np.random.seed(seed)

    sequence = dc.reverse_translate(dc.random_protein_sequence(50))

    # MAXIMIZE THE GC CONTENT

    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[dc.EnforceTranslation()],
        objectives=[dc.EnforceGCContent(target=1)],
        logger=None,
    )
    problem.optimize()

    # BRING THE GC CONTENT BACK TO 50%

    problem = dc.DnaOptimizationProblem(
        sequence=problem.sequence,
        constraints=[dc.EnforceTranslation()],
        objectives=[dc.EnforceGCContent(target=0.5)],
        logger=None,
    )
    problem.optimize()

    return problem.sequence


def experiment_2(seed=123):
    np.random.seed(seed)
    sequence = dc.reverse_translate(dc.random_protein_sequence(1000))
    problem = dc.DnaOptimizationProblem(
        sequence=sequence,
        constraints=[
            dc.EnforceTranslation(),
            dc.EnforceGCContent(mini=0.4, maxi=0.6, window=50),
        ],
        objectives=[dc.CodonOptimize(species="e_coli")],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()
    return problem.sequence


def create_test_sequences_files():
    """File test_sequences.csv has been created by running this function"""
    for i, experiment in enumerate([experiment_1, experiment_2]):
        sequences_by_seed = [
            (seed, experiment(seed)) for seed in (123, 456, 789)
        ]
        with open("test_sequences_%d.csv" % (i + 1), "w") as f:
            f.write(
                "\n".join(
                    [
                        "%d,%s" % (seed, sequence)
                        for (seed, sequence) in sequences_by_seed
                    ]
                )
            )


def test_determinism_1():
    path = os.path.join("tests", "test_determinism", "test_sequences_1.csv")
    with open(path, "r") as f:
        for line in f.read().split("\n"):
            seed, sequence = line.split(",")
            assert sequence == experiment_1(int(seed))


def test_determinism_2():
    path = os.path.join("tests", "test_determinism", "test_sequences_2.csv")
    with open(path, "r") as f:
        for line in f.read().split("\n"):
            seed, sequence = line.split(",")
            assert sequence == experiment_2(int(seed))


if __name__ == "__main__":
    create_test_sequences_files()
