"""Misc. functions using DnaChisel which can be useful in other programs."""
import dnachisel as dc
from Bio import Restriction

def random_compatible_dna_sequence(sequence_length, constraints, probas=None,
                                   seed=None, max_random_iters=5000,
                                   logger='bar', **kwargs):
    """Produce a random sequence complying to some specifications.

    Parameters
    ----------

    sequence_length
      You guessed it.
    
    probas
      Either None for a fully random initial sequence, or a dict of the form
      {"A": 0.5, "T": 0.2, ...} to tune initial nucleotide representation
    
    constraints
      List of all DnaChisel specifications that will be applied as constraints.
    
    seed
      Optional seed for the random number generator, for reproducibility.
      
    max_random_iters
      Maximum number of random tries per location solving for the solver.
    
    logger
      Either 'bar' or None (no logger) or any proglog logger.

    """
    sequence = dc.random_dna_sequence(
        sequence_length, probas=probas, seed=seed)
    problem = dc.DnaOptimizationProblem(sequence, constraints=constraints,
                                        logger=logger)
    problem.max_random_iters = max_random_iters
    problem.resolve_constraints(**kwargs)
    return problem.sequence