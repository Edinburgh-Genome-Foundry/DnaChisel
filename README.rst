DNA Chisel
==========

DNAChisel is a Python library to modify the nucleotides of DNA sequences with respect to a set of
constraints and optimization objectives.

Example
-------

In the following example we consider a randomly-generated generated sequence, which
we optimize to have a GC content of 40%, under the constraints that (1) the sequence
should contain no restriction site for BsaI (GGTCTC), and (2) the local GC content on every
50-nucleotide window should remain between 30% and 70%.

.. code:: python
    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    canvas = DNACanvas(
        sequence=random_dna_sequence(10000),
        constraints=[NoPatternConstraint(enzyme_pattern("BsaI")),
                     GCContentConstraint(0.3, 0.7, gc_window=50)],
        objectives = [GCContentObjective(0.4)]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    canvas.solve_all_constraints_one_by_one()
    canvas.maximize_all_objectives_one_by_one()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    canvas.print_constraints_summary() # Will print success reports
    canvas.print_objectives_summary() # That will be good !

See also this other script for a more complete and meaningful example, featuring more advanced
constraints (preserving coding sequences, preserving arbitrary sequences,
avoiding homopolymers, etc.) and optimizations (such as codon optimization, or
conservative optimization). DNAChisel provides a framework in which it is also
very easy to implement new constraints and objectives.


Search strategies
-----------------

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive.
DNAChisel uses the following strategies to avoid exploring the whole search space:

- **Constraining of the mutation space:* no mutation can be done in segments of the sequence
  subject to a ``DoNotModify`` constraint, and in segments subject to an
  ``EnforceTranslation`` constraint only synonymous mutations of the codons are
  allowed.

- **Constraints solving before objective optimization**: DNAChisel currently enforces a
  resolution of problems in two steps: first solve the constraints and make sure
  that they all pass, then optimize the sequence with respect to the different
  objectives, while making sure that they all pass. While not to always yield
  optimal results, this heuristic gives generally very good results, and is more
  practical, as solving for the constraints first is generally very fast and directly
  informs on whether all constraints can be met.

- **Localized searches:** When DNAChisel finds that a constraint is not
  verified, and if the constraint breaches are localized on the
  sequence (for instance, a forbidden restriction site at a given location),
  then it will attempt to solve each breach separately
  by creating *localized* versions of the mutations space and constraints around
  the problematic region.
  It works the same for optimization objectives: localized objectives indicate
  on which segments of the sequence to focus the search.

- **A mix of exhaustive searches and random searches:** for each localized
  constraint problem, if the search space is small enough DNAChisel performs
  an exhaustive search (i.e. it tries every possible change of the sequence until
  all constraints are resolved), else DNAChisel performs a random search where
  if create random valid variations of the sequence until one meets all the
  constraints. The optimization of objectives functions in a similar way.



  Installation
  -------------

  You can install DNAChisel through PIP
  ::
      sudo pip install dnachisel

  Alternatively, you can unzip the sources in a folder and type
  ::
      sudo python setup.py install




Contribute
----------

DNAChisel is an open-source library originally written at the Edinburgh Genome Foundry by Zulko_.
It is released on Github under the MIT licence, everyone is welcome to contribute.

All examples
-------------

Reference
---------
