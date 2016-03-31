DNA Chisel
==========

DnaChisel is a Python library to modify the nucleotides of DNA sequences with respect to a set of
constraints and optimization objectives.

It can be used for many purposes, such as codon-optimizing the genes of a sequence
for a particular micro-organism, modifying a sequence to meet the constraints of
a DNA provider while preserving genes and other sensible patterns, or inserting
a pattern in a sequence using only synonymous mutations.

Example of use
---------------

In this basic example we optimize a sequence with respect to the following constraints and objectives:

- **Constraint 1:** The sequence should contain no restriction site for BsaI (GGTCTC).
- **Constraint 2:** The local GC content of every 50-nucleotide subsequence should be between 30% and 70%.
- **Objective 1:** The sequence's  GC content should be 40% (or as close as possible)

Here is the Python code to solve the problem with DnaChisel:
::
    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    canvas = DnaCanvas(
        sequence=random_dna_sequence(10000),
        constraints=[NoPatternConstraint(enzyme_pattern("BsaI")),
                     GCContentConstraint(0.3, 0.7, gc_window=50)],
        objectives = [GCContentObjective(0.4)]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    canvas.solve_all_constraints_one_by_one()
    canvas.maximize_all_objectives_one_by_one(max_random_iters=10000)

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    canvas.print_constraints_summary() # Will print success reports
    canvas.print_objectives_summary() # That will be good !

For a more complete and meaningful example, see also this other script, in which
a plasmid is codon-optimized and tweaked so as to verify constraints imposed by
a DNA synthesis company.

DnaChisel implements advanced constraints such as the preservation of coding
sequences,  or the inclusion or exclusion of advanced patterns, as well as
some common biological objectives (such as codon optimization, GC content), but it
is also very easy to implement new constraints and objectives.


Search strategies
-----------------

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive.
DnaChisel uses the following strategies to avoid exploring the whole search space:

- **Constraining of the mutation space:** no mutation can be done in segments of the sequence
  subject to a ``DoNotModify`` constraint, and in segments subject to an
  ``EnforceTranslation`` constraint only synonymous mutations of the codons are
  allowed.

- **Constraints solving before objective optimization**: DnaChisel currently enforces a
  resolution of problems in two steps: first solve the constraints and make sure
  that they all pass, then optimize the sequence with respect to the different
  objectives, while making sure that they all pass. While not to always yield
  optimal results, this heuristic gives generally very good results, and is more
  practical, as solving for the constraints first is generally very fast and directly
  informs on whether all constraints can be met.

- **Localized searches:** When DnaChisel finds that a constraint is not
  verified, and if the constraint breaches are localized on the
  sequence (for instance, a forbidden restriction site at a given location),
  then it will attempt to solve each breach separately
  by creating *localized* versions of the mutations space and constraints around
  the problematic region.
  It works the same for optimization objectives: localized objectives indicate
  on which segments of the sequence to focus the search.

- **A mix of exhaustive searches and random searches:** for each localized
  constraint problem, if the search space is small enough DnaChisel performs
  an exhaustive search (i.e. it tries every possible change of the sequence until
  all constraints are resolved), else DnaChisel performs a random search where
  if create random valid variations of the sequence until one meets all the
  constraints. The optimization of objectives functions in a similar way.


Installation
-------------

You can install DnaChisel through PIP
::
    sudo pip install dnachisel

Alternatively, you can unzip the sources in a folder and type
::
    sudo python setup.py install




Contribute
----------

DnaChisel is an open-source library originally written at the Edinburgh Genome Foundry by Zulko_.
It is released on Github under the MIT licence, everyone is welcome to contribute.
