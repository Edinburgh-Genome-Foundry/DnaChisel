DnaChisel
=========

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel
   :alt: Travis CI build status

(Documentation in progress)

DnaChisel (full documentation `here
<http://edinburgh-genome-foundry.github.io/DnaChisel/>`_) is a Python library to optimize
the nucleotides of DNA sequences with respect to a set of constraints and optimization objectives.


It can be used for codon-optimizing the genes of a sequence for a particular micro-organism,
modifying a sequence to meet the constraints of a DNA provider while preserving genes,
and other sensible patterns.

DnaChisel also provides much freedom to define optimization problems and model
new kinds of specifications, making it suitable for either automated sequence
design, or for complex custom design projects.

License = MIT
--------------

Bandwagon is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/Bandwagon>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

Installation
--------------

If you have pip installed, just type:

.. code:: python

    (sudo) pip install dnachisel

DnaChisel can be installed by unzipping the source code in one directory and using this command:

.. code:: python

    (sudo) python setup.py install


Example of use
---------------

In this basic example we optimize a sequence with respect to the following constraints and objectives:

- **Constraint 1:** The sequence should contain no restriction site for BsaI (GGTCTC).
- **Constraint 2:** The local GC content of every 50-nucleotide subsequence should be between 30% and 70%.
- **Objective 1:** The sequence's  GC content should be 40% (or as close as possible)

Here is the Python code to solve the problem with DnaChisel:

.. code:: python

    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000),
        constraints=[AvoidPattern(enzyme_pattern("BsaI")),
                     GCContentConstraint(0.3, 0.7, window=50)],
        objectives = [GCContentObjective(0.4)]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    canvas.resolve_constraints()
    canvas.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    print(problem.constraints_text_summary()) # Will print success reports
    print(problem.objectives_text_summary()) # That will be good !

For a more complete and meaningful example, see also this other script, in which
a plasmid is codon-optimized and tweaked so as to verify constraints imposed by
a DNA synthesis company.

DnaChisel implements advanced constraints such as the preservation of coding
sequences,  or the inclusion or exclusion of advanced patterns, as well as
some common biological objectives (such as codon optimization, GC content), but it
is also very easy to implement new constraints and objectives.

Installation
-------------

You can install DnaChisel through PIP
::
    sudo pip install dnachisel

Alternatively, you can unzip the sources in a folder and type
::
    sudo python setup.py install

To be able to generate plots and reports, run
::
    sudo pip install dna_features_viewer weasyprint

Contribute
----------

DnaChisel is an open-source software originally written at the `Edinburgh Genome Foundry
<http://www.genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (copyright Edinburgh Genome Foundry).
Everyone is welcome to contribute !
