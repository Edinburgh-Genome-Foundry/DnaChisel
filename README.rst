.. raw:: html

    <p align="center">
    <img alt="DNA Chisel Logo" title="DNA Chisel" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/title.png" width="450">
    <br /><br />
    </p>

DNA Chisel - a versatile sequence optimizer
============================================

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/DnaChisel/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/DnaChisel?branch=master


DNA Chisel (complete documentation `here <https://edinburgh-genome-foundry.github.io/DnaChisel/>`_) is a Python library for optimizing DNA sequences with respect to a set of constraints and optimization objectives. It comes with over 15 classes of sequence specifications which can be composed to codon-optimize genes, meet the constraints of a commercial DNA provider, avoid homologies between sequences, or all of this at once!

DNA Chisel also allows users to define their own specifications in Python, making the library suitable for a large range of automated sequence design applications, and complex custom design projects. It can be used as a Python library, a command-line interface, or a `web application <https://cuba.genomefoundry.org/sculpt_a_sequence>`_.

Example of use
---------------

Defining a problem via scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this basic example we generate a random sequence and optimize it so that

- It will be rid of BsaI sites.
- GC content will be between 30% and 70% on every 50bp window.
- The reading frame at position 500-1400 will be codon-optimized for *E. coli*.

Here is the code to achieve that:

.. code:: python

    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000),
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceGCContent(mini=0.3, maxi=0.7, window=50),
            EnforceTranslation((500, 1400))
        ],
        objectives=[CodonOptimize(species='e_coli', location=(500, 1400))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

DnaChisel implements advanced constraints such as the preservation of coding
sequences,  or the inclusion or exclusion of advanced patterns (see
`this page <https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html>`_
for an overview of available specifications), but it is also easy to implement
our own constraints and objectives as subclasses of ``dnachisel.Specification``.


Defining a problem via Genbank features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can also define a problem by annotating directly a Genbank as follows:

.. raw:: html

    <p align="center">
    <img alt="report" title="report" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/example_sequence_map.png" width="450">
    <br /><br />
    </p>

I this record:

- Constraints (colored in blue in the illustration) are features of type
  ``misc_feature`` with a prefix ``@`` followed
  by the name of the constraints and its parameters, which are the same as in
  python scripts.
- Optimization objectives (colored in yellow in the illustration) are features
  of type ``misc_feature`` with a prefix ``~`` followed by the name of the
  constraints and its parameters.

The file can be directly fed to the `web app <https://cuba.genomefoundry.org/sculpt_a_sequence>`_
or processed via the command line interface:

.. code:: bash

    # Output the result to "optimized_record.gb"
    dnachisel annotated_record.gb optimized_record.gb

Or via a Python script:

.. code:: python

    from dnachisel import DnaOptimizationProblem
    problem = DnaOptimizationProblem.from_record("my_record.gb")
    problem.optimize_with_report(target="report.zip")

By default, only the built-in specifications of DnaChisel can be used in the
annotations, however it is easy to add your own specifications to the Genbank
parser, and build applications supporting custom specifications on top of
DnaChisel.


Reports
~~~~~~~~

DnaChisel also implements features for verification and troubleshooting. For
instance by generating optimization reports:

.. code:: python

    problem.optimize_with_report(target="report.zip")

Here is an example of summary report:

.. raw:: html

    <p align="center">
    <img alt="report" title="report" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/report_screenshot.png" width="600">
    <br /><br />
    </p>




How it works
------------

DnaChisel hunts down every constraint breach and suboptimal region by
recreating local version of the problem around these regions. Each type of
constraint can be locally *reduced* and solved in its own way, to ensure fast
and reliable resolution.

Below is an animation of the algorithm in action:

.. raw:: html

    <p align="center">
    <img alt="DNA Chisel algorithm" title="DNA Chisel" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/dnachisel_algorithm.gif" width="800">
    <br />
    </p>

Installation
-------------

You can install DnaChisel through PIP:

.. code::

    sudo pip install dnachisel[reports]

The ``[reports]`` suffix will install some heavier libraries
(Matplotlib, PDF reports, sequenticon) for report generation,
you can omit it if you just want to use DNA chisel to edit sequences and
generate genbanks (for any interactive use, reports are highly recommended).

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

Optionally, also install Bowtie to be able to use ``AvoidMatches`` (which
removes short homologies with existing genomes). On Ubuntu:

.. code::

    sudo apt-get install ncbi-blast+


License = MIT
--------------

DnaChisel is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (¢ Edinburg Genome Foundry). Everyone is welcome to contribute !

More biology software
-----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Chisel is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
