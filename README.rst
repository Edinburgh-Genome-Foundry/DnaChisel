.. raw:: html

    <p align="center">
    <img alt="DNA Chisel Logo" title="DNA Chisel" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/title.png" width="450">
    <br /><br />
    </p>

DNA Chisel - a versatile sequence optimizer
===========================================

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/DnaChisel/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/DnaChisel?branch=master


DNA Chisel (complete documentation `here <https://edinburgh-genome-foundry.github.io/DnaChisel/>`_)
is a Python library for optimizing DNA sequences with respect to a set of
constraints and optimization objectives. It can also be used via a command-line
interface, or a `web application <https://cuba.genomefoundry.org/sculpt_a_sequence>`_.

The library comes with over 15 classes of sequence specifications which can be
composed to, for instance, codon-optimize genes, meet the constraints of a
commercial DNA provider, avoid homologies between sequences, tune GC content,
or all of this at once! Users can also define their own specifications using
Python, making the library suitable for a large range of automated sequence
design applications, and complex custom design projects.



Usage
-----

Defining a problem via scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example below will generate a random sequence and optimize it so that:

- It will be rid of BsaI sites (on both strands).
- GC content will be between 30% and 70% on every 50bp window.
- The reading frame at position 500-1400 will be codon-optimized for *E. coli*.

.. code:: python

    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000),
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceGCContent(mini=0.3, maxi=0.7, window=50),
            EnforceTranslation(location=(500, 1400))
        ],
        objectives=[CodonOptimize(species='e_coli', location=(500, 1400))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

    final_sequence = problem.sequence  # string
    final_record = problem.to_record(with_sequence_edits=True)


Defining a problem via Genbank features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can also define a problem by annotating directly a Genbank as follows:

.. raw:: html

    <p align="center">
    <img alt="report" title="report" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/example_sequence.png" width="450">
    <br /><br />
    </p>


Note that constraints (colored in blue in the illustration) are features of type
``misc_feature`` with a prefix ``@`` followed by the name of the constraints
and its parameters, which are the same as in python scripts. Optimization
objectives (colored in yellow in the illustration) use prefix ``~``. See
`the Genbank API documentation <https://edinburgh-genome-foundry.github.io/DnaChisel/genbank/genbank_api.html>`_
for more details.

Genbank files with specification annotations can be directly fed to the
`web application <https://cuba.genomefoundry.org/sculpt_a_sequence>`_
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
~~~~~~~

DnaChisel also implements features for verification and troubleshooting. For
instance by generating optimization reports:

.. code:: python
    problem = DnaOptimizationProblem(...)
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
------------

DNA Chisel requires Python 3, and can be installed via a pip command:

.. code::
    sudo pip install dnachisel     # <= minimal install without reports support
    sudo pip install dnachisel[reports] # <= full install with all dependencies

The full installation using ``dnachisel[reports]`` downloads heavier libraries
(Matplotlib, PDF reports, sequenticon) for report generation, but is highly
recommended to use DNA Chisel interactively via Python scripts.

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

Optionally, also install Bowtie to be able to use ``AvoidMatches`` (which
removes short homologies with existing genomes). On Ubuntu:

.. code::

    sudo apt-get install bowtie


License = MIT
-------------

DnaChisel is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (Copyright 2017 Edinburgh Genome Foundry). Everyone is welcome to contribute!

More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Chisel is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.

Related projects
----------------

(If you would like to see a DNA Chisel-related project advertized here, please open
an issue or propose a PR)

- `Benchling <https://www.benchling.com/>`_ uses DNA Chisel as part of its sequence
  optimization pipeline according to `this webinar video <https://www.youtube.com/watch?v=oIcz5fQgtS8&t=865s>`_.
- `dnachisel-dtailor-mode <https://github.com/Lix1993/dnachisel_dtailor_mode>`_ brings
  features from `D-tailor <https://academic.oup.com/bioinformatics/article/30/8/1087/254801>`_
  to DNA Chisel, in particular for the generation of large collection of sequences
  covering the objectives fitness landscape (i.e. with sequences with are good at
  some objectives and bad at others, and vice versa).
