DNA Chisel
==========

DNA Chisel (complete documentation `here <https://edinburgh-genome-foundry.github.io/DnaChisel/>`_)
is a Python library for optimizing DNA sequences with respect to a set of
constraints and optimization objectives. It comes with over 15 classes of
sequence specifications which can be composed to, for instance, codon-optimize
genes, meet the  constraints of a commercial DNA provider, avoid homologies
between sequences, tune GC content, or all of this at once!

DNA Chisel also allows users to define their own specifications in Python,
making the library suitable for a large range of automated sequence design
applications, and complex custom design projects. It can be used as a Python
library, a command-line interface, or a `web application <https://cuba.genomefoundry.org/sculpt_a_sequence>`_.


Example of use
---------------

.. code:: python

    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    some_sequence = random_dna_sequence(10000)
    problem = DnaOptimizationProblem(
        sequence=some_sequence,
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

Alternatively, DNA Chisel lets you define problems by annotating a Genbank file.
You can also define a problem by annotating directly a Genbank as follows:

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/example_sequence.png
   :alt: [logo]
   :align: center
   :width: 450px

See `this page <https://edinburgh-genome-foundry.github.io/DnaChisel/ref/builtin_specifications.html>`_
for an overview of available specifications.

Infos
-----

**PIP installation:**

.. code:: bash

  pip install dnachisel[reports]

(you can omit the ``[reports]`` suffix if you intend to use dnachisel only
for sequence optimization, without generating figures or PDF reports)

**Web documentation:**

`<https://edinburgh-genome-foundry.github.io/DnaChisel/>`_

**Github Page**

`<https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_

**Live demo**

`<http://cuba.genomefoundry.org/sculpt_a_sequence>`_

**License:** MIT, Copyright Edinburgh Genome Foundry

More biology software
-----------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

DNA Chisel is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
