DNA Chisel
==========

DNA Chisel is a Python library to optimize the nucleotides of DNA sequences with respect to a set of constraints and optimization objectives. It can be used for codon-optimizing the genes of a sequence for a particular organism, modifying a sequence to meet the constraints of a DNA provider while preserving genes, and much more.


DNA Chisel comes with more than 15 types of optimizations and constraints and allows users to define new specifications in Python, making the library suitable for a large range of automated sequence design applications, or complex custom design projects.

Example of use
---------------

.. code:: python

    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    random_sequence = random_dna_sequence(10000)
    problem = DnaOptimizationProblem(
        sequence=random_sequence,
        constraints=[AvoidPattern("BsaI_site"),
                     EnforceGCContent(mini=0.3, maxi=0.7, window=50)],
        objectives=[CodonOptimize(species='e_coli', location=(500, 1300))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

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
