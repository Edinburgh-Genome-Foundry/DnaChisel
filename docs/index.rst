
.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/DnaChisel/master/docs/_static/images/title.png
   :alt: [logo]
   :align: center
   :width: 500px

DNA Chisel
===========

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/DnaChisel
   :alt: Travis CI build status

DNA Chisel is a Python library to optimize the nucleotides of DNA sequences with respect to a set of constraints and optimization objectives. It can be used for codon-optimizing the genes of a sequence for a particular organism, modifying a sequence to meet the constraints of a DNA provider while preserving genes, and much more.

DNA Chisel comes with more than 15 types of optimizations and constraints and allows users to define new specifications in Python, making the library suitable for a large range of automated sequence design applications, or complex custom design projects.

Example of use
---------------

In this basic example we generate a random sequence and optimize it so that

- It will be rid of BsaI sites.
- GC content will be between 30% and 70% on every 50bp window.
- The reading frame at position 500-1300 will be codon-optimized for *E. coli*.

Here is the code to achieve that:

.. code:: python

    from dnachisel import *

    # DEFINE THE OPTIMIZATION PROBLEM

    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(10000),
        constraints=[AvoidPattern(enzyme_="BsaI"),
                     EnforceGCContent(mini=0.3, maxi=0.7, window=50)],
        objectives=[CodonOptimize(species='e_coli', location=(500, 1300))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

DnaChisel implements advanced constraints such as the preservation of coding
sequences,  or the inclusion or exclusion of advanced patterns, as well as
some common biological objectives (such as codon optimization, GC content), but it
is also very easy to implement new constraints and objectives.


Installation
-------------

You can install DnaChisel through PIP

.. code::

    sudo pip install dnachisel

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

To be able to generate plots and reports, run

.. code::

    sudo pip install dna_features_viewer weasyprint

License = MIT
--------------

DnaChisel is an open-source software originally written at the `Edinburgh Genome Foundry
<http://edinburgh-genome-foundry.github.io/home.html>`_ by `Zulko <https://github.com/Zulko>`_
and `released on Github <https://github.com/Edinburgh-Genome-Foundry/DnaChisel>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !


.. raw:: html

       <a href="https://twitter.com/share" class="twitter-share-button"
       data-text="DnaChisel - A Python library for DNA sequence design and optimization" data-size="large" data-hashtags="Bioprinting">Tweet
       </a>
       <script>!function(d,s,id){var js,fjs=d.getElementsByTagName(s)[0],p=/^http:/.test(d.location)?'http':'https';
       if(!d.getElementById(id)){js=d.createElement(s);js.id=id;js.src=p+'://platform.twitter.com/widgets.js';
       fjs.parentNode.insertBefore(js,fjs);}}(document, 'script', 'twitter-wjs');
       </script>
       <iframe src="http://ghbtns.com/github-btn.html?user=Edinburgh-Genome-Foundry&repo=dnachisel&type=watch&count=true&size=large"
       allowtransparency="true" frameborder="0" scrolling="0" width="152px" height="30px" margin-bottom="30px"></iframe>




.. toctree::
    :hidden:
    :maxdepth: 3

    self

.. toctree::
   :hidden:
   :caption: Reference
   :maxdepth: 3

   ref/ref

.. toctree::
   :caption: Examples

   examples/plasmid_optimization
   examples/non_unique_kmers_minimization
   examples/pattern_instertion


.. _Zulko: https://github.com/Zulko/
.. _Github: https://github.com/EdinburghGenomeFoundry/dnachisel
.. _PYPI: https://pypi.python.org/pypi/dnachisel
