


DNA Chisel Documentation
=========================


.. image:: _static/images/title.png
   :width: 500px
   :align: center


Dna Chisel is a Python library to modify the nucleotides of DNA sequences with respect to
constraints and optimization objectives. Potential use cases include:

- Codon-optimization of a sequence for a particular micro-organism.
- Modification of a sequence to meet the constraints of a DNA manufacturer.
- Insertion of patterns in a sequence through non-functionally altering modifications.


Example of use
---------------

Let us optimize a sequence with respect to the following objectives:

- The sequence should contain no restriction site for BsaI (GGTCTC).
- GC content should be between 30% and 70% on every 50-nucleotide subsequence.
- The sequence's  global GC content should be 40% (or as close as possible)
- there should be

Here is the Python code to solve the problem with DnaChisel:

.. code:python

    from dnachisel import *

    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(1500, seed=123),
        constraints=[
            AvoidPattern(enzyme="BsmBI"),
            AvoidChanges((100, 200)),
            EnforceTranslation((200, 1100, +1)),
            AvoidChanges((1100, 1200)),
            EnforcePattern('CCWGG', location=(600, 700))
        ],
        objectives = [
            CodonOptimize('e_coli', (200, 1100, +1)),
            EnforceGCContent(target=0.4, window=80)
        ]
    )

    # RESOLVE CONSTRAINTS AND OPTIMIZE
    problem.resolve_constraints()
    problem.optimize()

    # SAVE OPTIMIZED SEQUENCE
    problem.to_record("final_sequence.gb")

This prints the following result, indicating that all constraints pass in the end
and the objective has been (very well) optimized:

.. code:python

    ===> SUCCESS - all constraints evaluations pass
    NoPattern(GGTCTC (BsaI), None) Passed. Pattern not found !
    GCContent(min 0.30, max 0.70, gc_win 50, window None) Passed !


    ===> TOTAL OBJECTIVES SCORE: 0.00
    GCContentObj(0.40, global): scored -0.00E+00. GC content is 0.400 (0.400 wanted)



For a more complete and meaningful example, see also :ref:`this other script <plasmid-optimization>`,
in which a plasmid is codon-optimized and tweaked so as to verify constraints imposed by
a DNA synthesis company.

DnaChisel implements advanced constraints such as the preservation of coding
sequences,  or the inclusion or exclusion of advanced patterns, as well as
some common biological objectives (such as codon optimization, GC content), but it
is also very easy to implement new constraints and objectives.


Search strategies
-----------------

Long DNA sequences have a huge space of possible mutations
(just 20 nucleotides can form a trillion different sequences), therefore it is not
possible to solve a DNA optimization problem through an exhaustive search.
DnaChisel uses the following strategies to avoid exploring the whole search space:

- **Constraining of the mutation space:** Prior to searching, DnaChisel trims the
  possible mutations by analyzing the constraints of the problem. For instance
  a ``AvoidChanges(segment)`` constraint makes it impossible to mutate the nucleotides
  of the concerned DNA segment, and in segments subject to an
  ``EnforceTranslation`` constraint, only synonymous mutations of the codons are
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

.. raw:: html

       <a href="https://twitter.com/share" class="twitter-share-button"
       data-text="DnaChisel - A Python module for printing with living matter" data-size="large" data-hashtags="Bioprinting">Tweet
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
