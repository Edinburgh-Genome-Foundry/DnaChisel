.. _reference ::

DnaChisel Reference manual
==========================

Reports
--------

.. toctree::
   :maxdepth: 2

   optimization_reports
   constraints_reports

Built-in Specifications
-----------------------

.. toctree::
   :maxdepth:2

   builtin_specifications/AvoidPattern

Core classes
---------------
.. mermaid::

   graph TD;
      sequence --> P[DnaOptimizationProblem]
      o[objectives...] -->P
      c[constraints...] -->P
      s[Specifications] --> o
      s --> o
      s --> c
      s --> c
      P -->|solve, optimize| O[Optimized problem.sequence]


 .. toctree::
    :maxdepth: 1

    DnaOptimizationProblem
    Location
    MutationSpace
    Specification
    SpecEvaluation






Constraint reports
~~~~~~~~~~~~~~~~~~





Biotools
--------

Features annotations
~~~~~~~~~~~~~~~~~~~~~

.. mermaid::

   graph TD;
     a[annotate_record] -->|used by| d[annotate_differences]
     a -->|used by| o[annotate_pattern_occurences]

.. automodule:: dnachisel.biotools.features_annotations
  :members:

 .. automodule:: dnachisel.biotools.biotools
   :members:

SequencePattern
---------------

.. automodule:: dnachisel.biotools.patterns
   :members:
