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



Reports
~~~~~~~~~~~~~~~~~~~~

.. mermaid::

   graph TD;
      o[optimization_with_report] <--> D[DnaOptimizationProblem]
      o -->|constraints clash| N[write_no_solution_error]
      o -->|no constraints clash| W[write_optimization_report]

.. automodule:: dnachisel.reports.optimization_reports
   :members:
.. automodule:: dnachisel.reports.constraints_reports
  :members:


Constraint reports
~~~~~~~~~~~~~~~~~~



DnaOptimizationProblem
----------------------

.. automodule:: dnachisel.DnaOptimizationProblem
   :members:

Specification
-------------

.. automodule:: dnachisel.specifications.Specification
   :members:

Built-in Specifications
-----------------------

.. automodule:: dnachisel.specifications.builtin_specifications
   :members:



SpecEvaluation
--------------


.. mermaid::

   graph TD;
       P[DnaOptimizationProblem] --> SE[SpecEvaluation]
       S[Specification] --> SE
       SE -->|grouped into| SES[SpecEvaluations]
       SES --- OES[Objective evaluations]
       SES --- CES[Constraints evaluations]

.. automodule:: dnachisel.specifications.SpecEvaluation
   :members:



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
