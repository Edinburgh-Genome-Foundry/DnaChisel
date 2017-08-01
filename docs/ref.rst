.. reference ::

DnaChisel Reference manual
==========================

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

DnaOptimizationProblem
----------------------

.. autoclass:: dnachisel.DnaOptimizationProblem
   :members:

Specification
-------------


Built-in Specifications
-----------------------


SpecEvaluation
--------------

   graph TD;
       P[DnaOptimizationProblem] --> SE[SpecEvaluation]
       S[Specification] --> SE
       SE -->|grouped into| SES[SpecEvaluations]
       SES --- OES[Objective evaluations]
       SES --- CES[Constraints evaluations]

.. automodule:: dnachisel.SpecEvaluation
   :members:

Objectives
-----------

.. automodule:: dnachisel.objectives.Objective
   :members:

.. automodule:: dnachisel.objectives.objectives
   :members:



Patterns
---------

.. automodule:: dnachisel.biotools.patterns
   :members:


Biotools
--------

.. automodule:: dnachisel.biotools.biotools
  :members:
