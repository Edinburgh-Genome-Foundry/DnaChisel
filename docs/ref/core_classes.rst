

DnaOptimizationProblem
----------------------

.. automodule:: dnachisel.DnaOptimizationProblem
   :members:

Specification
-------------

.. automodule:: dnachisel.Specification
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

.. automodule:: dnachisel.SpecEvaluation
  :members:

Location
---------

.. automodule:: dnachisel.Location
   :members:

MutationSpace
--------------

.. automodule:: dnachisel.MutationSpace
   :members:

SequencePattern
----------------

.. automodule:: dnachisel.SequencePattern
  :members:
