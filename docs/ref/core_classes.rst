Core Classes
-------------

Classes dependencies
~~~~~~~~~~~~~~~~~~~~~

.. mermaid::
   :align: center

   graph TD;
      sequence --> P[DnaOptimizationProblem]
      o[objectives...] --> P
      c[constraints...] --> P
      s[Specifications] --> o
      s --> o
      s --> c
      s --> c
      P -->|solve, optimize| O[Optimized problem.sequence]

DnaOptimizationProblem
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.DnaOptimizationProblem
   :members:

Specification
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.Specification
   :members:

SpecEvaluation
~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: dnachisel.Location
   :members:

MutationSpace
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: dnachisel.MutationSpace
   :members:

SequencePattern
~~~~~~~~~~~~~~~~~~~~~~~

Specifications ``AvoidPattern``, ``EnforcePatternOccurence`` accept a Pattern
as argument. The pattern can be provided either as a string or as a class:

================== ================================= ====================
     Name                Pattern Class                    As string
================== ================================= ====================
Restriction site    ``EnzymeSitePattern('BsmBI')``    ``'BsmBI_site'``
Homopolymer         ``HomopolymerPattern("A", 6)``    ``'6xA'``
Direct repeats      ``RepeatedKmerPattern(3, 2)``     ``'3x2mer'``
DNA notation        ``DnaNotationPattern('ANNKA')``   ``'ANNKA'``
Regular Expression  ``SequencePattern('A[CG]*A')``    ``'A[CG]*A'``
================== ================================= ====================

.. automodule:: dnachisel.SequencePattern
   :members:
