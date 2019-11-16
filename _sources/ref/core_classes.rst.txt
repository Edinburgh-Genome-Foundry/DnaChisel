Core Classes
=============

.. contents::

DnaOptimizationProblem
-----------------------

.. autoclass:: dnachisel.DnaOptimizationProblem
   :members:

CircularDnaOptimizationProblem
------------------------------

.. autoclass:: dnachisel.CircularDnaOptimizationProblem
   :members:

Specification
-----------------------

.. autoclass:: dnachisel.Specification
   :members:

SpecEvaluation
-----------------------

.. .. mermaid::

..    graph TD;
..       P[DnaOptimizationProblem] --> SE[SpecEvaluation]
..       S[Specification] --> SE
..       SE -->|grouped into| SES[SpecEvaluations]
..       SES --- OES[Objective evaluations]
..       SES --- CES[Constraints evaluations]

.. automodule:: dnachisel.SpecEvaluation
   :members:

Location
-----------------------

.. automodule:: dnachisel.Location
   :members:

MutationSpace
-----------------------

.. automodule:: dnachisel.MutationSpace
   :members:

SequencePattern
-----------------------

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

SequencePattern class
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.SequencePattern
   :members:


DnaNotationPattern
~~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.DnaNotationPattern
   :members:

MotifPssmPattern
~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.MotifPssmPattern
   :members:

EnzymeSitePattern
~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.EnzymeSitePattern
   :members:

RepeatedKmerPattern
~~~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.RepeatedKmerPattern
   :members:

HomopolymerPattern
~~~~~~~~~~~~~~~~~~~

.. autoclass:: dnachisel.HomopolymerPattern
   :members: