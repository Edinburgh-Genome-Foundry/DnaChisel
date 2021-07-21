# Code explanation

This document walks you trough the DNA Chisel code. Please request changes if anything is unclear.

## Core classes

Here we explain the relationship between the main classes in DNA Chisel. The next section will explain how the code is organized.

- **DnaOptimizationProblem** is the class to define and solve an optimization problems. Its methods implement all the solver logics. Notable attributes:

  - Constraints (list of _Specifications_ that must be verified in the final sequence)
  - Objectives (list of _Specifications_ that must be as close to verified as possible in the final sequence)
  - A _MutationSpace_ storing the changes that can be made to the sequence

- **CircularDnaOptimizationProblem** is a variant of _DnaOptimizationProblem_ whose optimization algorithm assumes that the sequence is circular.
- **Specification** is a class to represent sequence design rules that can be used either as constraints or objectives. Its methods define what the specification looks at in a sequence, how it behaves on local problems, etc. The evaluation of a _DnaOptimizationProblem_ by a *Specification* produces a *SpecEvaluation*. Notable attributes of *Specifications*:

  - A _Location_ indicating where it applies in a sequence (absence of location implies that it applies to the whole sequence).
  - A _SequencePattern_ if the Specification is about inserting or removing a pattern.

- **SpecificationSet** (same file) is a special Specification class which actually groups together several _Specifications_. See [AllowPrimers](https://github.com/Edinburgh-Genome-Foundry/DnaChisel/blob/master/dnachisel/builtin_specifications/AllowPrimer.py) for an example.
- **Location** is a class representing the location of a segment of a sequence. It is used primarily to indicate where a _Specification_ applies in the sequence, but also in _SpecEvaluation_ to indicate constraint breaches, in _SequencePattern_ to indicate the location of a pattern, etc.
  Class methods describe how _Locations_ are created, extended, merged, intersected, imported from Genbank locations, etc.
- **SpecEvaluation** is a class describing the result of the evaluation of a *DnaOptimizationProblem* by a *Specification*. Contains a score, a message, a list of *Locations* of sub-optimal regions. Several evaluations can be grouped using classes *ProblemConstraintsEvaluations* and *ProblemObjectivesEvaluations*, which implement methods for printing or exporting as Genbank a set of evaluations.
- **SequencePattern** is a class for representing pattern. Methods define how to parse a pattern (such as "BsmBI_site") and find the pattern in a sequence.
- **MutationSpace** is a class to represent the possible mutations at different locations of the sequence for a given problem. It is initialized at the creation of the problem by the problem's constraints. A MutationSpace is basically a list of *MutationChoices* applying at various locations of the sequence. *MutationSpace* methods describe how to extract variants from the mutation space (this is used by DnaOptimizationProblem's solver to explore new sequence variants). A *MutationChoice* comprises a location and a set of sequence choices for this location, and the methods define how to extract an random choice, how to merge *MutationChoices* together, etc.

## Code organization

- **biotools/** contains many methods and data tables related to biology and sequence manipulation, either used in the core DNA Chisel classes, or very helpful when writing DNA Chisel scripts (see this folder's README for more).  
- **builtin specifications/** contains all the built-in _Specification_ subclasses, one per file. While most files define a directly usable *Specification* subclass (CodonOptimize, EnforceGCContent, etc.), some files encode generic subclasses meant to be in turn subclassed (*CodonSpecification*, *TerminalSpecification*)
- **DnaOptimizationProblem/** contains the code for the *DnaOptimizationProblem* and *CircularDnaOptimizationProblem* classes. As *DnaOptimizationProblem* implements the solver and is very big, the methods in this class have been regrouped into "mixins" in other files.
- **MutationSpace/** contains the implementation of the *MutationSpace* and *MutationChoice* classes.
- **reports/** contains methods to generate plots, PDF reports, etc. from a _DnaOptimizationProblem_ (before and after its optimization). It also contains assets (logo, stylesheet, template) for the PDF report.
- **Specification/** contains the implementation of the *Specification* class (some methods for this class are regrouped in file *FeatureRepresentationMixin.py* to make smaller files). The *SpecEvaluation/* subfolder contains the implementations of *SpecEvaluation*, *SpecEvaluations*, *ProblemConstraintEvaluations*, *ProblemObjectiveEvaluations*.
- **utils/** contains generic methods built on top of DNA Chisel that can be useful in other libraries/applications. For instance `random_compatible_dna_sequence()` which returns a random sequence verifying the given constraints.
- **Location_py** implements the *Location* class which is used everywhere to define sequence segments.
- **SequencePattern_py**  implements the *SequencePattern* class.
