# Changes

## v3.0.0

**API-breaking changes:**

- The strand on which to look for a pattern in AvoidPattern and
  EnforcePatternOccurence is not given by the location (strand 1, -1, or strand
  0 to look in both strands).
- AvoidNonUniqueSegments renamed in UniquifyKmers
- Complete rewrite of CodonOptimization, method names are now use_best_codon,


**New features**

- Support for different genetic codes in EnforceTranslation.
- Support for non-standard start codons in EnforceTranslation
- New AvoidMatches specification based on Bowtie.
- AvoidChanges and EnforceChanges can now tunable (EnforceChanges(60%),
  AvoidChanges(max_edits=2))
- New SequencePattern MotifPssmPattern.
  
**Misc.**

- Got rid of some tables like Genetic code, now taken from Biopython.
- Complete code restructuring by cutting large files into smaller ones.

## v2.0.0

- ``problem.n_mutations`` became ``problem.mutations_per_iteration`` for clarity
- ``'harmonized'`` option in ``CodonOptimize()`` replaced by ``harmonized_frequencies`` for clarity