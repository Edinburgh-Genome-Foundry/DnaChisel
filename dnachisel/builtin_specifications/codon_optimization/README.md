# Codon Optimization specifications

This submodule contains various classes implementing different flavours of codon
optimization that one can find in the literature.

- MaximizeCAI will replace every codon by the "best synonymous codon", i.e. the
  most frequent in the target species.
- MatchTargetCodonUsage will optimize the sequence so it has the same codon
  usage (=frequencies profile) as the target species.
- Harmonize RCA will make so that codons that were rare in the original sequence's
  "original species" will also be rare in the final sequence for the target species.
- CodonOptimize is a generic pseudo-class which uses a "mode" parameter to return
  a specification of one of the above classes.
