import constraints as cst
import objectives as obj
import patterns

PROVIDERS_CONSTRAINTS = {
    "Gen9": [
        cst.NoPatternConstraint(patterns.enzyme_pattern("BsaI")),
        cst.NoPatternConstraint(patterns.enzyme_pattern("AarI")),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("A", 9)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("T", 9)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("G", 6)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("C", 9)),
        cst.GCContentConstraint(0.4, 0.65),
        cst.GCContentConstraint(0.25, 0.80, gc_window=50)
    ]
}
