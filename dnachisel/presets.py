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
    ],
    "IDT": [
        cst.GCContentConstraint(0.25, 0.68),
        cst.GCContentConstraint(0.28, 0.76, gc_window=100),
        cst.GCContentConstraint(0.15, 0.90, gc_window=20),
        cst.TerminalGCContentConstraint(0.24, 0.76, window_size=30),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("A", 13)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("T", 13)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("G", 6)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("C", 6)),
        cst.NoPatternConstraint(patterns.repeated_kmers(3, n_repeats=5)),
        cst.NoPatternConstraint(patterns.repeated_kmers(2, n_repeats=9)),
        cst.NoHairpinsIDTConstraint(stem_size=20, hairpin_window=200)
    ],
    "Twist": [
        cst.GCContentConstraint(0.4, 0.65),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("A", 9)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("T", 9)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("G", 9)),
        cst.NoPatternConstraint(patterns.homopolymer_pattern("C", 9))
    ]
}
