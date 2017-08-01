"""Module with collections of pre-sets objectives."""

from .builtin_specifications import (AvoidPattern, EnforceGCContent,
                                     EnforceTerminalGCContent,
                                     AvoidHairpins)
from ..biotools import (homopolymer_pattern, repeated_kmers)

PROVIDERS_CONSTRAINTS = {
    "Gen9": [
        AvoidPattern(enzyme="BsaI"),
        AvoidPattern(enzyme="AarI"),
        AvoidPattern(homopolymer_pattern("A", 9)),
        AvoidPattern(homopolymer_pattern("T", 9)),
        AvoidPattern(homopolymer_pattern("G", 6)),
        AvoidPattern(homopolymer_pattern("C", 9)),
        EnforceGCContent(0.4, 0.65),
        EnforceGCContent(0.25, 0.80, gc_window=50)
    ],
    "IDT": [
        EnforceGCContent(0.25, 0.68),
        EnforceGCContent(0.28, 0.76, gc_window=100),
        EnforceGCContent(0.15, 0.90, gc_window=20),
        EnforceTerminalGCContent(window_size=30, gc_min=0.24, gc_max=0.76),
        AvoidPattern(homopolymer_pattern("A", 13)),
        AvoidPattern(homopolymer_pattern("T", 13)),
        AvoidPattern(homopolymer_pattern("G", 6)),
        AvoidPattern(homopolymer_pattern("C", 6)),
        AvoidPattern(repeated_kmers(3, n_repeats=5)),
        AvoidPattern(repeated_kmers(2, n_repeats=9)),
        AvoidHairpins(stem_size=20, hairpin_window=200)
    ],
    "Twist": [
        EnforceGCContent(0.4, 0.65),
        AvoidPattern(homopolymer_pattern("A", 9)),
        AvoidPattern(homopolymer_pattern("T", 9)),
        AvoidPattern(homopolymer_pattern("G", 9)),
        AvoidPattern(homopolymer_pattern("C", 9))
    ]
}
