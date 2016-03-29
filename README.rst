# DNA Chisel

In DNA Chisel you define constraints on a DNA sequence, such as:
- DNA coding sequences that must be preserved (only synonymous mutations are authorized)
- Restriction sites that must be preserved
-
- Remove patterns from a sequence, e.g. restriction sites
- Make patterns appear in a DNA sequence, e.g. a restriction site.
- Codon-optimize a sequence for some organism.



## Examples

canvas = DNACanvas(
    sequence = "ATGCGTGTGCGTATGCGTGTGTGCGTGATG",
    constraints = [
        Pattern("ATTCTT", window = [100, 200]),
        NoPattern("AGTC", window = [300, 600]),

        PreserveORF(start1, end1),
        PreserveORF(start2, end2),
        PreserveORF(start3, end3),

        GCPercent(min=0.4, max=0.6, window=100)
        GCPercent(min=0.7, max=0.6, window=100)
    ]
    objective = [
        CodonOptimization(start=, end=, organism=)

    ]
)
max_local_gc_content(percentage=20, window=60)
optimization= []

# Codon optimization


## Installation
