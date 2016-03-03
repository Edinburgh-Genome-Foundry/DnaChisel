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
        Pattern(start, end, pattern, number=1),
        NoPattern(start, end, pattern, number=0),

        PreserveCDS(start, end),
        PreserveCDS(start, end),
        PreserveCDS(start, end),

        LocalGCContent()
        GlobalGCContent()


    ]
    optimizations = [
        codon_optimization(start=, end=, organism=)
    ]

)
max_local_gc_content(percentage=20, window=60)
optimization= []

# Codon optimization


## Installation
