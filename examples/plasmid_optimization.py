import urllib
from dnachisel import *
from Bio import SeqIO

# DOWNLOAD THE PLASMID

url ="http://www.stevekellylab.com/constructs/pDex/pDex577-G.gb"
urllib.urlretrieve(url, "pDex577-G.gb")


# PARSE THE PLASMID FILE WITH BIOPYTHON, GET SEQUENCE AND GENES LOCATIONS

annotated_sequence = SeqIO.read(open("pDex577-G.gb"), "genbank")
sequence = str(annotated_sequence.seq)
CDS_list = [
    (int(f.location.start), int(f.location.end), int(f.location.strand))
    for f in annotated_sequence.features
    if f.type == "CDS"
]


# DEFINE THE CONSTRAINTS

GEN9_constraints = [
    NoPatternConstraint(enzyme_pattern("BsaI")),
    NoPatternConstraint(enzyme_pattern("AarI")),
    NoPatternConstraint(homopolymer_pattern("A",9)),
    NoPatternConstraint(homopolymer_pattern("T",9)),
    NoPatternConstraint(homopolymer_pattern("G",6)),
    NoPatternConstraint(homopolymer_pattern("C",9)),
    GCPercentConstraint(0.4, 0.65),
    GCPercentConstraint(0.25, 0.80, gc_window=50)
]
CDS_constraints = [
    EnforceTranslationConstraint(
        (start, end),
        sequence=sequence,
        strand=strand
    )
    for (start, end, strand) in CDS_list
]


# DEFINE OBJECTIVES

objectives =

# DEFINE AND SOLVE THE PROBLEM
canvas = DNACanvas(
    sequence=sequence,
    constraints= GEN9_constraints + CDS_constraints,
    objectives = objectives
)

print ("\n\n=== Status before optimization ===")
canvas.print_constraints_summary()
canvas.solve_all_constraints_one_by_one()
print ("\n\n=== Status after optimization ===")
canvas.print_constraints_summary()
