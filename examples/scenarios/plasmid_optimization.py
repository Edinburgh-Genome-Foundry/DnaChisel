"""Example of plasmid optimization with DnaChisel.

In this example, we download a plasmid from the web in GENBANK format and
modify the sequence with respect to the following constraints and objectives:

- For each coding sequence (CDS) found in the Genbank file:
    - Make sure the corresponding protein sequence remains unchanged.
    - Do not modify the 30-base-pair region upstream of the CDS (promoter)
    - Codon-optimize for E. coli (use the most frequent codons preferentially)

- Enforce the constraints given by the DNA provider Gen9:
    - Make sure there is no restriction site for BsaI and AarI in the sequence
    - Make sure there is no 9-homopolymers of A,T,C or any 6-homopolymers of G.
    - Make sure the GC content remains between 40 and 65 percent globally and
      between 25 and 80 percent over sliding 100-base-pair windows

- Aim at a global GC content of 51 percent like in E. coli

At the beginning, many of the constraints fail (there are BsaI sites, A- and T-
homopolymers, and the GC-content is out-of-bounds in many locations).
After a run of the constraint solver (~ 1 second) all constraints are solved.
And after a run of the optimization algorithm (~20 seconds) the objectives are
much improved.

The final sequence (with the original annotations) is exported to Genbank.
"""

from dnachisel import (
    DnaOptimizationProblem,
    AvoidPattern,
    AvoidChanges,
    EnforceTranslation,
    HomopolymerPattern,
    EnforceGCContent,
    CodonOptimize,
    load_record,
)
from io import StringIO
import urllib


# DOWNLOAD THE PLASMID FROM THE WEB (it is a 7kb plasmid with 3 genes)
url = "http://www.stevekellylab.com/constructs/pDex/pDex577-G.gb"
response = urllib.request.urlopen(url)
record_file = StringIO(response.read().decode("utf-8"))
record = load_record(record_file, fmt="genbank")


CDS_list = [
    (int(f.location.start), int(f.location.end), int(f.location.strand))
    for f in record.features
    if f.type == "CDS"
]


# DEFINE CONSTRAINTS

dna_provider_constraints = [
    AvoidPattern("BsaI_site"),
    AvoidPattern("AarI_site"),
    AvoidPattern("9xA"),
    AvoidPattern("9xT"),
    AvoidPattern(HomopolymerPattern("6xG")),
    AvoidPattern(HomopolymerPattern("6xC")),
    EnforceGCContent(0.4, 0.65),
    EnforceGCContent(0.25, 0.80, window=50),
]

CDS_constraints = []
for (start, end, strand) in CDS_list:
    if strand == 1:
        promoter_region = (start - 30, start - 1)
    else:
        promoter_region = (end + 1, end + 30)
    CDS_constraints += [
        AvoidChanges(promoter_region),
        EnforceTranslation((start, end, strand)),
    ]


# DEFINE OBJECTIVES

objectives = [EnforceGCContent(0.51, boost=10000)] + [
    CodonOptimize("e_coli", location=(start, end, strand))
    for (start, end, strand) in CDS_list
]


# DEFINE AND SOLVE THE PROBLEM

problem = DnaOptimizationProblem(
    sequence=record,
    constraints=dna_provider_constraints + CDS_constraints,
    objectives=objectives,
)

print("\n\n=== Initial Status ===")
print(problem.constraints_text_summary(failed_only=True))
print(problem.objectives_text_summary())

print("Now solving constraints...")
problem.resolve_constraints()
print(problem.constraints_text_summary(failed_only=True))

print("Now optimizing objectives...")

problem.optimize()

print("\n\n=== Status after optimization ===\n\n")
print(problem.objectives_text_summary())
print("Let us check again on the constraints:")
print(problem.constraints_text_summary(failed_only=True))
