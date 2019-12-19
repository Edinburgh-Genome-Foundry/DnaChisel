from urllib import request
from io import StringIO
from dnachisel import (
    DnaOptimizationProblem,
    load_record,
    Location,
    EnforceTranslation,
    EnforceGCContent,
    AvoidPattern,
)

print("DOWNLOADING AND PARSING THE GENBANK DATA...")

url = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    + "db=nucleotide&id=48994873&rettype=gb&retmode=txt"
)
genbank_data = request.urlopen(url).read().decode("utf-8")
genbank_record = load_record(StringIO(genbank_data), file_format="genbank")

print("INITIALIZING THE PROBLEM WITH CONSTRAINTS FOR EACH GENE...")

constraints = []
for feature in genbank_record.features:
    if feature.type == "gene" and len(feature.location.parts) == 1:
        location = Location.from_biopython_location(feature.location)
        if (len(location) % 3 == 0) and len(location) > 100:
            gene_constraints = [
                EnforceTranslation(location=location),
                AvoidPattern("BsmBI_site", location, strand='both'),
                EnforceGCContent(
                    mini=0.40, maxi=0.60, window=150, location=location
                ),
            ]
            constraints.extend(gene_constraints)
problem = DnaOptimizationProblem(genbank_record, constraints)

print("RESOLVING THE CONSTRAINTS...")

problem.logger.ignore_bars_under = 50
problem.resolve_constraints()
problem.to_record("ecoli_genes_optimization.gb")
