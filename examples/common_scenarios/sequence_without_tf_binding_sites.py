from dnachisel import DnaOptimizationProblem, AvoidPattern, random_dna_sequence
from urllib import request

# DOWNLOAD THE LIST OF TF BINDING SITES
url = "http://regulondb.ccg.unam.mx/menu/download/datasets/files/PSSMSet.txt"
data = request.urlopen(url).read().decode("utf-8")

# PARSE THE DATA LINE BY LINE TO OBTAIN A LIST OF TF BINDING SEQUENCES
tf_binding_sequences = [
    line for line in data.splitlines() if set() < set(line) <= set("ATGC")
]

# DEFINE AND SOLVE THE OPTIMIZATION PROBLEM
problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(50000),
    constraints=[AvoidPattern(pattern) for pattern in tf_binding_sequences],
)
problem.resolve_constraints()
problem.to_record("sequence_without_tf_binding_sites.gb")
