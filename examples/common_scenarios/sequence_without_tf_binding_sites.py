from dnachisel import DnaOptimizationProblem, AvoidPattern, random_dna_sequence
from urllib import request
import pandas as pd
from io import StringIO

# DOWNLOAD THE LIST OF TF BINDING SITES
url = "http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt"
data = request.urlopen(url).read().decode("utf-8")
df = pd.read_csv(
    StringIO(data), sep="\t", skiprows=(range(0, 46)), header=None
)  # First 46 lines are description

# OBTAIN A LIST OF TF BINDING SEQUENCES
tf_column = df[13]  # 14th column contains TF binding sites
tf_column.dropna(inplace=True)
tf_list = tf_column.to_list()
# According to the description, the binding sites are in uppercase, so we remove lowercase:
tf_binding_sequences = ["".join(ch for ch in tf if not ch.islower()) for tf in tf_list]
# Remove single-nucleotide TFs:
tf_binding_sequences = [tf for tf in tf_binding_sequences if len(tf) > 1]

# DEFINE AND SOLVE THE OPTIMIZATION PROBLEM
problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(50000),
    constraints=[AvoidPattern(pattern) for pattern in tf_binding_sequences],
)
problem.resolve_constraints()
problem.to_record("sequence_without_tf_binding_sites.gb")
