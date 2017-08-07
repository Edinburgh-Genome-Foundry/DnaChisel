"""Example of use for AvoidBLASTMatches

In this example we create a 1000bp random sequence, then edit out every match
with E. coli that is 14bp or longer.

"""
import os
import gzip
import shutil
from urllib.request import urlretrieve
import subprocess

from dnachisel import (DnaOptimizationProblem, random_dna_sequence,
                       AvoidBlastMatches)

# THIS CREATES THE BLAST DATABASE IF NOT ALREADY HERE

ecoli_genome_url = (
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845"
    "/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
)
genome_gz = os.path.join('downloaded_data', 'ecoli_genome.gz')
genome_fasta = os.path.join('downloaded_data', 'ecoli_genome.fa')
genome_blastdb = os.path.join('downloaded_data', 'ecoli_genome')

if not os.path.exists('downloaded_data'):
    os.mkdir('downloaded_data')
if not os.path.exists(os.path.join('downloaded_data', 'ecoli_genome.nsq')):

    urlretrieve(ecoli_genome_url, genome_gz)

    with open(genome_fasta, 'wb') as f_fasta:
        with gzip.open(genome_gz, 'rb') as f_gz:
            shutil.copyfileobj(f_gz, f_fasta)

    subprocess.Popen(["makeblastdb",
                      "-in", genome_fasta,
                      "-dbtype", "nucl",
                      "-out", genome_blastdb])


# DEFINE AND SOLVE THE PROBLEM

problem = DnaOptimizationProblem(
    sequence=random_dna_sequence(1000, seed=123),
    constraints=[
        AvoidBlastMatches(blast_db=genome_blastdb, min_align_length=13,
                          perc_identity=80)
    ]
)

print ("\nBefore optimization\n")
print (problem.constraints_text_summary())
problem.resolve_constraints(progress_bars=2, final_check=False)
print ("\nAfter optimization\n")
print (problem.constraints_text_summary())
