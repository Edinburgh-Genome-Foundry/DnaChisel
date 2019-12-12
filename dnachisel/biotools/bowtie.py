import subprocess
import tempfile
import os

def run_process(name, parameters):
    process = subprocess.run(
        parameters,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    if process.returncode:
        error = process.stderr.decode()
        raise OSError("%s failed:\n\n%s" % (name, error))
    return process.stdout

def create_bowtie_index_from_sequences(sequences, path):
    fasta_path = os.path.join(path, 'sequences.fa')
    bowtie_path = os.path.join(path, 'bowtie')
    with open(fasta_path, 'w') as f:
        f.write("\n".join([
            ">%d\n%s" % (i, sequence)
            for i, sequence in enumerate(sequences)
        ]))
    run_process("build-bowtie", [
        "bowtie-build", "-f", fasta_path, bowtie_path, "--quiet"
    ])
    return bowtie_path

def find_all_bowtie_matches(
    sequence, bowtie_index_path, match_length, max_mismatches=0
):
    """Return (short) matches between a sequence and a Bowtie index.
    
    The result is of the form [(start, end), n_mismatches)] where (start, end)
    indicates the position of the match in the sequence, and n_mismatches is
    the number of mismatches with the closest homology in the index.
    
    
    """

    # CREATE THE PARAMETERS

    parameters = ["bowtie"]
    parameters += ["--best", "-k", "1"]  # only return the best alignments
    parameters += ["-v", str(max_mismatches)]  # only allow that N mismatches
    parameters += [bowtie_index_path]
    parameters += ["--quiet", "--suppress", "2,3,4,5,6,7"]  # small output
    k = match_length
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    if k * len(kmers) < 10000:
        # Input the sequences directly
        tmp_fasta_path = None
        parameters += ["-c", ",".join(kmers)]
    else:
        # Write sequences to a file if too many.
        tmp_fasta_path = tempfile.mktemp(".fa")
        with open(tmp_fasta_path, "w") as f:
            entries = [">%d\n%s" % (i, s) for i, s in enumerate(kmers)]
            f.write("\n\n".join(entries))

        parameters += ["-f", tmp_fasta_path]

    # RUN THE PROCESS
    try:
        output = run_process("BOWTIE", parameters)
    except Exception as err:
        raise err
    finally:
        if tmp_fasta_path is not None:
            os.remove(tmp_fasta_path)
    output_records = [
        line.split("\t") for line in output.decode().split("\n") if len(line)
    ]
    return [
        ((int(index), int(index) + k), edits.count(":"))
        for index, edits in output_records
    ]