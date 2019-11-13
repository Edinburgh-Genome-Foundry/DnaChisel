from Bio.Blast import NCBIXML
import os
import tempfile
import subprocess
import time


def blast_sequence(
    sequence,
    blast_db=None,
    subject_sequences=None,
    subject=None,
    word_size=4,
    perc_identity=80,
    num_alignments=1000,
    ungapped=False,
    num_threads=3,
    culling_limit=None,
    e_value=None,
    task=None,
    dust="no",
):
    """Return a Biopython BLAST record of the given sequence BLASTed
    against the provided database.

    Parameters
    ----------

    sequence
      An ATGC sequence

    Examples
    --------

    >>> blast_record = blast_sequence("ATTGTGCGTGTGTGCGT", "blastdb/ecoli")
    >>> for alignment in blast_record.alignments:
    >>>     for hit in alignment.hsps:
    >>>         print (hit.identities)
    """

    xml_file, xml_name = tempfile.mkstemp(".xml")
    fasta_file, fasta_name = tempfile.mkstemp(".fa")
    with open(fasta_name, "w+") as f:
        f.write(">seq\n" + sequence)

    close_subject = False
    remove_subject = False

    if subject is not None:
        close_subject = True

    if subject_sequences is not None:
        close_subject = True
        remove_subject = True
        subject_file, subject = tempfile.mkstemp(".fa")
        if isinstance(subject_sequences[0], str):
            subject_sequences = [
                ("%06d" % i, seq) for i, seq in enumerate(subject_sequences)
            ]
        fasta_content = "\n".join(
            [">%s\n%s" % name_sequence for name_sequence in subject_sequences]
        )
        with open(subject, "w+") as f:
            f.write(fasta_content)

    def parameter_if_not_none(label, param):
        return [label, str(param)] if param else []

    command = [
        "blastn",
        "-out",
        xml_name,
        "-outfmt",
        "5",
        "-max_target_seqs",
        str(num_alignments),
        "-query",
        fasta_name,
        "-word_size",
        str(word_size),
        "-num_threads",
        str(num_threads),
        "-perc_identity",
        str(perc_identity),
    ]
    command += (
        (["-db", blast_db] if subject is None else ["-subject", subject])
        + parameter_if_not_none("-dust", dust)
        + parameter_if_not_none("-evalue", e_value)
        + parameter_if_not_none("-culling_limit", culling_limit)
        + parameter_if_not_none("-task", task)
    )
    if ungapped:
        command += ["-ungapped"]
    
    process = subprocess.run(command, stderr=subprocess.PIPE, close_fds=True)
    if process.returncode:
        raise OSError("BLAST failed: %s" % process.stderr)
    with open(xml_name, "r") as f:
        result = list(NCBIXML.parse(f))
    os.fdopen(xml_file, "w").close()
    os.fdopen(fasta_file, "w").close()
    os.remove(xml_name)
    os.remove(fasta_name)
    if close_subject:
        open(subject, "w").close()
        if remove_subject:
            os.remove(subject)
    if len(result) == 1:
        return result[0]
    else:
        return result
