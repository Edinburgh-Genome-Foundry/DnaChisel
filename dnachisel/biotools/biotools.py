import tempfile
import os
import subprocess
import time
import itertools
from copy import deepcopy

import numpy as np
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
from Bio import Restriction

from .biotables import (CODONS_TRANSLATIONS, NUCLEOTIDE_TO_REGEXPR,
                        COMPLEMENTS, CODONS_SEQUENCES, IUPAC_NOTATION)


def complement(dna_sequence):
    """Return the complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"TACGGC"``.

    Uses BioPython for speed.
    """
    if len(dna_sequence) <= 30:
        return "".join([COMPLEMENTS[nuc] for nuc in dna_sequence])
    # This alternative has overhead but is really fast on long sequences
    return str(Seq(dna_sequence).complement())


def reverse_complement(sequence):
    """Return the reverse-complement of the DNA sequence.

    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.

    Uses BioPython for speed.
    """
    return complement(sequence)[::-1]


def is_palyndromic(dna_sequence):
    """Return True if the DNA sequence is equal to its reverse complement."""
    return reverse_complement(dna_sequence) == dna_sequence


def random_dna_sequence(length, gc_share=None, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    proba
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)
    if gc_share is not None:
        g_or_c = gc_share / 2.0
        not_g_or_c = (1 - gc_share) / 2.0
        probas = {"G": g_or_c, "C": g_or_c, "A": not_g_or_c, "T": not_g_or_c}
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def load_record(filename, linear=True, name="unnamed", fmt='auto'):
    """Load a FASTA/Genbank/... record"""
    if fmt is not 'auto':
        record = SeqIO.read(filename, fmt)
    elif filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(('fa', 'fasta')):
        record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError('Unknown format for file: %s' % filename)
    record.linear = linear
    if name != "unnamed":
        record.id = name
        record.name = name.replace(" ", "_")[:20]
    return record

def random_protein_sequence(length, seed=None):
    """Return a random protein sequence "MNQTW...YL*" of the specified length.

    Parameters
    ----------

    length
      Length of the protein sequence (in number of amino-acids). Note that the
      sequence will always start with ``"M"`` and end with a stop codon ``"*"``
      with (length-2) random amino-acids in the middle

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility

    """
    if seed is not None:
        np.random.seed(seed)

    aa_list = list('ACEDGFIHKLNQPSRTWVY')
    aa_choices = np.random.choice(aa_list, length - 2)
    return "M" + "".join(aa_choices) + "*"


def reverse_translate(protein_sequence):
    """Return a DNA sequence which translates to the provided protein sequence

    Note: at the moment, the first valid codon found is used for each
    amino-acid (so it is deterministic but no codon-optimization is done).
    """
    return "".join([
        CODONS_SEQUENCES[aa][0]
        for aa in protein_sequence
    ])


def translate(dna_sequence, translation_table="Bacterial"):
    """Translate the DNA sequence into an amino-acids sequence "MLKYQT...".
    If ``translation_table`` is the name or number of  NCBI genetic table,
    Biopython will be used. See here for options:

    http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25

    ``translation_table`` can also be a dictionnary of the form
    ``{"ATT": "M", "CTC": "X", etc.}`` for more exotic translation tables


    """
    if isinstance(translation_table, dict):
        return "".join([
            translation_table[dna_sequence[i:i + 3]]
            for i in range(0, len(dna_sequence), 3)
        ])
    else:
        return str(Seq(dna_sequence).translate(table=translation_table))


def dna_pattern_to_regexpr(dna_pattern):
    """Return a regular expression pattern for the provided DNA pattern

    For instance ``dna_pattern_to_regexpr('ATTNN')`` returns
    ``"ATT[A|T|G|C][A|T|G|C]"``.
    """
    return "".join([
        NUCLEOTIDE_TO_REGEXPR[n]
        for n in dna_pattern
    ])


def windows_overlap(window1, window2):
    """Return the overlap span between two windows.

    Parameters
    ----------

    window1, window2
      Each window is a couple of the form (start, end) indicating the range of
      a segment of integers.

    Returns
    -------

    None
      In case the two windows do not overlap.
    [start, end]
      The coordinates of the overlap segment if there is one.
    """
    start1, end1 = window1
    start2, end2 = window2

    if start2 < start1:
        return windows_overlap(window2, window1)

    if start1 <= start2 <= end1:
        return [start2, min(end1, end2)]
    else:
        return None


def gc_content(sequence, window_size=None):
    """Compute global or local GC content.

    Parameters
    ----------

    sequence
      An ATGC DNA sequence (upper case!)

    window_size
      If provided, the local GC content for the different sliding windows of
      this size is returned, else the global GC content is returned.

    Returns
    --------

      A number between 0 and 1 indication the proportion
      of GC content. If window_size is provided, returns
      a list of len(sequence)-window_size values indicating
      the local GC contents (sliding-window method). The i-th value
      indicates the GC content in the window [i, i+window_size]
    """
    # The code is a little cryptic but speed gain is 300x
    # compared with pure-python string operations

    arr = np.fromstring(sequence + "", dtype="uint8")
    arr_GCs = (arr == 71) | (arr == 67)  # 67=C, 71=G

    if window_size is None:
        return 1.0 * arr_GCs.sum() / len(sequence)
    else:
        cs = np.cumsum(arr_GCs)
        a = cs[window_size - 1:]
        b = np.hstack([[0], cs[:-window_size]])
        return 1.0 * (a - b) / window_size

def blast_sequence(sequence, blast_db=None, subject_sequences=None,
                   subject=None, word_size=4, perc_identity=80,
                   num_alignments=1000, ungapped=False, num_threads=3,
                   culling_limit=None, e_value=None, use_megablast=True,
                   dust="no"):
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
                ("%06d" % i, seq)
                for i, seq in enumerate(subject_sequences)
            ]
        fasta_content = "\n".join([
            ">%s\n%s" % name_sequence
             for name_sequence in subject_sequences
        ])
        with open(subject, "w+") as f:
            f.write(fasta_content)

    def parameter_if_not_none(label, param):
        return [label, str(param)] if param else []

    command = [
        "blastn", "-out", xml_name,
        "-outfmt", "5",
        "-max_target_seqs", str(num_alignments),
        "-query", fasta_name,
        "-word_size", str(word_size),
        "-num_threads", str(num_threads),
        "-perc_identity", str(perc_identity)
    ]
    command += (
        (["-db", blast_db] if subject is None else ['-subject', subject]) +
        parameter_if_not_none("-dust", dust) +
        parameter_if_not_none("-evalue", e_value) +
        parameter_if_not_none("-culling_limit", culling_limit)
    )
    if use_megablast:
        command += ["-task", "megablast"]
    if ungapped:
        command += ["-ungapped"]

    p = subprocess.Popen(command, close_fds=True)
    out, err = p.communicate()
    p.wait()
    n_trials = 3
    for i in range(n_trials):
        try:
            with open(xml_name, "r") as f:
                res = list(NCBIXML.parse(f))
                os.fdopen(xml_file, 'w').close()
                os.fdopen(fasta_file, 'w').close()
                os.remove(xml_name)
                os.remove(fasta_name)
                if close_subject:
                    open(subject, 'w').close()
                    if remove_subject:
                        os.remove(subject)
                if len(res) == 1:
                    return res[0]
                else:
                    return res
            break
        except ValueError as err:
            if i == n_trials - 1:
                raise err
            time.sleep(0.1)

def subdivide_window(window, max_span):
    """Subdivide a window (start, end) into windows of size < max_span
    (start, i_1), (i_1, i_2), ... (i_n, end)"""
    start, end = window
    inds = list(range(start, end, max_span)) + [end]
    return list(zip(inds, inds[1:]))


def change_biopython_record_sequence(record, new_seq):
    """Return a version of the record with the sequence set to new_seq"""
    new_record = deepcopy(record)
    new_record.seq = Seq(new_seq, alphabet=DNAAlphabet())
    return new_record


def sequence_to_biopython_record(sequence, id='<unknown id>',
                                 name='<unknown name>', features=()):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(Seq(sequence, alphabet=DNAAlphabet()),
                     id=id, name=name, features=list(features))


def find_specification_in_feature(feature):
    """Analyse a Biopython feature to find a DnaChisel Specification in it.

    The specification should start with either "@" or "~", in the feature's
    field "label" or "note".
    """
    for labelfield in ["label", "note"]:
        if labelfield not in feature.qualifiers:
            continue
        potential_label = feature.qualifiers.get(labelfield, "_")
        if isinstance(potential_label, list):
            potential_label = potential_label[0]
        if (potential_label != '') and (potential_label[0] in "@~"):
            return potential_label
    return None

def sequences_differences_array(seq1, seq2):
    """Return an array [0, 0, 1, 0, ...] with 1s for sequence differences.

    seq1, seq2 should both be ATGC strings.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Only use on same-size sequences (%d, %d)" %
                         (len(seq1), len(seq2)))
    arr1 = np.fromstring(seq1, dtype="uint8")
    arr2 = np.fromstring(seq2, dtype="uint8")
    return arr1 != arr2


def sequences_differences(seq1, seq2):
    """Return the number of nucleotides that differ in the two sequences.

    seq1, seq2 should be strings of DNA sequences e.g. "ATGCTGTGC"
    """
    return sequences_differences_array(seq1, seq2).sum()


def sequences_differences_segments(seq1, seq2):
    """Return the list of segments on which sequence seq1 differs from seq2.

    The list is of the form [(start1, end1), (start2, end2), etc.]

    Parameters
    ----------

    seq1, seq2
      ATGC sequences to be compared
    """
    arr = 1 * sequences_differences_array(seq1, seq2)
    diffs = np.diff([0] + list(arr) + [0]).nonzero()[0]
    half = int(len(diffs) / 2)
    return [(diffs[2 * i], diffs[2 * i + 1]) for i in range(half)]

def group_nearby_indices(indices, max_gap=None, max_group_spread=None):
    """Return a list of groups of the different indices.

    Indices are considered from smaller to larger and placed into groups

    Parameters
    ----------
    max_gap
      Maximal allowed difference between two consecutive numbers of a group

    max_group_spread
      Maximal allowed difference between the smallest and largest elements
      of a group.
    """
    if len(indices) == 0:
        return []
    indices = sorted(indices)
    current_group = [indices[0]]
    groups = [current_group]
    for ind in indices[1:]:
        gap_small_enough = ((max_gap is None) or
                            (ind - current_group[-1] < max_gap))
        spread_small_enough = ((max_group_spread is None) or
                               (ind - current_group[0] < max_group_spread))
        if gap_small_enough and spread_small_enough:
            current_group.append(ind)
        else:
            current_group = [ind]
            groups.append(current_group)
    return groups

def group_nearby_segments(segments, max_start_gap=None, max_start_spread=None):
    """Return a list of groups of the different indices.

    Indices are considered from smaller to larger and placed into groups

    Parameters
    ----------
    max_gap
      Maximal allowed difference between two consecutive numbers of a group

    max_group_spread
      Maximal allowed difference between the smallest and largest elements
      of a group.
    """
    if len(segments) == 0:
        return []
    segments = sorted(segments)
    current_group = [segments[0]]
    groups = [current_group]
    for seg in segments[1:]:
        gap_small_enough = ((max_start_gap is None) or
                            (seg[0] - current_group[-1][0] < max_start_gap))
        spread_small_enough = ((max_start_spread is None) or
                               (seg[0] - current_group[0][0] <
                                   max_start_spread))
        if gap_small_enough and spread_small_enough:
            current_group.append(seg)
        else:
            current_group = [seg]
            groups.append(current_group)
    return groups

def codons_frequencies_and_positions(sequence):
    """Return dicts indicating codons frequencies and positions.

    Parameters
    ----------
    sequence
      A sequence with length divisible by 3 (like an ORF)
    
    Returns
    -------
    codons_frequencies
      A dict ``{'K': {'total': 60, 'AAA': 0.5, 'AAG': 0.4}}, 'M': {...} ...}``
      providing the frequency of each version of each codon in the given
      sequence.
    
    
    codons_positions
      A dict ``{'ATA': [2, 55, 1], 'ATG': [0] ...}`` indicating the position of
      the different codons in the sequence.

    """
    codons_positions = {cod: [] for cod in CODONS_TRANSLATIONS}
    for i in range(int(len(sequence) / 3)):
        codon = sequence[3 * i: 3 * (i + 1)]
        codons_positions[codon].append(i)
    # aa: amino-acid
    codons_frequencies = {aa: {'total': 0} for aa in CODONS_SEQUENCES}
    for codon, positions in codons_positions.items():
        count = len(positions)
        aa = CODONS_TRANSLATIONS[codon]
        codons_frequencies[aa][codon] = count
        codons_frequencies[aa]['total'] += count
    for aa, data in codons_frequencies.items():
        total = max(1, data['total'])
        for codon, value in data.items():
            if codon != 'total':
                data[codon] = 1.0 * value / total
    return codons_frequencies, codons_positions

def all_iupac_variants(iupac_sequence):
    """Return all unambiguous possible versions of the given sequence."""
    return[
        "".join(nucleotides)
        for nucleotides in itertools.product(*[
            IUPAC_NOTATION[n] for n in iupac_sequence
        ])
    ]
def list_common_enzymes(site_length=(6,), opt_temp=(37,), min_suppliers=1,
                        site_unlike=()):
    """Return a list of enzyme names with the given constraints.
    
    Parameters
    ----------

    site_length
      List of accepted site lengths (6, 4, ...)
    
    opt_temp
      List of accepted optimal temperatures for the enzyme
    
    min_suppliers
      Minimal number registered suppliers in the Biopython data. A minimum
      of 3 known suppliers returns the most common enzymes.

    site_unlike
      List of (ambiguous or unambiguous) DNA sequences that should NOT be
      recognized by the selected enzymes.
    """
    site_unlike = set([
        variant
        for enzyme in site_unlike
        for variant in all_iupac_variants(Restriction.__dict__[enzyme].site)
    ])
    def is_valid(enzyme_name):
        enzyme = Restriction.__dict__[enzyme_name]
        return (len(enzyme.site) in site_length
                and enzyme.opt_temp in opt_temp
                and len(enzyme.supplier_list()) >= min_suppliers
                and len(set(all_iupac_variants(enzyme.site))
                        .intersection(site_unlike)) == 0)
    return [
        enzyme_name
        for enzyme_name in Restriction.AllEnzymes.elements()
        if is_valid(enzyme_name)
    ]

# def crop_record(record, crop_start, crop_end, features_suffix=" (part)"):
#     """Return the cropped record with possibly cropped features.

#     Note that this differs from ``record[start:end]`` in that in the latter
#     expression, cropped features are discarded.

#     Parameters
#     ----------

#     record
#       A Biopython record

#     crop_start, crop_end
#       Start and end of the segment to be cropped.

#     features_suffix
#       All cropped features will have their label appended with this suffix.
#     """
#     features = []
#     for feature in record.features:
#         start, end = sorted([feature.location.start, feature.location.end])
#         new_start, new_end = max(start, crop_start), min(end, crop_end)
#         if new_end <= new_start:
#             continue
#         new_start, new_end = new_start - crop_start, new_end - crop_start

#         feature = deepcopy(feature)
#         feature.location = FeatureLocation(int(new_start), int(new_end),
#                                            feature.location.strand)
#         label = "".join(feature.qualifiers.get("label", ""))
#         feature.qualifiers["label"] = label + features_suffix
#         features.append(feature)

#     new_record = record[crop_start: crop_end]
#     new_record.features = features
#     return new_record