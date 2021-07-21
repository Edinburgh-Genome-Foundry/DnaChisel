"""Generic methods for reading/modifying Genbank/Biopython records"""

from copy import deepcopy

import numpy as np
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from Bio import SeqIO

try:
    from snapgene_reader import snapgene_file_to_seqrecord
except ImportError:

    def snapgene_file_to_seqrecord(*a, **k):
        """Please install the snapgene_reader library to use this function."""
        raise ImportError(
            "Please install snapgene_reader to import Snapgene .dna files"
        )


def load_record(filepath, linear=True, name="unnamed", file_format="auto"):
    """Load a FASTA/Genbank/Snapgene record.

    Note that reading Snapgene records requires the library snapgene_reader
    installed.
    """
    if file_format != "auto":
        record = SeqIO.read(filepath, file_format)
    elif filepath.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filepath, "genbank")
    elif filepath.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(filepath, "fasta")
    elif filepath.lower().endswith(".dna"):
        record = snapgene_file_to_seqrecord(filepath)
    else:
        raise ValueError("Unknown format for file: %s" % filepath)
    record.linear = linear
    if name != "unnamed":
        record.id = name
        record.name = name.replace(" ", "_")[:20]
    return record


def annotate_record(
    seqrecord, location="full", feature_type="misc_feature", margin=0, **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1).

    feature_type
      The type associated with the feature.

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type,
        )
    )


def annotate_differences(record, reference, feature_type="misc_feature", prefix="#"):
    """Annotate differences between two records in a new record.

    Returns a version of SeqRecord ``record`` where differences with the
    references are annotated as new features.

    Parameters
    ----------
    record
      The SeqRecord to be compared to the reference.

    reference
      The reference SeqRecord. Must be the same size as ``reference``.

    feature_type
      The type of the features added to mark differences.

    prefix
      Each new feature will be labeled "po" where p is the prefix and o the
      original sequence at the feature's location. For instance "#A" or "#TT".

    """
    seq1 = str(record.seq)
    seq2 = str(reference.seq)
    indices_diff = (
        np.frombuffer(seq1.encode(), dtype="uint8")
        - np.frombuffer(seq2.encode(), dtype="uint8")
    ).nonzero()[0]
    indices_diff = [int(e) for e in indices_diff]
    locations = [[indices_diff[0], indices_diff[0]]]
    for ind in indices_diff[1:]:
        if ind - locations[-1][-1] == 1:
            locations[-1][-1] = ind
        else:
            locations.append([ind, ind])
    new_record = deepcopy(record)
    for (start, end) in locations:
        annotate_record(
            new_record,
            location=(start, end + 1),
            feature_type=feature_type,
            label=prefix + seq2[start : end + 1],
        )
    return new_record


def annotate_pattern_occurrences(
    record, pattern, feature_type="misc_feature", prefix="!"
):
    """Return a new record annotated w. all occurences of pattern in sequence.

    Parameters
    -----------
    record
      A Biopython record.

    pattern
      A DnaChisel SequencePattern object (such as DnaPAttern).

    feature_type
      Type of the annotations in the returned record.
    """
    new_record = deepcopy(record)
    label = prefix + str(pattern)
    for location in pattern.find_matches(str(record.seq)):
        annotate_record(
            new_record,
            location=(location.start, location.end),
            feature_type=feature_type,
            label=label,
        )
    return new_record


def change_biopython_record_sequence(record, new_seq):
    """Return a version of the record with the sequence set to new_seq."""
    new_record = deepcopy(record)

    if has_dna_alphabet:
        seq = Seq(new_seq, alphabet=DNAAlphabet())
    else:
        seq = Seq(new_seq)

    new_record.seq = seq
    return new_record


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    if has_dna_alphabet:
        seq = Seq(sequence, alphabet=DNAAlphabet())
    else:
        seq = Seq(sequence)

    return SeqRecord(
        seq=seq,
        id=id,
        name=name,
        features=list(features),
        annotations={"molecule_type": "DNA"},
    )


def find_specification_label_in_feature(feature):
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
        if (potential_label != "") and (potential_label[0] in "@~"):
            return potential_label
    return None


def write_record(
    record,
    target,
    file_format="genbank",
    remove_locationless_features=True,
    max_name_length=20,
):
    """Write a record as genbank, fasta, etc. via Biopython, with fixes.

    Parameters
    ----------
    record
      A biopython record.

    target
      Path to a file or filelike object.

    file_format
      Format, either Genbank or fasta.

    remove_locationless_features
      If True, will remove all features whose location is None, to avoid a
      Biopython bug

    max_name_length
      The record's name will be truncated if longer than this (also here to
      avoid a Biopython bug).

    """
    record = deepcopy(record)
    if remove_locationless_features:
        record.features = [f for f in record.features if f.location is not None]
    record.name = record.name[:max_name_length]
    if has_dna_alphabet:
        if str(record.seq.alphabet.__class__.__name__) != "DNAAlphabet":
            record.seq.alphabet = DNAAlphabet()
    if hasattr(target, "open"):
        target = target.open("w")
    SeqIO.write(record, target, file_format)
