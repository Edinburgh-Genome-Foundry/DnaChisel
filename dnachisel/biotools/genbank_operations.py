"""Generic methods for reading/modifying Genbank/Biopython records"""

from copy import deepcopy

import numpy as np
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio import SeqIO


def load_record(filename, linear=True, name="unnamed", fmt="auto"):
    """Load a FASTA/Genbank/... record"""
    if fmt != "auto":
        record = SeqIO.read(filename, fmt)
    elif filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(("fa", "fasta")):
        record = SeqIO.read(filename, "fasta")
    else:
        raise ValueError("Unknown format for file: %s" % filename)
    record.linear = linear
    if name != "unnamed":
        record.id = name
        record.name = name.replace(" ", "_")[:20]
    return record


def annotate_record(
    seqrecord,
    location="full",
    feature_type="misc_feature",
    margin=0,
    **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The biopython seqrecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1)

    feature_type
      The type associated with the feature

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionnary that will be the Biopython feature's `qualifiers` attribute.
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


def annotate_differences(
    record, reference, feature_type="misc_feature", prefix="#"
):
    """Annotate differences between two records in a new record.

    Returns a version of SeqRecord ``record`` where differences with the
    references are annotated as new features.

    Parameters
    ----------
    record
      The SeqRecord to be compared to the reference

    reference
      The reference SeqRecord. Must be the same size as ``reference``

    feature_type
      The type of the features added to mark differences.

    prefix
      Each new feature will be labeled "po" where p is the prefix and o the
      original sequence at the feature's location. For instance "#A" or "#TT".

    """
    seq1 = str(record.seq)
    seq2 = str(reference.seq)
    indices_diff = (
        np.fromstring(seq1, dtype="uint8") - np.fromstring(seq2, dtype="uint8")
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
      A Biopython record

    pattern
      A DnaChisel SequencePattern object (such as DnaPAttern)

    feature_type
      Type of the annotations in the returned record

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
    """Return a version of the record with the sequence set to new_seq"""
    new_record = deepcopy(record)
    new_record.seq = Seq(new_seq, alphabet=DNAAlphabet())
    return new_record


def sequence_to_biopython_record(
    sequence, id="<unknown id>", name="<unknown name>", features=()
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    return SeqRecord(
        Seq(sequence, alphabet=DNAAlphabet()),
        id=id,
        name=name,
        features=list(features),
    )


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
        if (potential_label != "") and (potential_label[0] in "@~"):
            return potential_label
    return None
