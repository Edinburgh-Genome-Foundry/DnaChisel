from dnachisel import (
    random_dna_sequence,
    sequence_to_biopython_record,
    DnaNotationPattern,
)
from dnachisel.biotools import (
    annotate_record,
    annotate_differences,
    annotate_pattern_occurrences,
)


def random_record(length, seed=None):
    return sequence_to_biopython_record(random_dna_sequence(length, seed=seed))


def test_annotate_record():
    random_rec = random_record(100)
    annotate_record(
        random_rec, location=(10, 50), feature_type="CDS", label="my_label"
    )
    assert len(random_rec.features) == 1
    feature = random_rec.features[0]
    assert feature.location.start == 10
    assert feature.location.end == 50
    assert feature.type == "CDS"
    assert feature.qualifiers["label"] == "my_label"


def test_annotate_differences():
    def rotate_nucleotide(n):
        return "ATGC"[("ATGC".index(n) + 1) % 4]

    rec1 = random_record(20, seed=123)
    modified = [1, 4, 5, 6, 9, 10]
    seq2 = "".join(
        [
            (n if i not in modified else rotate_nucleotide(n))
            for i, n in enumerate(rec1.seq)
        ]
    )
    rec2 = sequence_to_biopython_record(seq2)
    annotated = annotate_differences(rec1, rec2)
    assert len(annotated.features) == 3


def test_annotate_pattern_occurrences():
    rec = sequence_to_biopython_record("ATCCAAAAATCCTTTTTTATCCAAAAAA")
    annotated = annotate_pattern_occurrences(rec, DnaNotationPattern("ATCC"))
    assert len(annotated.features) == 3
