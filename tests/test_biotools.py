import os
from Bio.Data import CodonTable
from Bio.Seq import Seq
from dnachisel.biotools import (
    dna_pattern_to_regexpr,
    change_biopython_record_sequence,
    subdivide_window,
    sequence_to_biopython_record,
    annotate_record,
    load_record,
    write_record,
    translate,
    list_common_enzymes,
    reverse_translate,
)

data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")


def test_load_and_write_record(tmpdir):
    record = load_record(os.path.join(data_dir, "example_sequence.gbk"))
    assert len(record.seq) == 4720
    target = os.path.join(str(tmpdir), "example.gb")
    write_record(record, target)
    record2 = load_record(target)
    assert len(record2.seq) == 4720


def test_dna_pattern_to_regexpr():
    assert dna_pattern_to_regexpr("ATW") == "AT[ATW]"


def test_subdivide_window():
    windows = subdivide_window((0, 10), max_span=3)
    assert windows == [(0, 3), (3, 6), (6, 9), (9, 10)]


def test_change_biopython_record_sequence():
    record = sequence_to_biopython_record("ATGCATGCATGC")
    annotate_record(record, (0, 5), label="my_label")
    new_record = change_biopython_record_sequence(record, "GGCCGGCCGGCCGGCC")
    assert len(new_record.features) == 1
    assert new_record.features[0].location == record.features[0].location


def test_list_common_enzymes():
    assert len(list_common_enzymes(min_suppliers=3)) == 63


def biopython_translate(sequence, table="Standard", cds=False):
    table = CodonTable.unambiguous_dna_by_name[table]
    return str(Seq(sequence).translate(table=table, cds=cds))


def test_translate_standard_table():
    sequence = "ATGATTTTAATC"
    tables = ["Standard", "Mold Mitochondrial", "Bacterial"]
    for t in tables:
        protein = translate(sequence, table=t)
        biopython_protein = biopython_translate(sequence, table=t)
        assert protein == biopython_protein


def test_translate_bacterial_valine():
    sequence = "GTGATTGTGTGA"
    assert translate(sequence) == "VIV*"

    assert biopython_translate(sequence, "Bacterial", cds=False) == "VIV*"
    assert translate(sequence, "Bacterial") == "VIV*"

    assert translate(sequence, "Bacterial", assume_start_codon=True) == "MIV*"
    assert biopython_translate(sequence, "Bacterial", cds=True) == "MIV"


def reverse_translate_with_tables():
    protein = "MAWWL*"
    for table, expected in [
        ("Standard", "ATGGCTTGGTGGTTATAA"),
        ("Mycoplasma", "ATGGCTTGATGATTATAA"),
    ]:
        rv_translation = reverse_translate(protein, table=table)
        assert rv_translation == expected
        assert translate(rv_translation, table=table) == protein
