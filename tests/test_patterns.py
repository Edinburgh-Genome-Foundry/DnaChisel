from dnachisel.SequencePattern import SequencePattern

def test_patterns_from_string():
    pattern = SequencePattern.from_string("6xT")
    assert pattern.expression == "TTTTTT"
    pattern = SequencePattern.from_string("BsmBI_site")
    assert pattern.expression == "CGTCTC"
    pattern = SequencePattern.from_string("5x2mer")
    assert pattern.expression == '([ATGC]{2})\\1{4}'