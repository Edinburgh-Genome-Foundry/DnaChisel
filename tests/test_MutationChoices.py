from dnachisel.MutationSpace import MutationChoice


def test_mutation_choice_extract_varying_region():
    choice = MutationChoice((5, 12), [
        'ATGCGTG',
        'AAAAATG',
        'AAATGTG',
        'ATGAATG',
    ])
    choices = choice.extract_varying_region()
    assert [c.segment for c in choices] == [(5, 6), (6, 10), (10, 12)]

    choice = MutationChoice((5, 12), [
        'TAGCGTG',
        'AAAAATG',
        'AAATGTG',
        'ATGAATG',
    ])
    choices = choice.extract_varying_region()
    assert [c.segment for c in choices] == [(5, 10), (10, 12)]

    choice = MutationChoice((5, 12), [
        'ATGCGTG',
        'AAAAACC',
        'AAATGTG',
        'ATGAATG',
    ])
    choices = choice.extract_varying_region()
    assert [c.segment for c in choices] == [(5, 6), (6, 12)]
