from dnachisel.MutationSpace import MutationChoice, MutationSpace


def test_mutation_choice_extract_varying_region():
    choice = MutationChoice((5, 12), ["ATGCGTG", "AAAAATG", "AAATGTG", "ATGAATG"])
    choices = choice.extract_varying_region()
    assert [c.segment for c in choices] == [(5, 6), (6, 10), (10, 12)]

    choice = MutationChoice((5, 12), ["TAGCGTG", "AAAAATG", "AAATGTG", "ATGAATG"])
    choices = choice.extract_varying_region()
    assert [c.segment for c in choices] == [(5, 10), (10, 12)]

    choice = MutationChoice((5, 12), ["ATGCGTG", "AAAAACC", "AAATGTG", "ATGAATG"])
    choices = choice.extract_varying_region()
    assert [c.segment for c in choices] == [(5, 6), (6, 12)]


def test_mutation_choice_io():
    choice = MutationChoice((5, 12), ["ATGCGTG", "AAAAACC", "AAATGTG", "ATGAATG"])
    assert "ATGCGTG-AAAAACC-AAATGTG-ATGAATG" in str(choice)


def test_mutation_space():
    loc1, seqs1 = (0, 2), ["AT", "TG"]
    loc2, seqs2 = (2, 5), ["TTC", "TTA", "TTT"]
    c1 = MutationChoice(loc1, seqs1)
    c2 = MutationChoice(loc2, seqs2)
    space = MutationSpace([c1, c1, c2, c2, c2])
    assert space.string_representation() == "AT|TTC\nTG|TTA\n  |TTT"

    loc, seq = space.pick_random_mutations(n_mutations=1, sequence="ATTTC")[0]
    assert loc in [loc1, loc2]
    assert seq in (seqs1 + seqs2)
