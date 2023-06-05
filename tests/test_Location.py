from dnachisel.Location import Location

import pytest

@pytest.mark.parametrize("location_1,location_2,expected_result", [
        (Location(0, 1), Location(0, 1), True),
        (Location(0, 1), Location(1, 1), False),
        (Location(0, 1), Location(0, 1, 1), False),
        (Location(0, 1, 1), Location(0, 1, 1), True),
        (Location(0, 1, 0), Location(0, 1, 0), True),
        (Location(0, 1), Location(1, 0), False),
        (Location(0, 1), Location(1, 0, 1), False),
        (Location(0, 1), Location(1, 0, 0), False),
    ]
)
def test___eq__(location_1: Location, location_2: Location, expected_result: bool):
    assert (location_1 == location_2) == expected_result

@pytest.mark.parametrize("location_1,location_2,expected_result", [
        (Location(0, 1), Location(0, 1), True),
        (Location(0, 1), Location(1, 1), False),
        (Location(1, 0), Location(0, 1), True),
        (Location(1, 0, 1), Location(1, 0), True),
        (Location(0, 1), Location(0, 1, 1), False),
        (Location(0, 1, 1), Location(0, 1), True),
        (Location(0, 1, 1), Location(0, 1, 1), True),
        (Location(0, 1, 0), Location(0, 1, 0), True),
        (Location(0, 1), Location(1, 0), False),
        (Location(0, 1), Location(1, 0, 1), False),
        (Location(0, 1), Location(1, 0, 0), False),
    ]
)
def test___ge__(location_1: Location, location_2: Location, expected_result: bool):
    assert (location_1 >= location_2) == expected_result

@pytest.mark.parametrize("location_1,location_2,expected_result", [
        (Location(0, 1), Location(0, 1), False),
        (Location(0, 1), Location(1, 1), True),
        (Location(1, 0), Location(0, 1), False),
        (Location(1, 0, 1), Location(1, 0), False),
        (Location(0, 1), Location(0, 1, 1), True),
        (Location(0, 1, 1), Location(0, 1), False),
        (Location(0, 1, 1), Location(0, 1, 1), False),
        (Location(0, 1, 0), Location(0, 1, 0), False),
        (Location(0, 1), Location(1, 0), True),
        (Location(0, 1), Location(1, 0, 1), True),
        (Location(0, 1), Location(1, 0, 0), True),
    ]
)
def test___lt__(location_1: Location, location_2: Location, expected_result: bool):
    assert (location_1 < location_2) == expected_result

@pytest.mark.parametrize("location_1,location_2,expected_result", [
        (Location(0, 1), Location(0, 1), True),
        (Location(0, 1), Location(0, 1, 0), True),
        (Location(0, 1, 0), Location(0, 1), True),
        (Location(0, 1), Location(1, 1), False),
        (Location(1, 0), Location(0, 1), False),
        (Location(1, 0, 1), Location(1, 0), False),
        (Location(0, 1), Location(0, 1, 1), False),
        (Location(0, 1, 1), Location(0, 1), False),
        (Location(0, 1, 1), Location(0, 1, 1), True),
        (Location(0, 1, 0), Location(0, 1, 0), True),
        (Location(0, 1), Location(1, 0), False),
    ]
)
def test__hash__(location_1, location_2, expected_result):
    assert (hash(location_1) == hash(location_2)) == expected_result
