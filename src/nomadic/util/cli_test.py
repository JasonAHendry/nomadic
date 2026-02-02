import pytest

from nomadic.util.cli import get_parameter_name_from_hint


@pytest.mark.parametrize(["hint", "expected"], [["-k/--minknow-dir", "minknow-dir"]])
def test_get_parameter_name_from_hint(hint, expected):
    assert get_parameter_name_from_hint(hint) == expected
