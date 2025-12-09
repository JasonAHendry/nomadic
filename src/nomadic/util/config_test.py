import pytest

from nomadic.util.config import get_config_value, set_config_value


@pytest.mark.parametrize(
    "d, keys, default, expected",
    [
        ({"a": {"b": {"c": 1}}}, ["a", "b", "c"], None, 1),
        ({"a": {"b": {}}}, ["a", "b", "c"], "not-found", "not-found"),
        ({}, ["a", "b"], None, None),
        ({}, ["a", "b"], "default", "default"),
    ],
)
def test_get_config_value_table(d, keys, default, expected):
    assert get_config_value(d, keys, default=default) == expected


@pytest.mark.parametrize(
    "keys, value, expected",
    [
        (["a", "b", "c"], 1, {"a": {"b": {"c": 1}}}),
        (["a"], 1, {"a": 1}),
    ],
)
def test_set_config_value_table(keys, value, expected):
    d = {}
    set_config_value(d, keys, value)
    assert get_config_value(d, keys) == value
    assert d == expected
