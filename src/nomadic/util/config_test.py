import pytest

from nomadic.util.config import (
    InvalidConfigError,
    get_command_defaults,
    get_config_value,
    set_config_value,
)


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


@pytest.mark.parametrize(
    "config, command, expected",
    [
        ({}, None, {}),
        ({}, "test", {}),
        ({"defaults": None}, "test", {}),
        ({"defaults": {"foo": "bar"}}, "test", {"foo": "bar"}),
        (
            {"defaults": {"foo": "bar"}, "test": {"defaults": {"this": "that"}}},
            "test",
            {"foo": "bar", "this": "that"},
        ),
    ],
)
def test_get_command_defaults(config, command, expected):
    assert get_command_defaults(config, command) == expected


@pytest.mark.parametrize(
    "config, command",
    [
        ({"defaults": 5}, None),
        ({"defaults": {}, "test": 5}, "test"),
    ],
)
def test_get_command_defaults_invalid_configs(config, command):
    with pytest.raises(InvalidConfigError):
        get_command_defaults(config, command)
