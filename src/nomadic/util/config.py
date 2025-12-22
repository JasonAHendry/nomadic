from functools import reduce
from typing import Optional

from yaml import dump, load

try:
    from yaml import CDumper as Dumper
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Dumper, Loader

default_config_path = ".config.yaml"


class InvalidConfigError(Exception):
    """Raised when the config is invalid"""

    pass


def load_config(config_path: str) -> dict:
    """
    Load configuration from a YAML file.

    Args:
        config_path (str): Path to the YAML configuration file.

    Returns:
        dict: Configuration data loaded from the file.
    """
    with open(config_path, "r") as f:
        config = load(f.read(), Loader=Loader)
    return config


def write_config(config: dict, config_path: str) -> None:
    """
    Write configuration to a YAML file.

    Args:
        config (dict): Configuration data to write.
    """
    with open(config_path, "w") as f:
        dump(config, f, Dumper=Dumper, default_flow_style=False)


def set_config_value(dict: dict, path: list[str], value):
    for key in path[:-1]:
        if key not in dict:
            dict[key] = {}
        dict = dict[key]

    dict[path[-1]] = value


def get_config_value(d: dict, keys: list, default=None):
    """
    Identify a nested value in a dictionary given a list of keys.

    Args:
        d (dict): The dictionary to search.
        keys (list): A list of keys representing the sequential path to the desired value.
        default: The value to return if the path does not exist.

    Returns:
        The value found at the nested path, or the default value if the path does not exist
    """
    return (
        reduce(lambda acc, k: acc.get(k, {}) if isinstance(acc, dict) else {}, keys, d)
        or default
    )


def get_command_defaults(config: dict, command: Optional[str]) -> dict:
    defaults = must_get_dict(config, "defaults", "defaults should be a dict")
    if command is not None:
        command_config = must_get_dict(
            config, command, f"{command} config should be a dict"
        )
        command_defaults = must_get_dict(
            command_config, "defaults", f"{command} defaults should be a dict"
        )
        defaults = defaults | command_defaults
    return defaults


def must_get_dict(d: dict, key: str, message: str) -> dict:
    value = d.get(key, {})
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise InvalidConfigError(message)
    return value
