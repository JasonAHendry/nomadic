from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

default_config_path = ".config.yaml"


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
