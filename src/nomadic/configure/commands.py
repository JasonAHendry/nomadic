import os
from pathlib import Path
import click

from nomadic.util.config import (
    load_config,
    default_config_path,
    set_config_value,
    write_config,
)
from nomadic.util.ssh import is_ssh_target
from nomadic.util.workspace import check_if_workspace


@click.group(short_help="Configure nomadic settings")
def configure():
    """
    Configure different nomadics functionality. This mostly sets standard options in '.config.yaml' that can be overwritten from the command line.
    """
    pass


@configure.command(short_help="Configure the nomadic share command.")
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the workspace where all input/output files (beds, metadata, results) are stored. "
    "The workspace directory simplifies the use of nomadic in that many arguments don't need to be listed "
    "as they are predefined in the workspace config or can be loaded from the workspace",
)
@click.option(
    "-t",
    "--target_dir",
    "target_dir",
    type=click.Path(
        dir_okay=True, file_okay=False
    ),  # Use path here to enable autocomplete, even though not only paths are allowed
    help=(
        "Path to target folder or an SSH target like user@host:/path. "
        "The shared files will go inside of that folder into a folder with the name of the workspace."
    ),
    prompt="Set the target into which to share",
)
def share(workspace_path: str, target_dir: str):
    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Workspace '{workspace_path}' does not exist or is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
        )

    config_path = os.path.join(workspace_path, default_config_path)
    if os.path.isfile(config_path):
        config = load_config(config_path)
    else:
        config = {}

    if is_ssh_target(target_dir):
        to_store = target_dir
    else:
        target_path = Path(target_dir)
        if not target_path.is_dir():
            raise click.BadParameter(
                param_hint="-t/--target_dir",
                message=f"Directory {target_dir} could not be found",
            )
        to_store = str(target_path.resolve())

    set_config_value(config, ["share", "defaults", "target_dir"], to_store)

    write_config(config, config_path)


@configure.command(short_help="Configure the nomadic backup command.")
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the workspace where all input/output files (beds, metadata, results) are stored. "
    "The workspace directory simplifies the use of nomadic in that many arguments don't need to be listed "
    "as they are predefined in the workspace config or can be loaded from the workspace",
)
@click.option(
    "-t",
    "--target_dir",
    "target_dir",
    type=click.Path(
        dir_okay=True, file_okay=False
    ),  # Use path here to enable autocomplete, even though not only paths are allowed
    help=(
        "Path to target folder or an SSH target like user@host:/path. "
        "The shared files will go inside of that folder into a folder with the name of the workspace."
    ),
    prompt="Set the target into which to share",
)
def backup(workspace_path: str, target_dir: str):
    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Workspace '{workspace_path}' does not exist or is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
        )

    config_path = os.path.join(workspace_path, default_config_path)
    if os.path.isfile(config_path):
        config = load_config(config_path)
    else:
        config = {}

    if is_ssh_target(target_dir):
        to_store = target_dir
    else:
        target_path = Path(target_dir)
        if not target_path.is_dir():
            raise click.BadParameter(
                param_hint="-t/--target_dir",
                message=f"Directory {target_dir} could not be found",
            )
        to_store = str(target_path.resolve())

    set_config_value(config, ["backup", "defaults", "target_dir"], to_store)

    write_config(config, config_path)
