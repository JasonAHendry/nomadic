import os
from pathlib import Path
from typing import Optional

import click
import click.shell_completion

from nomadic.util.ssh import is_ssh_target, remote_dir_exists
from nomadic.util.workspace import (
    Workspace,
    find_workspace_root,
    check_if_workspace_root,
)
from nomadic.util.config import (
    InvalidConfigError,
    get_command_defaults,
    load_config,
    default_config_path,
)

WORKSPACE_OPTION_KEY = "workspace"


def find_workspace(ctx, param, value) -> Optional[Workspace]:
    """Find the workspace

    If no value is given, checks if the current path is inside a workspace(by walking up the file tree), if not return None
    If a value is given, check if it is a workspace (don't walk up the tree for explicit paths), and if not return an error.
    """
    if value is None:
        root = find_workspace_root(Path("./"))
        if root is None:
            return None
    elif check_if_workspace_root(value):
        root = Path(value)
    else:
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Directory {value} is not a workspace. Please use nomadic start to create a new workspace, or point to a valid workspace.",
        )

    return Workspace(str(root))


def must_find_workspace(ctx, param, value) -> Workspace:
    workspace = find_workspace(ctx, param, value)
    if workspace is None:
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message="The current directory is not a workspace. Please use nomadic start to create a new workspace, or point to a valid workspace.",
        )
    return workspace


def workspace_option(optional=False):
    return click.option(
        "-w",
        "--workspace",
        WORKSPACE_OPTION_KEY,
        show_default="current directory",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        help="Path of the workspace where all input/output files (beds, metadata, results) are stored. "
        "The workspace directory simplifies the use of nomadic in that many arguments don't need to be listed "
        "as they are predefined in the workspace config or can be loaded from the workspace",
        callback=find_workspace if optional else must_find_workspace,
    )


def complete_experiment_name(ctx: click.Context, param, incomplete):
    """Complete experiment names based on existing metadatafiles in the workspace."""
    workspace: Optional[Workspace] = ctx.params.get(WORKSPACE_OPTION_KEY, None)
    if workspace is None:
        return []
    metadata_path = workspace.get_metadata_dir()

    experiments = []
    if os.path.exists(metadata_path):
        experiments.extend(list_experiment_names(metadata_path))
    shared_workspace = workspace.get_shared_workspace()
    if shared_workspace:
        shared_metadata_path = os.path.join(
            shared_workspace, Workspace(shared_workspace).get_metadata_dir()
        )
        if os.path.exists(shared_metadata_path):
            experiments.extend(list_experiment_names(shared_metadata_path))
    experiments = sorted(list(set(experiments)))  # Remove duplicates and sort
    return [
        click.shell_completion.CompletionItem(experiment)
        for experiment in experiments
        if experiment.startswith(incomplete)
    ]


def list_experiment_names(metadata_folder: str) -> list[str]:
    return [
        f.removesuffix(".csv").removesuffix(".xlsx")
        for f in os.listdir(metadata_folder)
        if f.endswith(".csv") or f.endswith(".xlsx")
    ]


def complete_bed_file(ctx: click.Context, param, incomplete):
    """Complete bed file options based on existing BED files in the workspace."""
    workspace = ctx.params.get(WORKSPACE_OPTION_KEY, None)
    if workspace is None:
        return [click.shell_completion.CompletionItem(incomplete, type="file")]
    assert isinstance(workspace, Workspace)
    result = []
    if os.path.exists(workspace.get_beds_dir()):
        panels = workspace.get_panel_names()
        result = [
            click.shell_completion.CompletionItem(panel)
            for panel in panels
            if panel.startswith(incomplete)
        ]
    if not result:
        return [click.shell_completion.CompletionItem(incomplete, type="file")]
    return result


def load_defaults_from_config(ctx: click.Context, command: Optional[str] = None):
    """Load configuration from the default config file if it exists."""
    workspace = ctx.params.get(WORKSPACE_OPTION_KEY, None)
    if workspace is None:
        return
    assert isinstance(workspace, Workspace)
    config_path = os.path.join(workspace.path, default_config_path)
    if os.path.isfile(config_path):
        config = load_config(config_path)
        if not config:
            # When empty dict or non, the config is empty
            return
        if not isinstance(config, dict):
            if not ctx.resilient_parsing:
                raise click.UsageError(
                    f"Invalid config at {config_path}: config is not a dict."
                )
            else:
                return
        try:
            defaults = get_command_defaults(config, command)
        except InvalidConfigError as e:
            if not ctx.resilient_parsing:
                raise click.UsageError(f"Invalid config at {config_path}: {e}")
            else:
                return

        if defaults:
            if not ctx.resilient_parsing:
                # Don't print defaults if parsing is resilient, as this is used for shell completion
                click.echo(f"Loaded defaults from {config_path}")
            ctx.default_map = defaults
            ctx.show_default = True


def load_default_function_for(command: str):
    def click_callback(ctx: click.Context, param, value):
        load_defaults_from_config(ctx, command=command)
        return value

    return click_callback


def validate_target(ctx, param, value) -> Path | str:
    if is_ssh_target(value):
        ok, msg = remote_dir_exists(value)
        if not ok:
            raise BadParameterWithSource(
                message=f"Remote target '{value}' does not exist: {msg}",
            )
        return value
    else:
        path = Path(value)
        if not path.exists():
            raise BadParameterWithSource(
                message=f"'{path}' does not exist.",
            )
        if not path.is_dir():
            raise BadParameterWithSource(
                message=f"'{path}' is not a directory.",
            )
        return path


class BadParameterWithSource(click.BadParameter):
    def format_message(self) -> str:
        if (
            self.param is not None
            and self.ctx is not None
            and self.param.name is not None
            and self.ctx.get_parameter_source(self.param.name)
            is click.core.ParameterSource.DEFAULT_MAP
        ):
            return "Invalid config value for {param_name}: {message}".format(
                param_name=self.param.name, message=self.message
            )

        if (
            self.param_hint is not None
            and self.ctx is not None
            and self.ctx.get_parameter_source(
                get_parameter_name_from_hint(self.param_hint)
            )
            is click.core.ParameterSource.DEFAULT_MAP
        ):
            return "Invalid config value for {param_name}: {message}".format(
                param_name=get_parameter_name_from_hint(self.param_hint),
                message=self.message,
            )

        return super().format_message()


def get_parameter_name_from_hint(hint: str) -> str:
    """Parse a parameter hint like '-k/--minknow-dir'"""
    index = hint.find("--")
    return hint[index + 2 :]
