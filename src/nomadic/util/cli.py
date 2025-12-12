import os
from pathlib import Path
from typing import Optional

import click
import click.shell_completion

from nomadic.util.ssh import is_ssh_target, remote_dir_exists
from nomadic.util.workspace import Workspace
from nomadic.util.config import load_config, default_config_path


def complete_experiment_name(ctx: click.Context, param, incomplete):
    """Complete experiment names based on existing metadatafiles in the workspace."""
    workspace_path = ctx.params.get("workspace_path", "./")
    metadata_path = Workspace(workspace_path).get_metadata_dir()
    experiments = []
    if os.path.exists(metadata_path):
        experiments.extend(list_experiment_names(metadata_path))
    shared_folder = Workspace(workspace_path).get_shared_folder()
    if shared_folder:
        shared_metadata_path = os.path.join(
            shared_folder, Workspace(shared_folder).get_metadata_dir()
        )
        if os.path.exists(shared_metadata_path):
            experiments.extend(list_experiment_names(shared_metadata_path))
    experiments = list(set(experiments))  # Remove duplicates
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
    workspace_path = ctx.params.get("workspace_path", "./")
    workspace = Workspace(workspace_path)
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
    config_path = os.path.join(
        ctx.params.get("workspace_path", "./"), default_config_path
    )
    if os.path.isfile(config_path):
        config = load_config(config_path)
        defaults = config.get("defaults", {})
        defaults = defaults | config.get(command, {}).get("defaults", {})
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
            raise click.BadParameter(
                message=f"Remote target '{value}' does not exist: {msg}",
            )
        return value
    else:
        path = Path(value)
        if not path.exists():
            raise click.BadParameter(
                message=f"'{path}' does not exist.",
            )
        if not path.is_dir():
            raise click.BadParameter(
                message=f"'{path}' is not a directory.",
            )
        return path
