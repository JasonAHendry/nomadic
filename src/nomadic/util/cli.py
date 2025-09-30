import os

import click
import click.shell_completion

from nomadic.util.workspace import Workspace
from nomadic.util.config import load_config, default_config_path


def complete_experiment_name(ctx: click.Context, param, incomplete):
    """Complete experiment names based on existing metadatafiles in the workspace."""
    workspace_path = ctx.params.get("workspace_path", "./")
    metadata_path = Workspace(workspace_path).get_metadata_dir()
    if os.path.exists(metadata_path):
        experiments = [
            f.removesuffix(".csv")
            for f in os.listdir(metadata_path)
            if f.endswith(".csv")
        ]
        return [
            click.shell_completion.CompletionItem(experiment)
            for experiment in experiments
            if experiment.startswith(incomplete)
        ]
    return []


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


def load_defaults_from_config(ctx: click.Context, param, value):
    """Load configuration from the default config file if it exists."""
    config_path = os.path.join(
        ctx.params.get("workspace_path", "./"), default_config_path
    )
    if os.path.isfile(config_path):
        defaults = load_config(config_path).get("defaults", None)
        if defaults is not None:
            if not ctx.resilient_parsing:
                # Don't print defaults if parsing is resilient, as this is used for shell completion
                click.echo(f"Loaded defaults from {config_path}")
            ctx.default_map = defaults
            ctx.show_default = True
    return value
