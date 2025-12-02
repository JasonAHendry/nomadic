from collections import defaultdict
from pathlib import Path

import click

from nomadic.util.cli import load_default_function_for, validate_target
from nomadic.util.sync import (
    print_sync_summary,
    sync_status,
    share_minknow_data,
    share_nomadic_workspace,
)
from nomadic.util.ssh import is_ssh_target
from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(short_help="Share (lightweight) summary of a workspace.")
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Path to the nomadic workspace you want to share",
)
@click.option(
    callback=load_default_function_for("share"),
    expose_value=False,
)
@click.option(
    "-t",
    "--target_dir",
    "target_dir",
    type=str,
    required=True,
    help=(
        "Path to target folder, local or an SSH target like user@host:/path. "
        "The shared files will go inside of that folder into a folder with the name of the workspace."
    ),
    callback=validate_target,
)
@click.option(
    "-k",
    "--minknow_dir",
    "minknow_base_dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default="/var/lib/minknow/data",
    show_default="/var/lib/minknow/data",
    help="Path to minknow output directory (default it usually sufficient)",
)
@click.option(
    "--include-minknow/--exclude-minknow",
    "include_minknow",
    default=True,
    show_default=True,
    help="Include/exclude minknow data",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Increase logging verbosity, will show files that are copied.",
)
@click.option(
    "-c",
    "--checksum",
    is_flag=True,
    default=False,
    help="Use checksum check instead of file size and modification time.",
)
@click.option(
    "--update/--overwrite",
    "update",
    is_flag=True,
    default=True,
    show_default=True,
    help="Update will not overwrite files that are newer.",
)
@click.pass_context
def share(
    ctx: click.Context,
    target_dir: Path | str,
    minknow_base_dir: Path,
    include_minknow: bool,
    workspace_path: Path,
    checksum: bool,
    update: bool,
    verbose: bool,
):
    """
    Share summary nomadic workspace and associated minknow data to another folder
    e.g. a cloud synchronised folder for sharing.
    """
    if not check_if_workspace(str(workspace_path)):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"'{workspace_path.resolve()}' is not a workspace.",
        )

    if (
        ctx.get_parameter_source("minknow_base_dir")
        is not click.core.ParameterSource.DEFAULT
    ):
        # Only check if the minknow dir exists if not used the default
        # This is because we might not have to use this folder, it is only there for a fallback lookup method.
        if not minknow_base_dir.exists():
            raise click.BadParameter(
                param_hint="-k/--minknow_dir",
                message=f"'{minknow_base_dir.resolve()}' does not exist.",
            )
        if not minknow_base_dir.is_dir():
            raise click.BadParameter(
                param_hint="-k/--minknow_dir",
                message=f"'{minknow_base_dir.resolve()}' is not a directory.",
            )

    workspace = Workspace(str(workspace_path))
    failure_reasons = defaultdict(list)

    # Add workspace name to shared dir so multiple workspaces can be shared to same location
    if isinstance(target_dir, Path):
        target_dir = Path(target_dir) / workspace.get_name()
    elif is_ssh_target(target_dir):
        target_dir = f"{target_dir.rstrip('/')}/{workspace.get_name()}"
    else:
        raise click.BadParameter(
            param_hint="-t/--target_dir",
            message=f"Target '{target_dir}' is not a valid local path or SSH target.",
        )

    share_nomadic_workspace(
        target_dir=target_dir,
        workspace=workspace,
        checksum=checksum,
        verbose=verbose,
        update=update,
    )

    if include_minknow:
        share_minknow_data(
            target_base_dir=target_dir,
            minknow_base_dir=minknow_base_dir,
            workspace=workspace,
            failure_reasons=failure_reasons,
            checksum=checksum,
            verbose=verbose,
            update=update,
        )
    else:
        click.echo("Skipping sharing minknow data as requested.")

    all_backed_up, status_by_exp = sync_status(
        target_dir, workspace, include_minknow=include_minknow, verbose=verbose
    )
    print_sync_summary(all_backed_up, status_by_exp, failure_reasons, include_minknow)
