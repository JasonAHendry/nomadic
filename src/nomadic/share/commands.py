from collections import defaultdict
from pathlib import Path

import click

from nomadic.util.rsync import (
    print_rsync_summary,
    rsync_status,
    share_minknow_data,
    share_nomadic_workspace,
)
from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(short_help="Share summary nomadic and minknow data to another folder")
@click.option(
    "-s",
    "--shared_dir",
    "shared_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    required=True,
    help="Path to shared folder",
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
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Path to the nomadic workspace you want to share",
)
@click.option(
    "--include-minknow/--exclude-minknow",
    "include_minknow",
    default=True,
    show_default=True,
    help="Include/exclude minknow data in the backup.",
)
@click.pass_context
def share(
    ctx: click.Context,
    shared_dir: Path,
    minknow_base_dir: Path,
    include_minknow: bool,
    workspace_path: Path,
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
    shared_dir = shared_dir / workspace.get_name()

    share_nomadic_workspace(
        target_dir=shared_dir,
        workspace=workspace,
    )

    if include_minknow:
        share_minknow_data(
            target_base_dir=shared_dir,
            minknow_base_dir=minknow_base_dir,
            workspace=workspace,
            failure_reasons=failure_reasons,
        )
    else:
        click.echo("Skipping sharing minknow data as requested.")

    all_backed_up, status_by_exp = rsync_status(
        shared_dir, workspace, include_minknow=include_minknow
    )
    print_rsync_summary(all_backed_up, status_by_exp, failure_reasons, include_minknow)
