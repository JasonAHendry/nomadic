from collections import defaultdict
from pathlib import Path

import click

from nomadic.util.rsync import (
    backup_minknow_data,
    backup_nomadic_workspace,
    print_rsync_summary,
    rsync_status,
)
from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(short_help="Backup a workspace.")
@click.option(
    "-b",
    "--backup_dir",
    "backup_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    required=True,
    help="Path to root backup folder. The backup will go into backup_dir/<workspace_name>.",
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Path to the nomadic workspace you want to back up.",
)
@click.option(
    "-k",
    "--minknow_dir",
    "minknow_base_dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default="/var/lib/minknow/data",
    show_default=True,
    help="Path to the base minknow output directory. Only needed if the files were moved.",
)
@click.option(
    "--include-minknow/--exclude-minknow",
    "include_minknow",
    default=True,
    show_default=True,
    help="Include/exclude minknow data in the backup.",
)
@click.pass_context
def backup(
    ctx: click.Context,
    backup_dir: Path,
    workspace_path: Path,
    include_minknow: bool,
    minknow_base_dir: Path,
):
    """
    Backup entire nomadic workspace and associated minknow data to a different folder e.g. on a local hard disk drive.
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
    workspace_name = workspace.get_name()
    backup_dir = backup_dir / workspace_name

    failure_reasons = defaultdict(list)

    backup_nomadic_workspace(target_dir=backup_dir, workspace=workspace)

    if include_minknow:
        backup_minknow_data(
            target_base_dir=backup_dir,
            minknow_base_dir=minknow_base_dir,
            workspace=workspace,
            failure_reasons=failure_reasons,
        )
    else:
        click.echo("Skipping minknow data backup as requested.")

    all_backed_up, status_by_exp = rsync_status(
        backup_dir, workspace, include_minknow=include_minknow
    )
    print_rsync_summary(all_backed_up, status_by_exp, failure_reasons, include_minknow)
