from collections import defaultdict
from pathlib import Path

import click

from nomadic.util.rsync import (
    copy_minknow_data,
    copy_nomadic_workspace,
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
def backup(
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
    workspace = Workspace(str(workspace_path))
    workspace_name = workspace.get_name()
    backup_dir = backup_dir / workspace_name

    failure_reasons = defaultdict(list)

    copy_nomadic_workspace(backup_dir, workspace_path, workspace_name)

    if include_minknow:
        copy_minknow_data(
            backup_dir, workspace_path, minknow_base_dir, workspace, failure_reasons
        )
    else:
        click.echo("Skipping minknow data backup as requested.")

    all_backed_up, status_by_exp = rsync_status(
        backup_dir, workspace, include_minknow=include_minknow
    )
    print_rsync_summary(all_backed_up, status_by_exp, failure_reasons, include_minknow)
