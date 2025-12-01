from collections import defaultdict
from pathlib import Path

import click

from nomadic.util.cli import load_default_function_for, validate_target
from nomadic.util.sync import (
    backup_minknow_data,
    backup_nomadic_workspace,
    print_sync_summary,
    sync_status,
)
from nomadic.util.ssh import is_ssh_target
from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(short_help="Backup a workspace.")
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
    callback=load_default_function_for("backup"),
    expose_value=False,
)
@click.option(
    "-t",
    "--target_dir",
    "target_dir",
    type=str,
    required=True,
    help=(
        "Path to root target backup folder, local or an SSH target like user@host:/path. "
        "The backup will go into <target>/<workspace_name>."
    ),
    callback=validate_target,
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
@click.pass_context
def backup(
    ctx: click.Context,
    target_dir: Path | str,
    workspace_path: Path,
    include_minknow: bool,
    minknow_base_dir: Path,
    verbose: bool,
    checksum: bool,
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

    if isinstance(target_dir, Path):
        target_dir = target_dir / workspace.get_name()
    elif is_ssh_target(target_dir):
        target_dir = f"{target_dir.rstrip('/')}/{workspace.get_name()}"
    else:
        raise click.BadParameter(
            param_hint="-t/--target_dir",
            message=f"Target '{target_dir}' is not a valid local path or SSH target.",
        )

    failure_reasons = defaultdict(list)

    backup_nomadic_workspace(
        target_dir=target_dir,
        workspace=workspace,
        checksum=checksum,
        verbose=verbose,
    )

    if include_minknow:
        backup_minknow_data(
            target_base_dir=target_dir,
            minknow_base_dir=minknow_base_dir,
            workspace=workspace,
            failure_reasons=failure_reasons,
            checksum=checksum,
            verbose=verbose,
        )
    else:
        click.echo("Skipping minknow data backup as requested.")

    all_backed_up, status_by_exp = sync_status(
        target_dir, workspace, include_minknow=include_minknow, verbose=verbose
    )
    print_sync_summary(all_backed_up, status_by_exp, failure_reasons, include_minknow)
