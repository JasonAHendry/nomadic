from collections import defaultdict
from pathlib import Path

import click

from nomadic.util.rsync import (
    print_rsync_summary,
    rsync_minknow_data,
    rsync_nomadic,
    rsync_status,
)
from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(
    short_help="Share summary nomadic workspace and associated minknow data to a cloud synchronized folder"
)
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
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    default="/var/lib/minknow/data",
    show_default="/var/lib/minknow/data",
    help="Path to minknow output directory",
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
def share(
    shared_dir: Path,
    minknow_base_dir: Path,
    include_minknow: bool,
    workspace_path: Path,
):
    """
    Share summary nomadic workspace and associated minknow data to another location.
    """
    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"'{workspace_path.resolve()}' is not a workspace.",
        )
    workspace = Workspace(str(workspace_path))
    workspace_name = workspace.get_name()
    shared_dir = shared_dir / workspace_name

    failure_reasons = defaultdict(list)

    rsync_nomadic(
        target_dir=shared_dir,
        workspace_path=workspace_path,
        workspace_name=workspace_name,
        additional_exclusions=["barcodes", "summary.fastq.csv"],
    )

    if include_minknow:
        rsync_minknow_data(
            shared_dir,
            workspace_path,
            minknow_base_dir,
            workspace,
            failure_reasons,
            exclusions=[
                "fastq_fail",
                "fastq_pass",
                "other_reports",
                "pod5",
                "sequencing_summary_*.txt",
            ],
        )
    else:
        click.echo("Skipping sharing minknow data as requested.")

    all_backed_up, status_by_exp = rsync_status(
        shared_dir, workspace, include_minknow=include_minknow
    )
    print_rsync_summary(all_backed_up, status_by_exp, failure_reasons, include_minknow)
