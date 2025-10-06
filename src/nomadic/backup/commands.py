from pathlib import Path

import click

from nomadic.util import minknow
from nomadic.util.dirs import produce_dir
from nomadic.util.rsync import selective_rsync
from nomadic.util.settings import load_settings
from nomadic.util.workspace import check_if_workspace


@click.command(
    short_help="Backup entire nomadic workspace and associated minknow data to a different folder e.g. on a local hard disk drive"
)
@click.option(
    "-b",
    "--backup_dir",
    "backup_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    required=True,
    help="Path to backup folder",
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Path to the nomadic workspace you want to back up",
)
def backup(backup_dir: Path, workspace_path: Path):
    """
    Backup nomadic workspace and associated minknow data to another location.
    """

    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"'{workspace_path.resolve()}' is not a workspace.",
        )

    click.echo(f"Backing up nomadic workspace from {workspace_path} to {backup_dir}")
    selective_rsync(
        source_dir=workspace_path,
        target_dir=backup_dir,
        recursive=True,
        progressbar=True,
    )
    click.echo("Done.")

    exp_dirs = [f.name for f in (workspace_path / "results").iterdir() if f.is_dir()]
    for folder in exp_dirs:
        click.echo(folder)
        json_file = workspace_path / "results" / folder / "metadata" / "settings.json"
        settings = load_settings(json_file)

        if settings.minknow_dir is not None:
            minknow_dir = Path(settings.minknow_dir).resolve()
        elif settings.fastq_dir is not None:
            # this could be the minknow or the fastq_pass dir so first resolve
            minknow_dir, _ = minknow.resolve_minknow_fastq_dirs(
                Path(settings.fastq_dir), experiment_name=folder
            )
        else:
            click.BadParameter(
                message=f"Unable to determine the minknow directory for experiment {folder}, skipping backup of minknow data for this experiment.",
            )
        source_dir = minknow_dir / folder
        target_dir = backup_dir / "results" / folder / "minknow"

        if not minknow_dir.exists():
            click.echo(f"   ERROR: {source_dir} does not exist, unable to backup...")
            continue
        if not target_dir.exists():
            produce_dir(target_dir)
        click.echo(f"   {source_dir} to {target_dir}")
        selective_rsync(
            source_dir=source_dir,
            target_dir=target_dir,
            recursive=True,
            progressbar=True,
        )

    click.echo("Backup complete.")
