from pathlib import Path

import click

from nomadic.util import minknow
from nomadic.util.dirs import produce_dir
from nomadic.util.rsync import selective_rsync
from nomadic.util.settings import load_settings
from nomadic.util.workspace import check_if_workspace, Workspace


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
@click.option(
    "-k",
    "--minknow_dir",
    "minknow_base_dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default="/var/lib/minknow/data",
    show_default=True,
    help="Path to the base minknow output directory.",
)
@click.option(
    "--include-minknow/--exclude-minknow",
    "include_minknow",
    default=True,
    show_default=True,
    help="Include/exclude minknow data in the backup",
)
def backup(
    backup_dir: Path,
    workspace_path: Path,
    include_minknow: bool,
    minknow_base_dir: Path,
):
    """
    Backup nomadic workspace and associated minknow data to another location.
    """

    if not check_if_workspace(str(workspace_path)):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"'{workspace_path.resolve()}' is not a workspace.",
        )
    workspace = Workspace(str(workspace_path))
    workspace_name = workspace.get_name()
    backup_dir = backup_dir / workspace_name

    click.echo(f"Backing up nomadic workspace ({workspace_name}) to {backup_dir}")

    selective_rsync(
        source_dir=workspace_path,
        target_dir=backup_dir,
        recursive=True,
        progressbar=True,
    )
    click.echo("Done.")

    if not include_minknow:
        click.echo("Skipping minknow data backup as requested.")
        return

    click.echo("Merging minknow data into experiments:")

    for exp in workspace.get_experiment_names():
        click.echo(exp)
        json_file = workspace_path / "results" / exp / "metadata" / "settings.json"
        settings = load_settings(str(json_file))

        if settings is None:
            click.echo("Unable to find settings file, trying to find via minknow_dir")
            minknow_dir = minknow_base_dir / exp
        elif settings.minknow_dir is None:
            click.echo(
                "Minknow dir not stored in settings, trying to find via minknow_dir"
            )
            minknow_dir = minknow_base_dir / exp
        else:
            minknow_dir = Path(settings.minknow_dir)

        source_dir = minknow_dir
        target_dir = backup_dir / "results" / exp / "minknow"

        if not source_dir.exists():
            click.echo(f"   ERROR: {source_dir} does not exist, unable to backup...")
            continue
        if not minknow.is_minknow_experiment_dir(source_dir):
            click.echo(
                f"   ERROR: {source_dir} does not look like a valid minknow experiment directory, unable to backup..."
            )
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
