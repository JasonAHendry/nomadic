from pathlib import Path
from collections import defaultdict

import click

from nomadic.util import minknow
from nomadic.util.dirs import produce_dir
from nomadic.util.rsync import selective_rsync
from nomadic.util.settings import load_settings
from nomadic.util.workspace import check_if_workspace, Workspace


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

    backup_nomadic(backup_dir, workspace_path, workspace_name)

    if include_minknow:
        backup_minknow_data(
            backup_dir, workspace_path, minknow_base_dir, workspace, failure_reasons
        )
    else:
        click.echo("Skipping minknow data backup as requested.")

    all_backed_up, status_by_exp = backup_status(
        backup_dir, workspace, include_minknow=include_minknow
    )
    print_summary(all_backed_up, status_by_exp, failure_reasons)


def backup_nomadic(backup_dir, workspace_path, workspace_name):
    click.echo(f"Backing up nomadic workspace ({workspace_name}) to {backup_dir}")

    selective_rsync(
        source_dir=workspace_path,
        target_dir=backup_dir,
        recursive=True,
        progressbar=True,
    )
    click.echo("Done.")


def backup_minknow_data(
    backup_dir: Path,
    workspace_path: Path,
    minknow_base_dir: Path,
    workspace: Workspace,
    failure_reasons: dict[str, list[str]],
):
    click.echo("Merging minknow data into experiments:")

    for exp in workspace.get_experiment_names():
        click.echo(exp)
        json_file = workspace_path / "results" / exp / "metadata" / "settings.json"
        settings = load_settings(str(json_file))

        if settings is None:
            click.echo(
                "Unable to find settings file, trying to find via minknow_dir..."
            )
            minknow_dir = minknow_base_dir / exp
        elif settings.minknow_dir is None:
            click.echo(
                "Minknow dir not stored in settings.json, trying to find via minknow_dir..."
            )
            minknow_dir = minknow_base_dir / exp
        else:
            minknow_dir = Path(settings.minknow_dir)

        source_dir = minknow_dir
        target_dir = backup_dir / "results" / exp / "minknow"

        if not source_dir.exists():
            click.echo(f"   ERROR: {source_dir} does not exist, unable to backup...")
            failure_reasons[exp].append("no minknow data found")
            continue
        if not minknow.is_minknow_experiment_dir(source_dir):
            failure_reasons[exp].append("invalid minknow data")
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


def backup_status(backup_dir: Path, workspace: Workspace, include_minknow: bool):
    """
    Calculate backup status for all experiments in the workspace.

    This uses a very simple heruitic of checking if the result dir is there and relies
    on rsync working correctly. It is just a sanity check.
    """
    all_backed_up = True
    status_by_exp = defaultdict(list)
    for exp in workspace.get_experiment_names():
        exp_dir = backup_dir / "results" / exp
        nomadic_backed_up = exp_dir.exists()
        status_by_exp[exp].append(nomadic_backed_up)
        if not nomadic_backed_up:
            all_backed_up = False
        if include_minknow:
            minknow_backed_up = (exp_dir / "minknow").exists()
            status_by_exp[exp].append(minknow_backed_up)
            if not minknow_backed_up:
                all_backed_up = False
    return all_backed_up, status_by_exp


def print_summary(all_backed_up, status_by_exp, failure_reasons):
    click.echo("")
    click.echo("--- Summary ---")
    click.echo("")
    for exp, statuses in status_by_exp.items():
        print_experiment_summary(
            exp, statuses=statuses, reasons=failure_reasons.get(exp)
        )
    if all_backed_up:
        click.echo(click.style("All experiments backed up successfully!", fg="green"))
    click.echo("")
    click.echo("----------------")


def print_experiment_summary(exp, *, statuses, reasons):
    """
    Prints the backup status line for a single experiment
    """
    if len(statuses) > 0:
        if statuses[0]:
            click.echo(click.style("✓", fg="green"), nl=False)
        else:
            click.echo(click.style("✗", fg="red"), nl=False)
    if len(statuses) > 1:
        if statuses[1]:
            click.echo(click.style("✓", fg="green"), nl=False)
        else:
            click.echo(click.style("✗", fg="red"), nl=False)
    click.echo(click.style(f" {exp} "), nl=False)
    if reasons is not None:
        click.echo(click.style(",".join(reasons), fg="red"), nl=True)
    else:
        click.echo()
