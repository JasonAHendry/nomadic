import subprocess
from collections import defaultdict
from pathlib import Path

import click

from nomadic.util import minknow
from nomadic.util.dirs import produce_dir
from nomadic.util.settings import load_settings, settings_filepath
from nomadic.util.workspace import Workspace


def selective_rsync(
    source_dir: Path,
    target_dir: Path,
    exclusions: list = None,
    recursive: bool = False,
    delete: bool = False,
    checksum: bool = False,
    verbose: bool = False,
    progressbar: bool = False,
):
    """Copies contents of a folder to a new location.

    Args:
        source_dir(Path): The path to the source folder
        target_dir(Path): The path to the target folder
        exclusions(list): A list of file patterns to exclude
        recursive(bool): Copy top-level files or entire directory
        delete(bool): Delete files in target that are not in source
        checksum(bool): Use checksum to determine if files have changed
        verbose(bool): Whether to output verbose rsync information
        progressbar(bool): Whether to show a progress bar
    """
    # Base command with recursive (r) and timestamp (t) options
    # r is needed to select all entries in the folder even if the rsync is not recursive
    # a folder exclusion is then added
    rsync_components = ["rsync", "-rt"]

    # Add variables
    if progressbar:
        rsync_components.append("--info=progress2")
    if verbose:
        rsync_components.append("-v")
    if delete:
        rsync_components.append("--delete")
    if not recursive:
        rsync_components.extend(["--exclude", "*/"])
    if checksum:
        rsync_components.append("--checksum")
    if exclusions is not None:
        for exclusion in exclusions:
            rsync_components.extend(["--exclude", exclusion])

    # Complete the list:
    rsync_components.extend([source_dir, target_dir])

    if verbose:
        rsync_feedback = [
            f"{f.name}" if isinstance(f, Path) else f for f in rsync_components
        ]
        click.echo(f"{' '.join(rsync_feedback)}")

    # Format the rsync command properly for bash to run it
    rsync_command = [
        f"{f.resolve()}/" if isinstance(f, Path) else f for f in rsync_components
    ]
    result = subprocess.run(rsync_command, text=True, check=True)
    if result.stdout:
        click.echo(f"stdout: {result.stdout}")
    if result.stderr:
        click.echo(f"stderr: {result.stderr}")


def copy_nomadic_workspace(
    target_dir: Path,
    workspace: Workspace,
    additional_exclusions: list[str] = None,
):
    workspace_path = Path(workspace.path).resolve()
    workspace_name = workspace.get_name()
    click.echo(f"Copying nomadic workspace ({workspace_name}) to {target_dir}")

    exclusions = ["**/.incremental/", "**/intermediate", ".work.log"]

    if additional_exclusions is not None:
        exclusions = exclusions + additional_exclusions

    selective_rsync(
        source_dir=workspace_path,
        target_dir=target_dir,
        recursive=True,
        progressbar=True,
        exclusions=exclusions,
    )
    click.echo("Done.")


def copy_minknow_data(
    target_dir: Path,
    minknow_base_dir: Path,
    workspace: Workspace,
    failure_reasons: dict[str, list[str]],
    exclusions: list[str] = None,
):
    click.echo(f"Copying minknow data to {target_dir}")
    workspace_path = Path(workspace.path).resolve()

    for i, exp in enumerate(workspace.get_experiment_names()):
        click.echo(f"{exp} ({i + 1}/{len(workspace.get_experiment_names())})")
        json_file = settings_filepath(workspace_path, exp)
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
        target_dir = get_minknow_target_dir(target_dir, exp)

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

        if exclusions is None:
            exclusions = []

        click.echo(f"   {source_dir} to {target_dir}")
        selective_rsync(
            source_dir=source_dir,
            target_dir=target_dir,
            recursive=True,
            progressbar=True,
            exclusions=exclusions,
        )


def rsync_status(backup_dir: Path, workspace: Workspace, include_minknow: bool):
    """
    Calculate rsync status for all experiments in the workspace.

    This uses a very simple heuristic of checking if the result dir is there and relies
    on rsync working correctly. It is just a sanity check.
    """
    all_backed_up = True
    status_by_exp = defaultdict(list)
    for exp in workspace.get_experiment_names():
        exp_dir = get_nomadic_target_dir(backup_dir, exp)
        nomadic_backed_up = exp_dir.exists()
        status_by_exp[exp].append(nomadic_backed_up)
        if not nomadic_backed_up:
            all_backed_up = False
        if include_minknow:
            minknow_backed_up = get_minknow_target_dir(backup_dir, exp).exists()
            status_by_exp[exp].append(minknow_backed_up)
            if not minknow_backed_up:
                all_backed_up = False
    return all_backed_up, status_by_exp


def print_rsync_summary(all_synced_up, status_by_exp, failure_reasons, include_minknow):
    click.echo("")
    click.echo("--- Summary ---")
    click.echo("")
    n_nomadic, n_minknow = 0, 0
    for exp, statuses in status_by_exp.items():
        print_rsync_experiment_summary(
            exp, statuses=statuses, reasons=failure_reasons.get(exp)
        )
        if len(statuses) > 0 and statuses[0]:
            n_nomadic += 1
        if len(statuses) > 1 and statuses[1]:
            n_minknow += 1
    click.echo(f"nomadic: {n_nomadic}/{len(status_by_exp)}")
    if include_minknow:
        click.echo(f"minknow: {n_minknow}/{len(status_by_exp)}")
    if all_synced_up:
        click.echo(
            click.style("All experiments synchronised successfully!", fg="green")
        )
    click.echo("")
    click.echo("----------------")


def print_rsync_experiment_summary(exp, *, statuses, reasons):
    """
    Prints the rsync status line for a single experiment
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


def get_minknow_target_dir(dir, exp):
    return dir / "minknow" / exp


def get_nomadic_target_dir(dir, exp):
    return dir / "results" / exp
