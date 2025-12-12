import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Optional

import click

from nomadic.util import minknow
from nomadic.util.dirs import produce_dir
from nomadic.util.settings import load_settings, settings_filepath
from nomadic.util.ssh import remote_dir_exists, ensure_remote_dir, is_ssh_target
from nomadic.util.workspace import Workspace


def selective_rsync(
    source_dir: Path,
    target_dir: Path | str,
    exclusions: Optional[list] = None,
    recursive: bool = False,
    delete: bool = False,
    checksum: bool = False,
    verbose: bool = False,
    progressbar: bool = False,
    update: bool = False,
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
    # --modify-window is added to help working with fat32
    rsync_components: list[str | Path] = ["rsync", "-rt", "--modify-window=2"]

    # Add variables
    if progressbar:
        rsync_components.append("--info=progress2")
    if verbose:
        rsync_components.append("-v")
    if delete:
        rsync_components.append("--delete")
    # Only update newer files on target (don't overwrite) when requested by caller
    if update:
        rsync_components.append("--update")
    if not recursive:
        rsync_components.extend(["--exclude", "*/"])
    if checksum:
        rsync_components.append("--checksum")
    if exclusions is not None:
        for exclusion in exclusions:
            rsync_components.extend(["--exclude", exclusion])

    # Complete the list:
    rsync_components.extend([source_dir, target_dir])

    # Format the rsync command properly for bash to run it
    rsync_command = [
        f"{f.resolve()}/" if isinstance(f, Path) else f for f in rsync_components
    ]

    if verbose:
        click.echo(f"{' '.join(rsync_command)}")

    result = subprocess.run(rsync_command, text=True, check=True)
    if result.stdout:
        click.echo(f"stdout: {result.stdout}")
    if result.stderr:
        click.echo(f"stderr: {result.stderr}")


def backup_nomadic_workspace(
    target_dir: Path | str,
    workspace: Workspace,
    *,
    checksum: bool = False,
    verbose: bool = False,
):
    workspace_path = Path(workspace.path).resolve()
    workspace_name = workspace.get_name()
    click.echo(f"Backing up nomadic workspace ({workspace_name}) to {target_dir}")

    exclusions = [".incremental/", "intermediate/", ".work.log"]

    selective_rsync(
        source_dir=workspace_path,
        target_dir=target_dir,
        recursive=True,
        progressbar=True,
        exclusions=exclusions,
        checksum=checksum,
        verbose=verbose,
    )
    click.echo("Done.")


def share_nomadic_workspace(
    target_dir: Path | str,
    workspace: Workspace,
    *,
    checksum: bool = False,
    verbose: bool = False,
    update: bool = True,
):
    workspace_path = Path(workspace.path).resolve()
    workspace_name = workspace.get_name()
    click.echo(f"Sharing nomadic workspace ({workspace_name}) to {target_dir}")

    exclusions = [
        ".incremental/",
        "intermediate/",
        ".work.log",
        "barcodes/",
    ]

    selective_rsync(
        source_dir=workspace_path,
        target_dir=target_dir,
        recursive=True,
        progressbar=True,
        exclusions=exclusions,
        checksum=checksum,
        verbose=verbose,
        update=update,
    )
    click.echo("Done.")


def backup_minknow_data(
    target_base_dir: Path | str,
    minknow_base_dir: Path,
    workspace: Workspace,
    *,
    failure_reasons: dict[str, list[str]],
    checksum: bool = False,
    verbose: bool = False,
):
    click.echo(f"Backing up minknow data to {target_base_dir}")
    sync_minknow_data(
        target_base_dir=target_base_dir,
        minknow_base_dir=minknow_base_dir,
        workspace=workspace,
        failure_reasons=failure_reasons,
        checksum=checksum,
        verbose=verbose,
    )


def share_minknow_data(
    target_base_dir: Path | str,
    minknow_base_dir: Path,
    workspace: Workspace,
    failure_reasons: dict[str, list[str]],
    *,
    checksum: bool = False,
    verbose: bool = False,
    update: bool = True,
):
    click.echo(f"Sharing minknow data to {target_base_dir}")
    sync_minknow_data(
        target_base_dir=target_base_dir,
        minknow_base_dir=minknow_base_dir,
        workspace=workspace,
        failure_reasons=failure_reasons,
        exclusions=[
            "fastq_fail/",
            "fastq_pass/",
            "other_reports/",
            "pod5/",
            "sequencing_summary_*.txt",
        ],
        checksum=checksum,
        verbose=verbose,
        update=update,
    )


def sync_minknow_data(
    target_base_dir: Path | str,
    minknow_base_dir: Path,
    workspace: Workspace,
    failure_reasons: dict[str, list[str]],
    exclusions: Optional[list[str]] = None,
    *,
    checksum: bool = False,
    verbose: bool = False,
    update: bool = False,
):
    """
    Sync minknow data. This contains the core logic of finding the minknow data, but allows for different inclusions/exclusions
    """
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
        target_dir = get_minknow_target_dir(target_base_dir, exp)

        if not source_dir.exists():
            click.echo(f"   ERROR: {source_dir} does not exist, unable to sync...")
            failure_reasons[exp].append("no minknow data found")
            continue
        if not minknow.is_minknow_experiment_dir(source_dir):
            failure_reasons[exp].append("invalid minknow data")
            click.echo(
                f"   ERROR: {source_dir} does not look like a valid minknow experiment directory, unable to sync..."
            )
            continue

        # Ensure target dir exists
        if isinstance(target_dir, Path):
            if not target_dir.exists():
                produce_dir(target_dir)
        elif is_ssh_target(target_dir):
            ok, msg = ensure_remote_dir(target_dir, verbose=verbose)
            if not ok:
                raise click.ClickException(
                    f"Unable to create/verify remote directory '{target_dir}': {msg}"
                )

        if exclusions is None:
            exclusions = []

        click.echo(f"   {source_dir} to {target_dir}")
        selective_rsync(
            source_dir=source_dir,
            target_dir=target_dir,
            recursive=True,
            progressbar=True,
            exclusions=exclusions,
            checksum=checksum,
            verbose=verbose,
            update=update,
        )


def sync_status(
    backup_dir: Path | str,
    workspace: Workspace,
    include_minknow: bool,
    verbose: bool,
):
    """
    Calculate sync status for all experiments in the workspace.

    This uses a very simple heuristic of checking if the result dir is there and relies
    on rsync working correctly. It is just a sanity check.
    """
    all_backed_up = True
    status_by_exp = defaultdict(list)
    for exp in workspace.get_experiment_names():
        # If backup_dir is a local Path, we can check existence. If it's a remote target (str) we can't determine status here.
        exp_dir = get_nomadic_target_dir(backup_dir, exp)
        if isinstance(exp_dir, Path):
            nomadic_backed_up = exp_dir.exists()
        elif is_ssh_target(exp_dir):
            nomadic_backed_up, msg = remote_dir_exists(exp_dir, verbose=False)
            if msg and verbose:
                click.echo(
                    f"Warning: unable to check exp backup status for experiment {exp}: {msg}"
                )
        else:
            nomadic_backed_up = None
        status_by_exp[exp].append(nomadic_backed_up)
        if nomadic_backed_up is False:
            all_backed_up = False
        if include_minknow:
            minknow_target = get_minknow_target_dir(backup_dir, exp)
            if isinstance(minknow_target, Path):
                minknow_backed_up = minknow_target.exists()
            elif is_ssh_target(minknow_target):
                minknow_backed_up, msg = remote_dir_exists(
                    minknow_target, verbose=verbose
                )
                if msg and verbose:
                    click.echo(
                        f"Warning: unable to check minknow backup status for experiment {exp}: {msg}"
                    )
            else:
                minknow_backed_up = None
            status_by_exp[exp].append(minknow_backed_up)
            if minknow_backed_up is False:
                all_backed_up = False
    return all_backed_up, status_by_exp


def print_sync_summary(all_synced_up, status_by_exp, failure_reasons, include_minknow):
    """
    Prints a summary of the sync, showing which of the experiments have been synces successfully and
    if not, what the failure reasons are.
    """
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
        click.echo(click.style("All experiments synced successfully!", fg="green"))
    click.echo("")
    click.echo("----------------")


def print_rsync_experiment_summary(exp, *, statuses, reasons):
    """
    Prints the rsync status line for a single experiment
    """
    if len(statuses) > 0:
        if statuses[0] is True:
            click.echo(click.style("✓", fg="green"), nl=False)
        elif statuses[0] is False:
            click.echo(click.style("✗", fg="red"), nl=False)
        else:
            # Unknown (remote targets)
            click.echo(click.style("?", fg="yellow"), nl=False)
    if len(statuses) > 1:
        if statuses[1] is True:
            click.echo(click.style("✓", fg="green"), nl=False)
        elif statuses[1] is False:
            click.echo(click.style("✗", fg="red"), nl=False)
        else:
            click.echo(click.style("?", fg="yellow"), nl=False)
    click.echo(click.style(f" {exp} "), nl=False)
    if reasons is not None:
        click.echo(click.style(",".join(reasons), fg="red"), nl=True)
    else:
        click.echo()


def get_minknow_target_dir(dir: Path | str, exp: str):
    if isinstance(dir, Path):
        return dir / "minknow" / exp
    # assume string remote target like user@host:/path
    return f"{dir.rstrip('/')}/minknow/{exp}"


def get_nomadic_target_dir(dir: Path | str, exp: str):
    if isinstance(dir, Path):
        return dir / "results" / exp
    return f"{dir.rstrip('/')}/results/{exp}"
