import subprocess
from pathlib import Path

import click


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
