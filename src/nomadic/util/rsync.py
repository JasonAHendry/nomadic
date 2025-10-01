import logging
import subprocess
from pathlib import Path

script_dir = Path(__file__).parent.resolve()
log = logging.getLogger(script_dir.stem)


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
    # Base command with compress (z), recursive (r)) and timestamp (t) options
    # r is needed to select all entries in the folder even if the rsync is not recursive
    # a folder exclusion is then added
    rsync_components = ["rsync", "-zrt"]

    # Progressbar or verbose (can't have both)
    if progressbar:
        rsync_components.append("--info=progress2")
    elif verbose:
        rsync_components.append("-v")

    # Add in other supplied details
    if delete:
        rsync_components.append("--delete")
    if not recursive:
        rsync_components.extend(["--exclude", "*/"])
    if checksum:
        rsync_components.append("--checksum")
    if exclusions:
        for exclusion in exclusions:
            rsync_components.extend(["--exclude", exclusion])

    # Complete the list:
    rsync_components.extend([source_dir, target_dir])

    # Give user feedback on the rsync command being run
    rsync_feedback = [
        f"{f.name}" if isinstance(f, Path) else f for f in rsync_components
    ]
    print(f"{' '.join(rsync_feedback)}")

    # Format the rsync command properly for bash to run it
    rsync_command = [
        f"{f.resolve()}/" if isinstance(f, Path) else f for f in rsync_components
    ]
    result = subprocess.run(rsync_command, text=True, check=True)
    if result.stdout:
        print(f"stdout: {result.stdout}")
    if result.stderr:
        print(f"stderr: {result.stderr}")
