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

    # delete only works if recursive is True
    if recursive:
        if delete:
            rsync_components.append("--delete")
    else:
        rsync_components.extend(["--exclude", "*/"])

    if checksum:
        rsync_components.append("--checksum")

    # Add in specific exclusions if given
    if exclusions:
        for exclusion in exclusions:
            rsync_components.extend(["--exclude", exclusion])

    # Complete the list:
    rsync_components.extend([source_dir, target_dir])

    # Give user feedback on the rsync command being run
    rsync_feedback = [
        f"{f.name}" if isinstance(f, Path) else f for f in rsync_components
    ]
    # TODO: Ensure logging is working properly
    log.debug(f"{' '.join(rsync_feedback)}")

    try:
        # Format the rsync command properly for bash to run it
        rsync_command = [
            f"{f.resolve()}/" if isinstance(f, Path) else f for f in rsync_components
        ]
        result = subprocess.run(rsync_command, text=True, check=True)
        if result.stdout:
            log.debug(f"stdout: {result.stdout}")
        if result.stderr:
            log.warning(f"stderr: {result.stderr}")
    except Exception as e:
        log.error(f"An unexpected error occurred: {e}")
