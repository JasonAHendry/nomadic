import glob
import warnings
from pathlib import Path
from typing import Optional, Tuple

import click


def resolve_fastq_dir(fastq_dir_glob: str) -> Optional[str]:
    """
    Return the default FASTQ directory for MinKNOW experiments.
    This is typically where Minknow/Guppy writes the FASTQ files.

    """
    fastq_pass_dirs = sorted(glob.glob(fastq_dir_glob))
    if len(fastq_pass_dirs) > 1:
        warnings.warn(
            f"Found {len(fastq_pass_dirs)} 'fastq_pass' directories,"
            " This suggests more than one experiment with this name has"
            " been run from MinKNOW. Will proceed with most recent."
        )
    if len(fastq_pass_dirs) < 1:
        warnings.warn(
            f"Found no 'fastq_pass' directories in '{fastq_dir_glob}',"
            " If you just started minknow, this can mean the files have not been created yet."
            " If minknow is already running for more than 10 minutes, check if your experiment"
            " name matches minknow's experiment name."
        )
        return None

    return fastq_pass_dirs[-1]


def is_fastq_dir(path: Path) -> bool:
    return path.is_dir() and path.name == "fastq_pass"


def create_fastq_dir_glob(minknow_dir: Path, experiment_name: str) -> str:
    return str(minknow_dir / experiment_name / "*" / "*" / "fastq_pass")


def is_minknow_dir(path: Path) -> bool:
    expected_folders = {"persistence", "reads", "queued_reads", "intermediates"}
    if any(d.name in expected_folders for d in path.glob("*") if d.is_dir()):
        return True
    else:
        return False


def extract_minknow_dir_from_fastq_dir(fastq_dir: Path) -> Path:
    minknow_dir = fastq_dir.parents[3]
    if is_minknow_dir(minknow_dir):
        return minknow_dir
    else:
        click.BadParameter(
            "Unable to determine the parent MinKNOW directory from the provided fastq directory."
        )


def resolve_minknow_fastq_dirs(
    path: Path, experiment_name: str
) -> Optional[Tuple[Path, Path]]:
    """
    This function looks to see if the supplied path resembles a minknow data folder or a
    specific fastq_pass folder from a specific experiment
    """
    if is_minknow_dir(path):
        # should be base path of minknow data, build fastq glob with experiment name.
        fastq_dir_glob = create_fastq_dir_glob(path, experiment_name)
        fastq_dir = resolve_fastq_dir(fastq_dir_glob)
        minknow_dir = path
    elif is_fastq_dir(path):
        fastq_dir = path
        # determine the minknow dir from the fastq dir
        minknow_dir = extract_minknow_dir_from_fastq_dir(path)
    else:
        raise click.BadParameter(
            message=f"{path} does not look like a valid MinKNOW output directory or a valid `fastq_pass` directory.",
        )
    return minknow_dir, fastq_dir
