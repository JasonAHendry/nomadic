import glob
import platform
import warnings
from pathlib import Path
from typing import Optional, Tuple


class MinknowPathError(Exception):
    pass


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
    return path.is_dir() and any(path.glob("*.fastq*"))


def create_fastq_dir_glob(minknow_dir: Path) -> str:
    return str(minknow_dir / "*" / "*" / "fastq_pass")


def is_minknow_base_dir(path: Path) -> bool:
    expected_folders = {"persistence", "reads", "queued_reads", "intermediates"}
    has_expected_folders = any(
        d.name in expected_folders for d in path.glob("*") if d.is_dir()
    )
    if has_expected_folders:
        return True

    # Sometimes the expected folders are missing if the experiments where copied out
    # Check if there are any minknow experiments in the folder
    return any(is_minknow_experiment_dir(d) for d in path.glob("*") if d.is_dir())


def is_minknow_experiment_dir(path: Path) -> bool:
    expected_folders = {"pod5", "fastq_pass", "fastq_fail"}
    return any(d.name in expected_folders for d in path.glob("*/*/*") if d.is_dir())


def resolve_minknow_fastq_dirs(
    minknow_path: Path, experiment_name: str
) -> Tuple[Path, str]:
    """
    This function looks to see if the supplied path resembles a minknow data folder or a
    specific fastq_pass folder from a specific experiment
    """
    if not minknow_path.exists():
        raise MinknowPathError(
            f"{minknow_path} does not exist.",
        )

    if is_minknow_base_dir(minknow_path):
        minknow_dir = minknow_path / experiment_name
    elif is_minknow_experiment_dir(minknow_path):
        minknow_dir = minknow_path
    else:
        raise MinknowPathError(
            f"{minknow_path} does not look like a valid MinKNOW output directory. "
            f"Please ensure it either points to a minknow experiment folder containing a fastq_pass folder, "
            "or a folder containing minknow experiments.",
        )

    fastq_dir = create_fastq_dir_glob(minknow_dir)

    return minknow_dir, fastq_dir


def default_data_dir() -> Path:
    """
    Return the default data directory for MinKNOW experiments.
    This is typically where Minknow/Guppy writes the FASTQ files.

    see: https://nanoporetech.com/support/software/MinKNOW/post-run-options/what-folder-are-my-reads-in
    """
    system = platform.system()
    if system == "Darwin":  # MacOS
        return Path("/Library/MinKNOW/data")
    elif system == "Windows":  # Windows
        # We don't really support windows, but for completeness sake
        return Path("C:\\data\\")
    else:  # Unix/Linux
        standard_path = Path("/var/lib/minknow/data")
        integrated_devices_path = Path("/data")  # e.g. gridION

        if (
            not standard_path.is_dir()
            and integrated_devices_path.is_dir()
            and is_minknow_base_dir(integrated_devices_path)
        ):
            # we are on an integrated device
            return integrated_devices_path

        return standard_path
