import glob
import warnings
import os
from typing import Optional


def fastq_dir(fastq_dir_glob: str) -> Optional[str]:
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


def is_fastq_dir(path: str) -> bool:
    return (
        path.endswith("fastq_pass") or path.endswith("fastq_pass/")
    ) and os.path.isdir(path)


def fastq_dir_glob(data_dir: str, experiment_name: str) -> str:
    if not data_dir.endswith(experiment_name) and not data_dir.endswith(
        f"{experiment_name}/"
    ):
        data_dir = os.path.join(data_dir, experiment_name)
    return f"{data_dir}/*/*/fastq_pass"
