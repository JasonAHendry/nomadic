import glob
import warnings


def default_fastq_dir(experiment_name: str) -> str:
    """
    Return the default FASTQ directory for MinKNOW experiments.
    This is typically where Minknow/Guppy writes the FASTQ files.

    """
    fastq_pass_dirs = glob.glob(
        f"/var/lib/minknow/data/{experiment_name}/*/*/fastq_pass"
    )
    if len(fastq_pass_dirs) > 1:
        warnings.warn(
            f"Found {len(fastq_pass_dirs)} 'fastq_pass' directories,"
            " This suggests more than one experiment with this name has"
            " been run from MinKNOW. Will proceed with most recent."
        )

    return fastq_pass_dirs[-1]
