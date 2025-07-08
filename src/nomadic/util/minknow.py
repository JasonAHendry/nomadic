def default_fastq_dir(experiment_name: str) -> str:
    """
    Return the default FASTQ directory for MinKNOW experiments.
    This is typically where Minknow/Guppy writes the FASTQ files.
    """
    return f"/var/lib/minknow/data/{experiment_name}/fastq_pass"
