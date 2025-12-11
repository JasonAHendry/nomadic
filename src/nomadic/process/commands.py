import os
from pathlib import Path
from shutil import rmtree

import click

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.realtime.commands import (
    find_metadata_file,
    find_minknow_fastq_dirs,
    determine_output_path,
    find_region_file,
    find_workspace,
)
from nomadic.util.cli import (
    complete_bed_file,
    complete_experiment_name,
    load_default_function_for,
)


@click.command(
    short_help="(Re)process data from a completed run.",
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the workspace where all input/output files (beds, metadata, results) are stored. "
    "The workspace directory simplifies the use of nomadic in that many arguments don't need to be listed "
    "as they are predefined in the workspace config or can be loaded from the workspace",
)
@click.argument(
    "experiment_name",
    type=str,
    shell_complete=complete_experiment_name,
)
@click.option(
    callback=load_default_function_for("process"),
    expose_value=False,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=True, file_okay=False),
    show_default="<workspace>/results/<experiment_name>",
    help="Path to the output directory where results of this experiment will be stored. Usually the default of storing it in the workspace should be enough.",
)
@click.option(
    "-k",
    "--minknow_dir",
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default="/var/lib/minknow/data",
    show_default=True,
    help="Path to the minknow output directory. Can be either the base directory, e.g. /var/lib/minknow/data, or the directory of the experiment, e.g. /var/lib/minknow/data/<experiment_name>.",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=click.Path(file_okay=False, dir_okay=True),
    help="Path or glob to the fastq files. This should only be used when the full minknow dir can not be provided, as some features likes backup will not work. Prefer using --minknow_dir. If --fastq_dir is provided, --minknow_dir is ignored.",
)
@click.option(
    "-m",
    "--metadata_csv",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="Path to metadata CSV file containing barcode and sample information.",
    show_default="<workspace>/metadata/<experiment_name>.csv",
)
@click.option(
    "-b",
    "--region_bed",
    type=click.Path(),
    required=False,
    help="Path to BED file specifying genomic regions of interest or name of panel, e.g. 'nomads8' or 'nomadsMVP'.",
    shell_complete=complete_bed_file,
)
@click.option(
    "-r",
    "--reference_name",
    type=click.Choice(REFERENCE_COLLECTION),
    help="Choose a reference genome to be used in real-time analysis.",
)
@click.option(
    "-c",
    "--caller",
    help="Call biallelic SNPs in real-time with the indicated variant caller. If this flag is omitted, no variant calling is performed.",
    default=None,
    type=click.Choice(["bcftools", "delve"]),
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output directory if it exists.",
)
@click.option(
    "--resume",
    is_flag=True,
    default=False,
    help="Resume processing a previous experiment if the output directory already exists. This is necessary to pick of processing of an experiment that was aborted.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Increase logging verbosity. Helpful for debugging.",
)
def process(
    experiment_name,
    output,
    workspace_path,
    minknow_dir,
    fastq_dir,
    metadata_csv,
    region_bed,
    reference_name,
    caller,
    overwrite,
    resume,
    verbose,
):
    """
    (Re)Process data that was produced by MinKNOW
    """
    workspace = find_workspace(workspace_path)
    output = determine_output_path(experiment_name, output, workspace)
    metadata_csv = find_metadata_file(experiment_name, metadata_csv, workspace)
    region_bed = find_region_file(region_bed, workspace)
    minknow_dir, fastq_dir = find_minknow_fastq_dirs(
        experiment_name, minknow_dir, fastq_dir
    )

    if os.path.exists(output):
        if overwrite:
            rmtree(output)
        elif not resume:
            raise click.BadParameter(
                message=f"Output directory '{output}' already exists. Please choose a different output directory with -o/--output, or delete the existing directory if you want to overwrite it.",
            )

    from ..realtime.main import main

    # Calling the main function from realtime with realtime=False for now as they are very similar
    # If they diverge more in the future, we should split them up
    main(
        experiment_name,
        output,
        workspace_path,
        fastq_dir,
        minknow_dir,
        metadata_csv,
        region_bed,
        reference_name,
        caller,
        verbose,
        with_dashboard=False,
        realtime=False,
    )
