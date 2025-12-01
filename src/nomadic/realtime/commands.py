import os
from pathlib import Path
from shutil import rmtree

import click

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util import minknow
from nomadic.util.workspace import (
    Workspace,
    check_if_workspace,
    looks_like_a_bed_filepath,
)
from nomadic.util.exceptions import MetadataFormatError
from nomadic.util.cli import (
    complete_experiment_name,
    complete_bed_file,
    load_default_function_for,
)


@click.command(
    short_help="Run analysis in real-time.",
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
    callback=load_default_function_for("realtime"),
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
    "--call",
    is_flag=True,
    default=False,
    help="Perform preliminary variant calling of biallelic SNPs in real-time. (Deprecated, use --caller instead)",
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
    help="Resume a previous experiment run if the output directory already exists. Only use if you want to force resuming an already started experiment. "
    "Not needed in interactive mode as this will be prompted",
)
@click.option(
    "--dashboard/--no-dashboard",
    default=True,
    help="Whether to start the web dashboard to monitor the run.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Increase logging verbosity. Helpful for debugging.",
)
def realtime(
    experiment_name,
    output,
    workspace_path,
    minknow_dir,
    fastq_dir,
    metadata_csv,
    region_bed,
    reference_name,
    call,
    caller,
    overwrite,
    resume,
    dashboard,
    verbose,
):
    """
    Analyse data being produced by MinKNOW while sequencing is ongoing
    """
    workspace = get_workspace(workspace_path)
    output = get_output_path(experiment_name, output, workspace)
    metadata_csv = get_metadata_path(experiment_name, metadata_csv, workspace)
    region_bed = get_region_path(region_bed, workspace)
    minknow_dir, fastq_dir = get_minknow_fastq_dirs(
        experiment_name, minknow_dir, fastq_dir
    )

    validate_reference(reference_name)

    if os.path.exists(output):
        if overwrite:
            rmtree(output)
        elif not resume:
            choice = click.prompt(
                f"Output directory {output} already exists. Do you want to resume (y) a previous experiment run? If not, you can restart (r) the run, which will delete the existing output directory.",
                type=click.Choice(["y", "n", "r"], case_sensitive=False),
            )
            if choice == "n":
                raise click.Abort()
            elif choice == "r":
                click.confirm(
                    f"Are you sure you want to delete the existing output directory {output} and restart the experiment? This will delete all existing results in that directory.",
                    abort=True,
                )
                rmtree(output)

    if not caller and call:
        caller = "bcftools"

    from .main import main

    try:
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
            with_dashboard=dashboard,
            realtime=True,
        )
    except MetadataFormatError as e:
        raise click.BadParameter(
            param_hint="-m/--metadata_csv",
            message=str(e),
        ) from e


def get_minknow_fastq_dirs(experiment_name, minknow_dir, fastq_dir):
    if fastq_dir is None:
        return minknow.resolve_minknow_fastq_dirs(minknow_dir, experiment_name)
    else:
        # If fastq_dir is manually given, we assume there is no minknow dir
        return None, fastq_dir


def validate_reference(reference_name):
    if reference_name is None:
        raise click.BadParameter(
            param_hint="-r/--reference_name",
            message="Reference genome must be specified. Use -r/--reference_name to select a reference genome.",
        )
    elif reference_name not in REFERENCE_COLLECTION:
        raise click.BadParameter(
            param_hint="-r/--reference_name",
            message=f"Reference genome '{reference_name}' is not available. Available references: {', '.join(REFERENCE_COLLECTION.keys())}.",
        )


def get_region_path(region_bed, workspace):
    if not os.path.isfile(region_bed):
        if not looks_like_a_bed_filepath(region_bed):
            # Assume it's a panel name
            region_bed = workspace.get_bed_file(region_bed)
            if not os.path.isfile(region_bed):
                raise click.BadParameter(
                    message=f"Region BED file not found at {region_bed}. Possible panel names: {workspace.get_panel_names()}.",
                )
        else:
            raise click.BadParameter(
                message=f"Region BED file not found at {region_bed}.",
            )

    return region_bed


def get_metadata_path(experiment_name, metadata_csv, workspace):
    if not metadata_csv:
        metadata_csv = workspace.get_metadata_csv(experiment_name)
        if not os.path.isfile(metadata_csv):
            raise click.BadParameter(
                message=f"Metadata CSV file not found at {metadata_csv}. Did you create your metadata file in `{workspace.get_metadata_dir()}` and does the name match `{experiment_name}`?",
            )

    return metadata_csv


def get_output_path(experiment_name, output, workspace):
    if output is None:
        output = workspace.get_output_dir(experiment_name)
    return output


def get_workspace(workspace_path):
    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Workspace '{workspace_path}' does not exist or is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
        )

    workspace = Workspace(workspace_path)
    return workspace
