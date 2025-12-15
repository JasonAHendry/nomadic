import os
from pathlib import Path
from shutil import rmtree
from typing import Optional

import click

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util import minknow
from nomadic.util.cli import (
    BadParameterWithSource,
    complete_bed_file,
    complete_experiment_name,
    load_default_function_for,
    workspace_option,
)
from nomadic.util.exceptions import MetadataFormatError
from nomadic.util.workspace import (
    Workspace,
    looks_like_a_bed_filepath,
)


@click.command(
    short_help="Run analysis in real-time.",
)
@workspace_option(optional=True)
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
    "--metadata_path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="Path to metadata file (CSV or XLSX (Excel)) containing barcode and sample information.",
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
    experiment_name: str,
    output,
    workspace: Optional[Workspace],
    minknow_dir,
    fastq_dir,
    metadata_path,
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
    are_none = [
        value
        for value in [
            (output, "--output"),
            (metadata_path, "--metadata_path"),
            (region_bed, "--region_bed", reference_name, "--reference_name"),
        ]
        if value[0] is None
    ]
    if len(are_none) > 0 and workspace is None:
        # if any of those command line options is not set explicitly, we need a workspace to resolve them.
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Current directory is not a workspace. Please use nomadic start to create a new workspace, navigate to your workspace, or set {are_none[0][1]} explicit.",
        )

    output = determine_output_path(experiment_name, output, workspace)
    metadata_path = find_metadata_file(experiment_name, metadata_path, workspace)
    region_bed = find_region_file(region_bed, workspace)
    minknow_dir, fastq_dir = find_minknow_fastq_dirs(
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
            fastq_dir,
            minknow_dir,
            metadata_path,
            region_bed,
            reference_name,
            caller,
            verbose,
            with_dashboard=dashboard,
            realtime=True,
        )
    except MetadataFormatError as e:
        raise click.BadParameter(
            param_hint="-m/--metadata_path",
            message=str(e),
        ) from e


def find_minknow_fastq_dirs(
    experiment_name: str, minknow_dir: Path, fastq_dir
) -> tuple[Optional[Path], str]:
    """Finds the minknow and fastqdir folders to use"""
    if fastq_dir is None:
        return minknow.resolve_minknow_fastq_dirs(minknow_dir, experiment_name)
    else:
        # If fastq_dir is manually given, we assume there is no minknow dir
        return None, fastq_dir


def validate_reference(reference_name: str):
    if reference_name is None:
        raise BadParameterWithSource(
            param_hint="-r/--reference_name",
            message="Reference genome must be specified. Use -r/--reference_name to select a reference genome.",
        )
    elif reference_name not in REFERENCE_COLLECTION:
        raise BadParameterWithSource(
            param_hint="-r/--reference_name",
            message=f"Reference genome '{reference_name}' is not available. Available references: {', '.join(REFERENCE_COLLECTION.keys())}.",
        )


def find_region_file(region_bed: str, workspace: Workspace) -> str:
    """Tries to find the region bed file to use and signals to the user the right feedback if it can't be found."""
    if not os.path.isfile(region_bed):
        if not looks_like_a_bed_filepath(region_bed):
            # Assume it's a panel name
            region_bed = workspace.get_bed_file(region_bed)
            if not os.path.isfile(region_bed):
                raise BadParameterWithSource(
                    param_hint="-b/--region_bed",
                    message=f"Region BED file not found at {region_bed}. Possible panel names: {workspace.get_panel_names()}.",
                )
        else:
            raise BadParameterWithSource(
                param_hint="-b/--region_bed",
                message=f"Region BED file not found at {region_bed}.",
            )

    return region_bed


def find_metadata_file(
    experiment_name: str, metadata_path: Optional[str], workspace: Workspace
) -> str:
    """Finds the metadata file to use and gives feedback to the user if it can't be found"""
    if metadata_path is None:
        files = [
            workspace.get_metadata_csv(experiment_name),
            workspace.get_metadata_xlsx(experiment_name),
        ]

        shared_folder = workspace.get_shared_folder()
        if shared_folder is not None:
            click.echo(f"Found shared folder ({shared_folder})...")
            shared_workspace = Workspace(shared_folder)
            # Currently not checking if it actually is a workspace, to not require some of the folders that are not needed here
            files.extend(
                [
                    shared_workspace.get_metadata_csv(experiment_name),
                    shared_workspace.get_metadata_xlsx(experiment_name),
                ]
            )

        for file in files:
            if os.path.isfile(file):
                metadata_path = file
                break

    if metadata_path is None or not os.path.isfile(metadata_path):
        msg = f"Metadata file not found. Did you create your metadata file in `{workspace.get_metadata_dir()}`"
        if shared_folder:
            msg += f", or in `{shared_folder}`"
        msg += f", and does the name match `{experiment_name}`?"
        raise click.BadParameter(message=msg)
    return metadata_path


def determine_output_path(
    experiment_name: str, output: Optional[str], workspace: Workspace
):
    if output is None:
        output = workspace.get_output_dir(experiment_name)

    return output
