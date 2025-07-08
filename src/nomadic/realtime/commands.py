import os

import click
from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.workspace import (
    check_if_workspace,
    Workspace,
    looks_like_a_bed_filepath,
)
from nomadic.util.minknow import default_fastq_dir


@click.command(short_help="Run analysis in real-time.")
@click.argument(
    "experiment_name",
    type=str,
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the the workspace where all the files will be stored",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    show_default="/var/lib/minknow/data/<experiment_name>/fastq_pass",
    help="Path to `fastq_pass` directory produced by MinKNOW or Guppy.",
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
    required=True,
    help="Path to BED file specifying genomic regions of interest or name of panel, e.g. 'nomads8' or 'nomadsMVP'.",
)
@click.option(
    "-r",
    "--reference_name",
    type=click.Choice(REFERENCE_COLLECTION),
    default="Pf3D7",
    show_default=True,
    help="Choose a reference genome to be used in real-time analysis.",
)
@click.option(
    "-c",
    "--call",
    is_flag=True,
    default=False,
    help="Perform preliminary variant calling of biallelic SNPs in real-time.",
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
    workspace_path,
    fastq_dir,
    metadata_csv,
    region_bed,
    reference_name,
    call,
    verbose,
):
    """
    Analyse data being produced by MinKNOW while sequencing is ongoing
    """
    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Workspace '{workspace_path}' does not exist or is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
        )

    workspace = Workspace(workspace_path)

    if not metadata_csv:
        metadata_csv = workspace.get_meta_data_csv(experiment_name)
        if not os.path.isfile(metadata_csv):
            raise click.BadParameter(
                message=f"Metadata CSV file not found at {metadata_csv}. Does it match the experment name?",
            )

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

    if fastq_dir is None:
        fastq_dir = default_fastq_dir(experiment_name)
        if not os.path.isdir(fastq_dir):
            raise click.BadParameter(
                message=f"FASTQ directory not found at {fastq_dir}. Does {experiment_name} match the minknow experiment name and is minknow already running long enough?",
            )
    from .main import main

    main(
        experiment_name,
        workspace_path,
        fastq_dir,
        metadata_csv,
        region_bed,
        reference_name,
        call,
        verbose,
    )
