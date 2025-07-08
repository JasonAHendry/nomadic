import os

import click
from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.workspace import (
    check_if_workspace,
    Workspace,
    looks_like_a_bed_filepath,
)
from nomadic.util.minknow import default_fastq_dir
from nomadic.util.config import load_config, default_config_path


def set_workspace(ctx, param, workspace_path):
    """Set the workspace path in the context object so we can use it later."""
    ctx.ensure_object(dict)
    ctx.obj["workspace_path"] = workspace_path
    return workspace_path


def load_defaults_from_config(ctx, param, value):
    """Load configuration from the default config file if it exists."""
    config_path = os.path.join(ctx.obj.get("workspace_path", "./"), default_config_path)
    if os.path.isfile(config_path):
        defaults = load_config(config_path).get("defaults", None)
        if defaults is not None:
            click.echo(f"Loaded defaults from {config_path}")
            ctx.default_map = defaults
            ctx.show_default = True


@click.command(
    short_help="Run analysis in real-time.",
)
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
    callback=set_workspace,
)
@click.option(
    callback=load_defaults_from_config,
    expose_value=False,
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
                message=f"Metadata CSV file not found at {metadata_csv}. Did you create your metadata file in `{workspace.get_metadata_dir()}` and does the name match `{experiment_name}`?",
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
