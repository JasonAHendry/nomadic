import os

import click

from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util import minknow
from nomadic.util.config import default_config_path, load_config
from nomadic.util.workspace import (
    Workspace,
    check_if_workspace,
    looks_like_a_bed_filepath,
)


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
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the workspace where all input/output files (beds, metadata, results) are stored. "
    "The workspace directory simplifies the use of nomadic in that many arguments don't need to be listed "
    "as they are predefined in the workspace config or can be loaded from the workspace",
    callback=set_workspace,
)
@click.option(
    callback=load_defaults_from_config,
    expose_value=False,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    show_default="<workspace>/results/<experiment_name>",
    help="Path to the output directory where results of this experiment will be stored. Usually the default of storing it in the workspace should be enough.",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    default="/var/lib/minknow/data",
    show_default=True,
    help="Path to `fastq_pass` output directory of minknow (e.g. `/var/lib/minknow/data/<experiment_name>/.../.../fastq_pass`), or to the general output directory of minknow (e.g. `/var/lib/minknow/data`)",
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
    "--resume",
    is_flag=True,
    default=False,
    help="Resume a previous experiment run if the output directory already exists. Only use if you want to force resuming an already started experiment. "
    "Not needed in interactive mode as this will be prompted",
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
    fastq_dir,
    metadata_csv,
    region_bed,
    reference_name,
    call,
    resume,
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

    if output is None:
        output = workspace.get_output_dir(experiment_name)

    if not metadata_csv:
        metadata_csv = workspace.get_metadata_csv(experiment_name)
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

    if not minknow.is_fastq_dir(fastq_dir):
        # should be base path of minknow data, build fastq glob with experiment name.
        fastq_dir = minknow.fastq_dir_glob(fastq_dir, experiment_name)

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

    if os.path.exists(output) and not resume:
        click.confirm(
            f"Output directory {output} already exists. Do you want to resume a previous experiment run? If starting a new experiment, please restart with a different experiment name or output directory.",
            abort=True,
        )

    from .main import main

    main(
        experiment_name,
        output,
        workspace_path,
        fastq_dir,
        metadata_csv,
        region_bed,
        reference_name,
        call,
        verbose,
    )
