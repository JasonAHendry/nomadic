import click
from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.workspace import check_if_workspace


@click.command(short_help="Run analysis in real-time.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'.",
)
@click.option(
    "-w",
    "--workspace",
    default="./",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the the workspace where all the files will be stored. If not given, users default workspace will be used, set by init",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Path to `fastq_pass` directory produced by MinKNOW or Guppy.",
)
@click.option(
    "-m",
    "--metadata_csv",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="Path to metadata CSV file containing barcode and sample information.",
)
@click.option(
    "-b",
    "--region_bed",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="Path to BED file specifying genomic regions of interest.",
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
    expt_name,
    workspace,
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
    from .main import main

    if not check_if_workspace(workspace):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"Workspace '{workspace}' does not exist or is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
        )

    main(
        expt_name,
        workspace,
        fastq_dir,
        metadata_csv,
        region_bed,
        reference_name,
        call,
        verbose,
    )
