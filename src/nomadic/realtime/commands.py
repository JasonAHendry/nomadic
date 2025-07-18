import click
from nomadic.download.references import REFERENCE_COLLECTION


@click.command(short_help="Run analysis in real-time.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'.",
)
@click.option(
    "-f",
    "--fastq_dir",
    type=str,
    required=True,
    help="Path to `fastq_pass` directory produced by MinKNOW or Guppy.",
)
@click.option(
    "-m",
    "--metadata_csv",
    type=str,
    required=True,
    help="Path to metadata CSV file containing barcode and sample information.",
)
@click.option(
    "-b",
    "--region_bed",
    type=str,
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
    expt_name, fastq_dir, metadata_csv, region_bed, reference_name, call, verbose
):
    """
    Analyse data being produced by MinKNOW while sequencing is ongoing
    """
    from .main import main

    main(expt_name, fastq_dir, metadata_csv, region_bed, reference_name, call, verbose)
