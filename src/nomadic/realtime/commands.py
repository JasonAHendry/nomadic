import click


@click.command(short_help="Run nomadic in real-time.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'."
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
def realtime(expt_name, fastq_dir, metadata_csv, region_bed):
    from .main import main

    main(expt_name, fastq_dir, metadata_csv, region_bed)

