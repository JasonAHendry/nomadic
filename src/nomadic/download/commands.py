import click
from .references import REFERENCE_COLLECTION


@click.command(short_help="Download reference genomes.")
@click.option(
    "-r",
    "--reference_name",
    type=click.Choice(REFERENCE_COLLECTION),
    default=None,
    help="Choose a reference genome to download by name.",
)
@click.option("--all", is_flag=True, help="Download all reference genomes available.")
def download(reference_name: str, all: bool = False) -> None:
    """
    Download target reference genome by specifying a `reference_name`, or download
    all available genomes with the `--all` flag

    """

    from .main import main

    main(reference_name, all)
