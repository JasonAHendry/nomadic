import os
from importlib.resources import files
from dataclasses import dataclass

import click

from nomadic.util.config import default_config_path, write_config
from nomadic.util.workspace import Workspace


# --------------------------------------------------------------------------------
# Workspace definitions for different organisms
#
# To add an organism:
# - Define an Organism inside _organism = [...]
# - Ensure the reference genome is available: src/nomadic/download/references.py
# - Add BED(s) to src/nomadic/start/data/beds/<organism_name>/*.bed
# --------------------------------------------------------------------------------


@dataclass
class Organism:
    name: str
    reference: str
    default_bed: str
    caller: str


_organisms = [
    Organism("pfalciparum", "Pf3D7", "nomadsMVP", "delve"),
    Organism("agambiae", "AgPEST", "nomadsIR", "bcftools"),
]
ORGANISM_COLLECTION = {organism.name: organism for organism in _organisms}


# --------------------------------------------------------------------------------
# Workspace creation (downloading / copying files)
#
# --------------------------------------------------------------------------------


def setup_organism(
    workspace: Workspace,
    organism: Organism,
) -> None:
    """
    Setup a workspace for a specific organism
    """
    from nomadic.download.main import main as download_reference

    download_reference(organism.reference)
    copy_bed_files(workspace, organism_name=organism.name)
    copy_example_metadata(workspace)

    click.echo(f"Setting reference genome: {organism.reference}")
    click.echo(f"Setting default BED file: {organism.default_bed}")
    click.echo(f"Setting default variant caller: {organism.caller}")

    defaults = {
        "reference_name": organism.reference,
        "region_bed": organism.default_bed,
        "caller": organism.caller,
    }
    write_config(
        config={"defaults": defaults},
        config_path=os.path.join(workspace.path, default_config_path),
    )


def copy_example_metadata(workspace):
    click.echo("Copying example metadata file.")
    example_metadata_file = files("nomadic.start").joinpath(
        "data", "0000-00-00_example.csv"
    )
    data = example_metadata_file.read_text()
    dest_path = os.path.join(workspace.get_metadata_dir(), "0000-00-00_example.csv")
    with open(dest_path, "w") as text_file:
        text_file.write(data)


def copy_bed_files(workspace: Workspace, *, organism_name):
    bed_files = files("nomadic.start").joinpath("data", "beds", organism_name).iterdir()
    click.echo("Copying amplicon BED files.")
    for bed_file in bed_files:
        data = bed_file.read_text()
        dest_path = os.path.join(workspace.get_beds_dir(), bed_file.name)
        with open(dest_path, "w") as text_file:
            text_file.write(data)


# --------------------------------------------------------------------------------
# Main command
#
#
# --------------------------------------------------------------------------------


@click.argument(
    "organism",
    type=click.Choice(ORGANISM_COLLECTION, case_sensitive=False),
    required=True,
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="nomadic",
    type=click.Path(exists=False),
    show_default=True,
    help="Path to workspace.",
)
@click.command(short_help="Start a workspace.")
def start(organism, workspace_path) -> None:
    """
    Get started with nomadic.

    This command will help you set up a new workspace for a specific organism.

    \b
    Currently supported organisms:
    - Plasmodium falciparum (pfalciparum)
    - Anopheles gambiae (agambiae)

    """

    click.echo(f"Workspace will be created at: {workspace_path}")
    workspace = Workspace.create_from_directory(workspace_path)

    if organism not in ORGANISM_COLLECTION:  # this should be handled by click.
        raise RuntimeError(
            f"Organism {organism} is not available. Choose from {', '.join(ORGANISM_COLLECTION)}."
        )

    setup_organism(workspace=workspace, organism=ORGANISM_COLLECTION[organism])

    click.echo(
        f"You can now enter your workspace with `cd {workspace_path}`"
        " and run `nomadic realtime <experiment_name>` to start real-time analysis."
    )
