import enum
import os
from contextlib import chdir
from importlib.resources import files

import click

from nomadic.download.main import main as download_reference
from nomadic.download.references import REFERENCE_COLLECTION
from nomadic.util.config import default_config_path, write_config
from nomadic.util.workspace import Workspace, init_workspace


class Organism(enum.Enum):
    pfalciparum = enum.auto()


@click.argument(
    "organism",
    type=click.Choice(Organism),
    required=True,
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="nomadic",
    type=click.Path(exists=False),
    help="Path of the the workspace where all the files will be stored",
)
@click.option(
    "-r",
    "--reference_name",
    type=click.Choice(REFERENCE_COLLECTION),
    help="Override the reference genome to be used in real-time analysis.",
)
@click.command(short_help="Help setting up nomadic.")
def start(organism, workspace_path, reference_name) -> None:
    """
    Get started with nomadic.

    This command will help you set up a new workspace for nomadic for a specific organism.
    """

    click.echo(f"Workspace will be created at: {workspace_path}")
    init_workspace(workspace_path)
    workspace = Workspace(workspace_path)

    if organism == Organism.pfalciparum:
        setup_pfalciparum(workspace, reference_name=reference_name)

        click.echo(
            "You can now enter your workspace with `cd nomadic` and run `nomadic realtime <experiment_name>` to start real-time analysis."
        )


def setup_pfalciparum(workspace, *, reference_name):
    click.echo("Setting up workspace for Plasmodium falciparum.")

    if reference_name is None:
        reference_name = "Pf3D7"

    # workaround for that we can not set the root in download_reference
    with chdir(workspace.path):
        download_reference(reference_name)
        pass

    copy_bed_files(workspace, Organism.pfalciparum.name)

    copy_example_metadata(workspace)

    default_bed = "nomadsMVP"
    call = True

    click.echo(f"Setting reference genome: {reference_name}")
    click.echo(f"Setting default BED file: {default_bed}")
    click.echo(f"Setting default variant calling: {call}")

    defaults = {
        "reference_name": reference_name,
        "region_bed": default_bed,
        "call": call,
    }

    write_config(
        {"defaults": defaults}, os.path.join(workspace.path, default_config_path)
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


def copy_bed_files(workspace: Workspace, organism_name):
    bed_files = files("nomadic.start").joinpath("data", "beds", organism_name).iterdir()
    click.echo("Copying bed files")
    for bed_file in bed_files:
        data = bed_file.read_text()
        dest_path = os.path.join(workspace.get_beds_dir(), bed_file.name)
        with open(dest_path, "w") as text_file:
            text_file.write(data)
