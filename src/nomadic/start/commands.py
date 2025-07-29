import enum
import os
from importlib.resources import files

import click
import click.shell_completion

from nomadic.util.config import default_config_path, write_config
from nomadic.util.workspace import Workspace


class Organism(enum.Enum):
    pfalciparum = enum.auto()


# Autocomplete is currently not working for enums, see https://github.com/pallets/click/issues/3015
def complete_organism(ctx: click.Context, param, incomplete):
    """Complete organism names based on the Organism enum."""
    return [
        click.shell_completion.CompletionItem(organism.name)
        for organism in Organism
        if organism.name.casefold().startswith(incomplete.casefold())
    ]


@click.argument(
    "organism",
    type=click.Choice(Organism, case_sensitive=False),
    required=True,
    shell_complete=complete_organism,
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

    Currently supported organisms:

      - Plasmodium falciparum (pfalciparum)
    """

    click.echo(f"Workspace will be created at: {workspace_path}")
    workspace = Workspace.create_from_directory(workspace_path)

    if organism == Organism.pfalciparum:
        setup_pfalciparum(workspace)
    else:
        RuntimeError(
            "Organism is not available."
        )  # I am pretty sure it is impossible to enter this code block

    click.echo(
        f"You can now enter your workspace with `cd {workspace_path}` and run `nomadic realtime <experiment_name>` to start real-time analysis."
    )


def setup_pfalciparum(workspace):
    click.echo("Setting up workspace for Plasmodium falciparum.")

    reference_name = "Pf3D7"

    # To speed up initial run, only import here when needed
    from nomadic.download.main import main as download_reference

    download_reference(reference_name)

    copy_bed_files(workspace, organism_name=Organism.pfalciparum.name)

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


def copy_bed_files(workspace: Workspace, *, organism_name):
    bed_files = files("nomadic.start").joinpath("data", "beds", organism_name).iterdir()
    click.echo("Copying amplicon BED files.")
    for bed_file in bed_files:
        data = bed_file.read_text()
        dest_path = os.path.join(workspace.get_beds_dir(), bed_file.name)
        with open(dest_path, "w") as text_file:
            text_file.write(data)
