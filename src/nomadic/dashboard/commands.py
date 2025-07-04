import click


@click.command(short_help="Just run the dashboard.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=True,
    help="Name of previous experiment, e.g. '2023-05-12_exptA'.",
)
def dashboard(expt_name):
    """
    Launch the dashboard without performing real-time analysis,
    used to view results of a previous experiment.

    """

    from .main import main

    main(expt_name)
