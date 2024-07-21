import click


@click.command(short_help="Just run the dashboard.")
@click.option(
    "-e",
    "--expt_name",
    type=str,
    required=False,
    default="0000-00-00_example",
    help="Name of the experiment, used as output directory name. E.g. '2023-05-12_exptA'.",
)
def dashboard(expt_name):
    """
    Just launch the dashboard without performing real-time analysis

    Useful to view results of a previously completed experiment.

    """

    from .main import main

    main(expt_name)
