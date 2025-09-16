import click


@click.command(
    short_help="Summarize a set of experiments.",
)
@click.argument(
    "experiment_dirs",
    type=click.Path(exists=True),
    nargs=-1,  # allow multiple arguments; gets passed as tuple
)
@click.option("-n", "--summary_name", type=str, default="", help="Name of summary")
def summarize(experiment_dirs: tuple[str], summary_name: str):
    """
    Summarize a set of experiments to evaluate quality control and
    mutation prevalence

    """

    from .main import main

    main(expt_dirs=experiment_dirs, summary_name=summary_name)
