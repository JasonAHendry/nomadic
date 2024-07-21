import click
from nomadic.realtime.commands import realtime
from nomadic.download.commands import download
from nomadic.dashboard.commands import dashboard


@click.group()
def cli():
    """
    Run NOMADIC in real-time

    """
    pass


cli.add_command(download)
cli.add_command(realtime)
cli.add_command(dashboard)
