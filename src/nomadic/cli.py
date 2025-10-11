from collections import OrderedDict

import click

from nomadic.backup.commands import backup
from nomadic.dashboard.commands import dashboard
from nomadic.download.commands import download
from nomadic.process.commands import process
from nomadic.realtime.commands import realtime
from nomadic.start.commands import start


# From: https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
class OrderedGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super(OrderedGroup, self).__init__(name, commands, **attrs)
        #: the registered subcommands by their exported names.
        self.commands = commands or OrderedDict()

    def list_commands(self, ctx):
        return self.commands


@click.group(cls=OrderedGroup)
@click.version_option(message="%(prog)s-v%(version)s")
def cli():
    """
    Mobile sequencing and analysis in real-time

    """
    pass


cli.add_command(start)
cli.add_command(download)
cli.add_command(realtime)
cli.add_command(process)
cli.add_command(dashboard)
cli.add_command(backup)
