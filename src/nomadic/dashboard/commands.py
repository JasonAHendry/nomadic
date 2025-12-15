import difflib
import os
from typing import Optional

import click
import click.shell_completion

from nomadic.util.cli import WORKSPACE_OPTION_KEY, workspace_option
from nomadic.util.exceptions import MetadataFormatError
from nomadic.util.workspace import Workspace


def complete_experiment_name(ctx: click.Context, param, incomplete):
    """Complete experiment names based on existing folders in results directory."""
    workspace: Optional[Workspace] = ctx.params.get(WORKSPACE_OPTION_KEY, None)
    results_dir = workspace.get_results_dir() if workspace is not None else None
    result = []
    if results_dir is not None and os.path.exists(results_dir):
        experiments = [
            dir
            for dir in os.listdir(results_dir)
            if os.path.isdir(os.path.join(results_dir, dir))
        ]
        result = [
            click.shell_completion.CompletionItem(experiment)
            for experiment in experiments
            if experiment.startswith(incomplete)
        ]
    if not result:
        # If no experiments found, return the incomplete input as a directory suggestion as we also accept paths
        return [click.shell_completion.CompletionItem(incomplete, type="dir")]
    return result


def verify_experiment_exists(workspace: Workspace, expt_name: str) -> None:
    """
    When running just the dashboard, we want to verify the expected experiment exists
    Otherwise, we raise exceptions

    """

    expts_found = os.listdir(workspace.get_results_dir())
    if not expts_found:
        raise click.BadParameter("No experiments found in 'results' directory.")

    if expt_name not in expts_found:
        exception_str = (
            f"Found {len(expts_found)} experiments, but none matched '{expt_name}'."
        )
        close_match = difflib.get_close_matches(expt_name, expts_found, cutoff=0.8, n=1)
        if close_match:
            exception_str += f" Did you mean: '{close_match[0]}'?"
        raise click.BadParameter(exception_str)


@click.command(short_help="Just run the dashboard.")
@workspace_option(optional=True)
@click.argument(
    "experiment",
    type=str,
    shell_complete=complete_experiment_name,
)
def dashboard(workspace: Optional[Workspace], experiment):
    """
    Launch the dashboard without performing real-time analysis,
    used to view results of a previous experiment.

    EXPERIMENT can be the name of an previous experiment, located in <workspace>/results/<experiment_name>,
    or a path to a directory containing the results of a previous experiment.

    """

    if workspace is not None:
        input_dir = workspace.get_output_dir(experiment)
    else:
        # Normalize the path to for example remove trailing slashes
        input_dir = os.path.normpath(experiment)
        if "/" not in experiment and not os.path.exists(input_dir):
            # Probably tried to find experiment by name, so warn that workspace does not exist
            raise click.BadParameter(
                param_hint="-w/--workspace",
                message="The current directory is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
            )

    if not os.path.exists(input_dir):
        # Normalize the path to for example remove trailing slashes
        input_dir = os.path.normpath(experiment)
        if not os.path.exists(input_dir):
            if "/" in experiment or workspace is None:
                # Probably tried to find experiment by path, so warn that input directory does not exist
                raise click.BadParameter(
                    param_hint="experiment_name",
                    message=f"Input folder '{input_dir}' does not exist",
                )
            else:
                # Probably tried to find experiment by name, so warn that experiment does not exist
                verify_experiment_exists(workspace, experiment)

    from .main import main

    try:
        main(input_dir)
    except MetadataFormatError as e:
        raise click.ClickException(
            message=str(e),
        ) from e
    except FileNotFoundError as e:
        raise click.ClickException(
            message=str(e),
        ) from e
