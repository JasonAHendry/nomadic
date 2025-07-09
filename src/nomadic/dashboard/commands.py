import os
import difflib

import click

from nomadic.util.workspace import Workspace, check_if_workspace


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
@click.argument(
    "experiment_name",
    type=str,
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the the workspace where all the files will be stored",
)
def dashboard(workspace_path, experiment_name):
    """
    Launch the dashboard without performing real-time analysis,
    used to view results of a previous experiment.

    EXPERIMENT_NAME can be the name of an previous experiment, located in <workspace>/results/<experiment_name>,
    or a path to a directory containing the results of an experiment.

    """

    if check_if_workspace(workspace_path):
        workspace = Workspace(workspace_path)

        input_dir = workspace.get_output_dir(experiment_name)
    else:
        input_dir = experiment_name
        if "/" not in experiment_name and not os.path.exists(input_dir):
            # Probably tried to find experiment by name, so warn that workspace does not exist
            raise click.BadParameter(
                param_hint="-w/--workspace",
                message=f"Workspace '{workspace_path}' does not exist or is not a workspace. Please use nomadic start to create a new workspace, or navigate to your workspace",
            )

    if not os.path.exists(input_dir):
        input_dir = experiment_name
        if not os.path.exists(input_dir):
            if "/" in experiment_name:
                # Probably tried to find experiment by path, so warn that input directory does not exist
                raise click.BadParameter(
                    param_hint="experiment_name",
                    message=f"Input folder '{input_dir}' does not exist",
                )
            else:
                # Probably tried to find experiment by name, so warn that experiment does not exist
                verify_experiment_exists(workspace, experiment_name)

    from .main import main

    main(input_dir)
