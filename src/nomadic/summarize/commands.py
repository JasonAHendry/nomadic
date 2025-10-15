from pathlib import Path

import click
import dash

from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(
    short_help="Summarize a set of experiments.",
)
@click.argument(
    "experiment_dirs",
    type=click.Path(exists=True),
    nargs=-1,  # allow multiple arguments; gets passed as tuple
)
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path of the workspace where all input/output files (beds, metadata, results) are stored. "
    "The workspace directory simplifies the use of nomadic in that many arguments don't need to be listed "
    "as they are predefined in the workspace config or can be loaded from the workspace",
)
@click.option(
    "-m",
    "--metadata_csv",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    help="Path to the master metadata CSV file.",
    show_default="<workspace>/metadata/<summary_name>.csv",
)
@click.option("-n", "--summary_name", type=str, help="Name of summary")
@click.option(
    "--dashboard/--no-dashboard",
    default=True,
    help="Whether to start the web dashboard to monitor the run.",
)
def summarize(
    experiment_dirs: tuple[str],
    summary_name: str,
    workspace_path: str,
    metadata_csv: Path,
    dashboard: bool,
):
    """
    Summarize a set of experiments to evaluate quality control and
    mutation prevalence

    """
    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"'{workspace_path}' is not a workspace.",
        )
    workspace = Workspace(workspace_path)

    if summary_name is None:
        summary_name = workspace.get_name()

    if metadata_csv is None:
        metadata_csv = Path(workspace.get_master_metadata_csv(summary_name))

    if not metadata_csv.exists():
        raise click.BadParameter(
            param_hint="-m/--metadata_csv",
            message=f"Master metadata file '{metadata_csv}' does not exist.",
        )

    if len(experiment_dirs) == 0:
        experiment_dirs = workspace.get_experiment_dirs()

    from .main import main

    main(
        expt_dirs=experiment_dirs,
        summary_name=summary_name,
        meta_data_path=metadata_csv,
        show_dashboard=dashboard,
    )
