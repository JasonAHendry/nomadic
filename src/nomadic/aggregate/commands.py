from pathlib import Path

import click
import pandas as pd

from nomadic.util.dirs import produce_dir
from nomadic.util.experiment import legacy_summary_files, summary_files
from nomadic.util.workspace import Workspace, check_if_workspace


@click.command(short_help="Aggregate nomadic summary data from a workspace")
@click.option(
    "-w",
    "--workspace",
    "workspace_path",
    default="./",
    show_default="current directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Path to the nomadic workspace you want to summarise data for",
)
def aggregate(workspace_path: Path):
    """
    Aggregate nomadic summary data for a workspace.
    """

    if not check_if_workspace(workspace_path):
        raise click.BadParameter(
            param_hint="-w/--workspace",
            message=f"'{workspace_path.resolve()}' is not a workspace.",
        )

    click.echo(
        f"Aggregating summary data from nomadic workspace: {workspace_path.resolve().name}"
    )
    workspace = Workspace(str(workspace_path))

    # Identify current and legacy filenames
    sumfiles_current = summary_files._asdict()
    sumfiles_legacy = legacy_summary_files._asdict()
    sumfiles_exclude = ["fastqs_processed", "depth_profiles"]
    sumfiles = {
        key: [sumfiles_current[key], sumfiles_legacy[key]]
        for key in sumfiles_current
        if key not in sumfiles_exclude
    }

    aggregate_dir = workspace_path.resolve() / "aggregate"
    produce_dir(aggregate_dir)

    for filetype, summary in sumfiles.items():
        df = aggregate_summary_file(workspace, summary)
        if df.shape[0] == 0:
            click.echo(click.style("   x", fg="red"), nl=False)
            click.echo(click.style(f"   {filetype} (no data)", fg="red"))
            continue

        click.echo(click.style("   ✓", fg="green"), nl=False)
        click.echo(click.style(f"   {filetype} (aggregated)", fg="green"))
        df.to_csv(aggregate_dir / f"all_{filetype}.csv", index=False)

    click.echo(
        f"All data aggregated to {aggregate_dir.parent.name}/{aggregate_dir.name}"
    )


def aggregate_summary_file(
    workspace: Workspace, summary_fns: list[str]
) -> pd.DataFrame:
    """
    Extract and aggregate a summary file into a df with expid as a new col
    """
    results_dir = Path(workspace.get_results_dir())
    aggregate_df = pd.DataFrame()

    for summary_fn in set(summary_fns):
        for file in results_dir.glob(f"**/{summary_fn}"):
            temp_df = pd.read_csv(file)
            temp_df["exp_id"] = file.parent.name
            aggregate_df = pd.concat([aggregate_df, temp_df])

    return aggregate_df
