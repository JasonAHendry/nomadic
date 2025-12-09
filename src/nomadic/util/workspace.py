import os
import click
from pathlib import Path

from nomadic.util.dirs import produce_dir


def check_if_workspace(path: str) -> bool:
    """
    Check if the given path is a valid Nomadic workspace.
    A valid workspace contains a 'results', 'beds' and 'metadata' directory.
    """
    return (
        os.path.exists(path)
        and os.path.isdir(path)
        and os.path.exists(os.path.join(path, "results"))
        and os.path.exists(os.path.join(path, "beds"))
        and os.path.exists(os.path.join(path, "metadata"))
    )


def looks_like_a_bed_filepath(path: str) -> bool:
    """
    Check if the given path looks like a file path to a bed file.
    A file path typically contains a file extension or a directory structure.
    """
    return ".bed" in path or "/" in path


class Workspace:
    def __init__(self, workspace_path: str):
        self.path = workspace_path

    @classmethod
    def create_from_directory(cls, workspace_path: str):
        """
        Create a new workspace from a directory path
        """
        if os.path.exists(workspace_path):
            raise click.ClickException(
                f"A workspace with the name '{workspace_path}' already exists!"
            )

        produce_dir(
            workspace_path, "results"
        )  # NB: this creates intermediate directories
        produce_dir(workspace_path, "beds")
        produce_dir(workspace_path, "metadata")

        return cls(workspace_path)

    def __str__(self):
        return f"Nomadic Workspace at: {self.path}"

    def __repr__(self):
        return f"Workspace(path={self.path})"

    def get_name(self):
        """
        Get the name of the workspace (the last part of the path).
        """
        return Path(self.path).resolve().name

    def get_results_dir(self):
        """
        Get the results directory of the workspace.
        """
        return os.path.join(self.path, "results")

    def get_output_dir(self, experiment_name: str):
        """
        Get the output directory for a given experiment name.
        """
        return os.path.join(self.get_results_dir(), experiment_name)

    def get_summary_dir(self, summary_name: str):
        """
        Get the summary directory of the workspace.
        """
        return os.path.join(self.path, "summaries", summary_name)

    def get_beds_dir(self):
        """
        Get the beds directory of the workspace.
        """
        return os.path.join(self.path, "beds")

    def get_metadata_dir(self):
        """
        Get the metadata directory of the workspace.
        """
        return os.path.join(self.path, "metadata")

    def get_metadata_csv(self, experiment_name: str):
        """
        Get the path to the metadata CSV file for a given experiment.
        """
        return os.path.join(self.get_metadata_dir(), f"{experiment_name}.csv")

    def get_master_metadata_csv(self, summary_name: str):
        """
        Get the path to the master metadata CSV file for summaries.
        """
        return os.path.join(self.get_metadata_dir(), f"{summary_name}.csv")

    def get_summary_settings_file(self, summary_name: str):
        """
        Get the path to the master metadata CSV file for summaries.
        """
        return os.path.join(self.get_metadata_dir(), f"{summary_name}.yaml")

    def get_bed_file(self, panel_name: str):
        """
        Get the path to the BED file for a given panel name.
        """
        return os.path.join(self.get_beds_dir(), f"{panel_name}.amplicons.bed")

    def get_panel_names(self):
        """
        Get a list of available panel names in the workspace.
        """
        return [
            file.removesuffix(".amplicons.bed")
            for file in os.listdir(self.get_beds_dir())
            if file.endswith(".bed")
        ]

    def get_experiment_names(self):
        """
        Get a list of available experiment names in the workspace.
        """
        if not os.path.exists(self.get_results_dir()):
            return []

        return [
            name
            for name in os.listdir(self.get_results_dir())
            if os.path.isdir(os.path.join(self.get_results_dir(), name))
        ]

    def get_experiment_dirs(self):
        """
        Get a list of available experiment directories in the workspace.
        """
        if not os.path.exists(self.get_results_dir()):
            return []

        return [
            os.path.join(self.get_results_dir(), name)
            for name in self.get_experiment_names()
        ]
