from itertools import chain
import os
from pathlib import Path
from typing import Optional

import click

from nomadic.util.config import get_config_value, load_config, default_config_path
from nomadic.util.dirs import produce_dir


def find_workspace_root(path: Path) -> Optional[Path]:
    for parent in chain([path], path.resolve().parents):
        if check_if_workspace_root(parent):
            return parent
    return None


def check_if_workspace_root(path: Path) -> bool:
    """
    Check if the given path is a valid Nomadic workspace.
    A valid workspace contains a 'results', 'beds' and 'metadata' directory.
    """
    return (
        path.is_dir()
        and (path / "results").is_dir()
        and (path / "beds").is_dir()
        and (path / "metadata").is_dir()
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

    def get_metadata_xlsx(self, experiment_name: str):
        """
        Get the path to the metadata XLSX file for a given experiment.
        """
        return os.path.join(self.get_metadata_dir(), f"{experiment_name}.xlsx")

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

    def get_config_path(self):
        return os.path.join(self.path, default_config_path)

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

    def get_shared_folder(self) -> str | None:
        """
        Get the shared folder path from the workspace configuration, if it exists.
        """
        config_path = self.get_config_path()
        if os.path.exists(config_path) and os.path.isfile(config_path):
            shared_folder = get_config_value(
                load_config(config_path), ["share", "defaults", "target_dir"]
            )
        else:
            shared_folder = None
        if (
            shared_folder is not None
            and isinstance(shared_folder, str)
            and os.path.isdir(shared_folder)
        ):
            return shared_folder

        return None
