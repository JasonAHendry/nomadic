import os

from nomadic.util.dirs import produce_dir


def init_workspace(path: str):
    if os.path.exists(path):
        raise RuntimeError(f"folder/file {path} already exists")

    produce_dir(path, "results")
    produce_dir(path, "beds")
    produce_dir(path, "metadata")


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
    def __init__(self, path: str):
        self.path = path

    def __str__(self):
        return f"Nomadic Workspace at {self.path}"

    def __repr__(self):
        return f"Workspace(path={self.path})"

    def get_results_dir(self):
        """
        Get the results directory of the workspace.
        """
        return os.path.join(self.path, "results")

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

    def get_meta_data_csv(self, experiment_name: str):
        """
        Get the path to the metadata CSV file for a given experiment.
        """
        return os.path.join(self.get_metadata_dir(), f"{experiment_name}.csv")

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
