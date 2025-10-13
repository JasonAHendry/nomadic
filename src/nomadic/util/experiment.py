import os
import shutil
from pathlib import Path
from typing import NamedTuple

from nomadic.util.dirs import produce_dir
from nomadic.util.metadata import MetadataTableParser
from nomadic.util.regions import RegionBEDParser


class SummaryFiles(NamedTuple):
    """
    Named tuple to hold paths to summary files.
    """

    fastqs_processed: str
    read_mapping: str
    region_coverage: str
    depth_profiles: str
    variants: str


# Currently used summary file names
default_config_path = "config/defaults.json"
summary_files = SummaryFiles(
    fastqs_processed="summary.fastqs_processed.csv",
    read_mapping="summary.read_mapping.csv",
    region_coverage="summary.region_coverage.csv",
    depth_profiles="summary.depth_profiles.csv",
    variants="summary.variants.csv",
)

# Legacy summary file names for backward compatibility
legacy_summary_files = SummaryFiles(
    fastqs_processed="summary.fastq.csv",
    read_mapping="summary.bam_flagstats.csv",
    region_coverage="summary.bedcov.csv",
    depth_profiles="summary.depth.csv",
    variants="summary.variants.csv",
)


class ExperimentDirectories:
    """
    Put all the information about experimental
    directory structure in a single place

    Advantage is any future changes to how the directories
    are organised can be managed here

    Best practice is that objects/functions *do not* directly depend on
    this object, better to use its members as arguments to other
    objects or functions

    TODO:
    - For specific analysis dirs, can have two dictionaries that grow
      - They are forced to be in particular relation to the existing dirs
      - But can be added to or reduced depending on pipeline

    """

    def __init__(
        self,
        output_dir: str,
        metadata: MetadataTableParser,
        regions: RegionBEDParser = None,
        approach_name: str = "",
    ):
        """
        Initialise all the required directories

        """

        self.expt_dir = produce_dir(output_dir)

        self.metadata_dir = produce_dir(self.expt_dir, "metadata")
        self._setup_metadata_dir(metadata, regions)

        # This enables different Guppy versions, barcoding strategies, &c
        self.approach_name = approach_name
        self.approach_dir = produce_dir(self.expt_dir, approach_name)

        self.barcodes_dir = produce_dir(self.approach_dir, "barcodes")
        self._barcode_dirs = {
            b: produce_dir(self.barcodes_dir, b) for b in metadata.barcodes
        }

    def get_barcode_dir(self, barcode_name: str):
        """
        Get the path to a particular `barcode_name`

        """

        return self._barcode_dirs[barcode_name]

    def get_settings_file(self) -> str:
        """Get the path to the setting file for the experiment"""
        return os.path.join(self.metadata_dir, "settings.json")

    def get_summary_files(self) -> SummaryFiles:
        return get_summary_files(Path(self.approach_dir))

    def _setup_metadata_dir(
        self, metadata: MetadataTableParser, regions: RegionBEDParser
    ) -> None:
        """
        Move metadata CSV and regions BED into the metadata directory,
        and store their paths as attributes
        """
        if metadata is not None:
            self.metadata_csv = f"{self.metadata_dir}/{os.path.basename(metadata.csv)}"
            if not os.path.exists(self.metadata_csv):
                metadata.df.to_csv(self.metadata_csv, index=False)

        if regions is not None:
            self.regions_bed = f"{self.metadata_dir}/{os.path.basename(regions.path)}"
            if not os.path.exists(self.regions_bed):
                shutil.copy(regions.path, self.regions_bed)


def get_summary_files(exp_path: Path) -> SummaryFiles:
    if not exp_path.exists():
        raise FileNotFoundError(f"Experiment path does not exist: {exp_path}")
    if (exp_path / legacy_summary_files.fastqs_processed).exists():
        # Use legacy summary files if the old format exists
        return SummaryFiles(
            fastqs_processed=str(exp_path / legacy_summary_files.fastqs_processed),
            read_mapping=str(exp_path / legacy_summary_files.read_mapping),
            region_coverage=str(exp_path / legacy_summary_files.region_coverage),
            depth_profiles=str(exp_path / legacy_summary_files.depth_profiles),
            variants=str(exp_path / legacy_summary_files.variants),
        )
    else:
        return SummaryFiles(
            fastqs_processed=str(exp_path / summary_files.fastqs_processed),
            read_mapping=str(exp_path / summary_files.read_mapping),
            region_coverage=str(exp_path / summary_files.region_coverage),
            depth_profiles=str(exp_path / summary_files.depth_profiles),
            variants=str(exp_path / summary_files.variants),
        )
