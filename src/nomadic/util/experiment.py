import glob
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


def check_complete_experiment(expt_dir: str) -> None:
    """
    Check if an experiment is complete; in reality, it would be nice, at this point, to load an object
    that represents all the files I'd want to work with, e.g. the experiment directories class
    """

    if not os.path.isdir(expt_dir):
        raise FileNotFoundError(f"Experiment directory {expt_dir} does not exist.")

    # We can use this for now, but of course this is getting messy
    _ = find_metadata(expt_dir)
    _ = find_regions(expt_dir)

    used_summary_files = None
    for file_format in [summary_files, legacy_summary_files]:
        if not os.path.exists(f"{expt_dir}/{file_format.fastqs_processed}"):
            continue

        used_summary_files = file_format
        for file in used_summary_files:
            if not os.path.exists(f"{expt_dir}/{file}"):
                raise FileNotFoundError(f"Missing '{file}' file in {expt_dir}.")

    if not used_summary_files:
        raise FileNotFoundError(f"Could not find any summary files in {expt_dir}.")

    # TODO: for now, using this for VCF
    if not os.path.exists(f"{expt_dir}/vcfs"):
        raise FileNotFoundError(f"Could not find VCF directory in {expt_dir}.")


def find_metadata(input_dir: str) -> MetadataTableParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location

    TODO this function is probably not needed anymore, could combine it with get_metadata_csv
    """

    csv = get_metadata_csv(expt_dir=input_dir)
    return MetadataTableParser(csv)


def get_metadata_csv(expt_dir: str) -> str:
    """
    Get the metadata CSV file
    """
    # In most cases, should match experiment name
    metadata_csv = f"{expt_dir}/metadata/{os.path.basename(expt_dir)}.csv"
    if os.path.exists(metadata_csv):
        return metadata_csv
    metadata_csv = glob.glob(f"{expt_dir}/metadata/*.csv")
    if len(metadata_csv) == 1:
        return metadata_csv[0]
    raise ValueError(
        f"Found {len(metadata_csv)} *.csv files in '{expt_dir}/metadata', cannot determine which is metadata."
    )


def find_regions(input_dir: str) -> RegionBEDParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location

    TODO: Bad duplication from above, can write inner function
    """

    metadata_dir = os.path.join(input_dir, "metadata")
    beds = [
        f"{metadata_dir}/{file}"
        for file in os.listdir(metadata_dir)
        if file.endswith(".bed") and not file.endswith(".lowcomplexity_mask.bed")
    ]  # TODO: what about no-suffix files?

    if len(beds) != 1:  # Could alternatively load and LOOK
        raise FileNotFoundError(
            f"Expected one region BED file (*.bed) at {metadata_dir}, but found {len(beds)}."
        )

    return RegionBEDParser(beds[0])
