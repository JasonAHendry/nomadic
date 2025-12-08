import glob
import os
import json
import shutil
from pathlib import Path
from typing import NamedTuple

from nomadic.util.dirs import produce_dir
from nomadic.util.metadata import MetadataTableParser, ExtendedMetadataTableParser
from nomadic.util.regions import RegionBEDParser


# --------------------------------------------------------------------------------
# Handle summary file names: legacy and current
#
# --------------------------------------------------------------------------------


class SummaryFiles(NamedTuple):
    """Define summary file names / paths"""

    fastqs_processed: str
    read_mapping: str
    region_coverage: str
    depth_profiles: str
    variants: str


# Currently used summary file names
DEFAULT_CONFIG_PATH = "config/defaults.json"
SUMMARY_NAMES = SummaryFiles(
    fastqs_processed="summary.fastqs_processed.csv",
    read_mapping="summary.read_mapping.csv",
    region_coverage="summary.region_coverage.csv",
    depth_profiles="summary.depth_profiles.csv",
    variants="summary.variants.csv",
)

# Legacy summary file names for backward compatibility
SUMMARY_NAMES_LEGACY = SummaryFiles(
    fastqs_processed="summary.fastq.csv",
    read_mapping="summary.bam_flagstats.csv",
    region_coverage="summary.bedcov.csv",
    depth_profiles="summary.depth.csv",
    variants="summary.variants.csv",
)


def get_summary_files(expt_dir: Path) -> SummaryFiles:
    """
    Determine whether the summary files are use the legacy or current names,
    and return a SummaryFiles object with the appropriate file names

    """

    if not expt_dir.exists():
        raise FileNotFoundError(f"Experiment path does not exist: {expt_dir}")

    if (expt_dir / SUMMARY_NAMES_LEGACY.read_mapping).exists():
        # Detect legacy format using *one* of the differentiating file names
        format_used = SUMMARY_NAMES_LEGACY
    else:
        format_used = SUMMARY_NAMES
    return SummaryFiles(*[str(expt_dir / field) for field in format_used])


# --------------------------------------------------------------------------------
# Define experiment directories
#
# --------------------------------------------------------------------------------


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


# --------------------------------------------------------------------------------
# Checks on experiment outputs
#
# --------------------------------------------------------------------------------


def find_metadata(
    expt_dir: str, Parser: MetadataTableParser = MetadataTableParser
) -> MetadataTableParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location and load it

    """

    # In most cases, should match experiment name
    csv = f"{expt_dir}/metadata/{os.path.basename(expt_dir)}.csv"
    if os.path.exists(csv):
        return Parser(csv)

    csv = glob.glob(f"{expt_dir}/metadata/*.csv")
    if len(csv) == 1:
        return Parser(csv[0])

    raise ValueError(
        f"Found {len(csv)} *.csv files in '{expt_dir}/metadata', cannot determine which is metadata."
    )


def find_regions(expt_dir: str) -> RegionBEDParser:
    """
    Given an experiment directory, search for the metadata CSV file in thee
    expected location

    """

    bed = [
        f
        for f in glob.glob(f"{expt_dir}/metadata/*.bed")
        if f.endswith(".bed") and not f.endswith(".lowcomplexity_mask.bed")
    ]

    if len(bed) == 1:
        return RegionBEDParser(bed[0])

    raise FileNotFoundError(
        f"Expected one region BED file (*.bed) at '{expt_dir}/metadata', but found {len(bed)}."
    )


class ExperimentOutputChecker:
    """Check that all the relevant output files of an experiment present"""

    def __init__(self, expt_dir: str):
        # Confirm experiment directory exists
        if not os.path.isdir(expt_dir):
            raise FileNotFoundError(f"Experiment directory {expt_dir} does not exist.")
        self.expt_dir = expt_dir

        # Load metadata
        parser = find_metadata(expt_dir, Parser=ExtendedMetadataTableParser)
        self.metadata = parser.df
        self.metadata.insert(0, "expt_name", os.path.basename(self.expt_dir))

        # Load regions
        self.regions = find_regions(expt_dir)

        # Other file checks
        self._check_summary_files()
        self._check_vcf_files()
        self._load_settings()
        self._set_variant_caller()

    def _check_summary_files(self):
        # Define the summary file paths
        self.summary_files = get_summary_files(Path(self.expt_dir))

        # Check if they are the legacy files, pretty ugly
        self.summary_files_legacy = self.summary_files.read_mapping.endswith(
            "bam_flagstats.csv"
        )

        # Check they exist
        for file in self.summary_files:
            if "depth" in file:
                # depth files are optional, TODO: not so robust
                continue
            if "fastq" in file:
                # fastq files are optional
                continue
            if not os.path.exists(file):
                raise FileNotFoundError(f"Missing '{file}' file in {self.expt_dir}.")

    def _check_vcf_files(self):
        vcf_dir = f"{self.expt_dir}/vcfs"
        self.vcf_dir_exists = os.path.exists(vcf_dir)
        self.vcf_complete_exists = os.path.exists(f"{vcf_dir}/summary.variants.vcf.gz")
        self.vcf_filtered_exists = os.path.exists(
            f"{vcf_dir}/summary.variants.filtered.annotated.vcf.gz"
        )

    def _load_settings(self):
        settings_path = f"{self.expt_dir}/metadata/settings.json"
        if not os.path.exists(settings_path):
            self.settings = None
        else:
            self.settings = json.load(open(settings_path, "r"))

    def _set_variant_caller(self):
        if self.settings is None:
            self.caller = "bcftools"
        else:
            self.caller = self.settings["caller"]
